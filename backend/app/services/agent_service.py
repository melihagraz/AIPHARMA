"""LLM-driven chat agent (Phase 2).

A thin orchestration layer over the existing production pipeline. Claude
sees a small, targeted toolset and decides when to call which tool — the
backend wraps these tools as plain async functions over the existing
service clients (REINVENT 4 / ADMET-AI / Vina / ChEMBL / pipeline).

Why a manual loop instead of the SDK's tool runner: we need to stream
intermediate tool_use / tool_result events to the SSE channel as they
happen so the UI can show "🛠️ generating molecules…" live. The SDK
runner is great when you only want the final answer, but here progress
visibility is the whole UX.

Public surface:
    - ``stream_chat(request) -> AsyncIterator[ChatStreamEvent]``
        the only entrypoint used by routes_agent.

Configuration:
    - ``settings.anthropic_api_key`` (required at request time; empty
      string ⇒ raises a clean error instead of failing inside the SDK)
    - ``settings.agent_model`` (default ``claude-opus-4-7``)
    - ``settings.agent_max_tokens``
    - ``settings.agent_max_iterations`` (loop safety stop)
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from typing import Any, AsyncIterator
from uuid import uuid4

try:
    import anthropic
except ImportError:  # optional LLM feature; import API even if the SDK is absent
    class _AnthropicUnavailable:
        class APIStatusError(Exception):
            pass
        class APIConnectionError(Exception):
            pass
        class AsyncAnthropic:
            def __init__(self, **kwargs):
                raise RuntimeError("Install anthropic to enable the optional chat agent")
    anthropic = _AnthropicUnavailable()

from app.core.config import get_settings
from app.db.session import SessionLocal
from app.models.audit_log import AuditLog
from app.models.job import Job
from app.repositories.job_repository import JobRepository
from app.schemas.agent import (
    ChatDone,
    ChatError,
    ChatRequest,
    ChatStreamEvent,
    ChatTextDelta,
    ChatThinkingDelta,
    ChatToolResult,
    ChatToolUse,
)
from app.schemas.job import JobCreateRequest
from app.services import (
    admet_client,
    chembl_client,
    docking_client,
    pocket_client,
    production_pipeline,
    reinvent_client,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# AuditLog persistence — make agent runs first-class citizens.
#
# Before this hook, agent tool calls were ephemeral: results lived in the chat
# bubble only, vanishing when the session closed. Meanwhile the UI's "Step 3 /
# 4 / 5" buttons hit /api/jobs/start-* endpoints that persisted to the same
# AuditLog the Analysis Results page reads. So a researcher who spent 12
# minutes letting the agent run a 4-round optimisation lost the entire trail
# the moment they refreshed.
#
# These helpers mirror the exact payload shape the UI writer uses (action +
# results_json.data) so the existing ResultsModal branches render agent runs
# without any front-end changes. ``actor_id="agent"`` distinguishes them from
# button-click runs in case we ever want to filter by origin.
#
# Failure mode is "log and shrug" — DB persistence must never break the chat
# response itself. The tool's main payload still goes back to the LLM.
# ---------------------------------------------------------------------------


def _safe_receptor_label(target: str | None, prefix: str) -> str:
    """Build a pdb_id/receptor_name value that respects the 255-char column."""
    base = (target or "AGENT").strip().upper() or "AGENT"
    label = f"{prefix}-{base}-{uuid4().hex[:8]}"
    return label[:255]


def _persist_agent_audit(
    *,
    action: str,
    receptor_label: str,
    candidate_count: int,
    payload: dict,
) -> dict[str, Any]:
    """Write one AuditLog row mirroring the UI's writer shape.

    Returns ``{"job_id": ..., "persisted": True}`` on success or
    ``{"persisted": False, "persist_error": "..."}`` on failure — the caller
    merges this into the tool result so the LLM can mention the saved
    job_id to the user (and so a future agent turn can reference it).
    """
    db = SessionLocal()
    try:
        repo = JobRepository(db)
        job = repo.create(
            job_id=str(uuid4()),
            payload=JobCreateRequest(
                tenant_id="default",
                project_id="agent-chat",
                pdb_id=receptor_label,
                receptor_name=receptor_label,
                mode="agonist",
                candidate_count=max(1, candidate_count),
            ),
        )
        pinned = db.query(Job).filter(Job.id == job.id).first()
        if pinned is not None:
            pinned.status = "completed"
            pinned.updated_at = datetime.utcnow()
        db.add(
            AuditLog(
                id=str(uuid4()),
                actor_id="agent",
                action=action,
                resource_type="job_output",
                resource_id=job.id,
                payload=json.dumps(payload, ensure_ascii=True, default=str),
            )
        )
        db.commit()
        return {"job_id": job.id, "persisted": True}
    except Exception as exc:  # noqa: BLE001
        try:
            db.rollback()
        except Exception:  # noqa: BLE001
            pass
        logger.warning("Agent run persist failed (%s): %s", action, exc)
        return {"persisted": False, "persist_error": str(exc)}
    finally:
        db.close()


# ---------------------------------------------------------------------------
# System prompt — kept stable so prompt-cache hits across turns
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """You are AI Pharma Assistant, the chat layer over a production drug-discovery pipeline.

Your job: take a researcher's free-text request (often Turkish), pick the right tools, run them, then present the result clearly. You have no opinions about a target's clinical merit; you execute.

Available capabilities (via tools):
  • Generate target-aware molecules with REINVENT 4 (uses fine-tuned prior if available, else base + ChEMBL Tanimoto re-ranking).
  • Score ADMET (52 endpoints) with ADMET-AI.
  • Dock against pre-bundled receptors (EGFR, JAK2, BACE1, KDR, MET, AKT1, LCK, KIT) with AutoDock Vina.
  • Run the full production pipeline (REINVENT → ChEMBL → ADMET → Vina) end-to-end.
  • Look up known potent actives in ChEMBL.
  • List which targets currently have a fine-tuned prior or a docking receptor.

Operating rules:
  • Respond in the user's language (Turkish or English). Default to Turkish if the user used Turkish.
  • For multi-step requests prefer ``run_full_pipeline`` over chaining the smaller tools — it is faster and produces ranked candidates with all signals attached.
  • Use the smaller tools when the user explicitly asks for one stage (e.g. "sadece üret, dock etme").
  • Vina docking is expensive (~5–30s per ligand). Never dock more than 30 ligands without confirming.
  • When the user names a target (e.g. "EGFR"), pass it as ``target_name`` — the system resolves ChEMBL + fine-tuned-prior + Vina receptor matches automatically.
  • Number values: present binding affinities in kcal/mol with one decimal, ADMET probabilities as 0.00–1.00 with two decimals.
  • Lead-optimisation loop: when the user (or you) identified a hit in a previous run and the request is "optimise around it" / "analoglarını üret" / "second round" / "focused round" / "decorate scaffold", call ``generate_analogs`` with the seed SMILES. Then volunteer to chain ``predict_admet`` + ``dock_molecules`` on the top analogs without re-asking — that is the canonical second-round workflow.
  • Persistence: ``run_full_pipeline``, ``generate_analogs``, ``predict_admet``, and ``dock_molecules`` now write their results to the same database the UI's Analysis Results page reads. Each tool returns ``saved_as_job_id`` on success — append a short "(saved as <8-char-id>...)" trailing parenthetical to your reply when you cite a result, so the user can find the run in the UI later. Do NOT show the full UUID, just the first 8 characters.

Output format (very important — keep responses tight and scannable):
  • Lead with one short sentence (≤ 15 words) summarising the headline result.
  • Then a 4–6 line bullet list: pipeline used, generated/valid/docked counts, reference score, the take-away. No header above this list.
  • Then a markdown table with the top 3–5 candidates: columns = Rank, Vina (kcal/mol), QED, MW, hERG, Tanimoto, SMILES. ALWAYS print the full SMILES string — never truncate, never use "...". The UI handles wrapping. At most 5 rows.
  • Then a short analysis of 2–4 bullets: which candidate wins, which to drop, what's risky. Be specific about hERG > 0.5 = cardiac risk, QED < 0.3 = poor drug-likeness.
  • End with a single suggested next action sentence (no header). One sentence, not a paragraph.

Language consistency (critical — pick ONE language per response and stick to it):
  • If responding in Turkish: use Turkish labels — table header "Vina" (NOT "VİNA"; use ASCII), section label "Yorum:" if you need one, units in TR.
  • If responding in English: use English labels — table header "Vina", section label "Analysis:" if you need one, units in EN.
  • NEVER mix: do not write a Turkish "Yorum:" label above an English analysis or vice versa.
  • Use plain ASCII for column headers — no Turkish capitalisation quirks (İ, Ş, etc.). Always "Vina", "Tanimoto", not "VİNA", "TAN.".
  • Do not use H1 (#) — H2 (##) only when needed.
  • Do not add success badges (✅) or decorative pins (📌). Skip empty filler lines.

You are not a medical advisor. Do not invent data: if a tool fails or returns empty, say so plainly.
"""


# ---------------------------------------------------------------------------
# Tool schemas (sent to Claude). Names mirror the dispatch table below.
# ---------------------------------------------------------------------------

TOOLS: list[dict[str, Any]] = [
    {
        "name": "generate_molecules",
        "description": (
            "Generate de novo small molecules with REINVENT 4. If a target_name is "
            "provided, automatically uses a fine-tuned prior (if one exists) and "
            "scores each molecule against ChEMBL known actives via ECFP4 Tanimoto. "
            "Returns SMILES + RDKit descriptors (mw, logp, qed, tpsa) + max_tanimoto."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_name": {
                    "type": "string",
                    "description": "Target name like 'EGFR'. Optional — leave empty for generic generation.",
                },
                "batch_size": {
                    "type": "integer",
                    "minimum": 50,
                    "maximum": 1000,
                    "default": 200,
                },
                "seed": {
                    "type": "integer",
                    "description": "Optional RNG seed for reproducibility.",
                },
            },
            "required": [],
        },
    },
    {
        "name": "predict_admet",
        "description": (
            "Predict 52 ADMET endpoints (hERG, DILI, AMES, Caco2, PPB, Bioavailability, "
            "etc.) for a list of SMILES using ADMET-AI."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                    "maxItems": 200,
                },
            },
            "required": ["smiles"],
        },
    },
    {
        "name": "dock_molecules",
        "description": (
            "Dock SMILES against a pre-bundled Vina receptor and return binding "
            "affinity in kcal/mol (more negative = stronger binder)."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                    "maxItems": 30,
                },
                "target": {
                    "type": "string",
                    "description": "One of: EGFR, JAK2, BACE1, KDR, MET, AKT1, LCK, KIT.",
                },
                "exhaustiveness": {"type": "integer", "minimum": 1, "maximum": 32, "default": 8},
            },
            "required": ["smiles", "target"],
        },
    },
    {
        "name": "run_full_pipeline",
        "description": (
            "End-to-end production discovery: REINVENT 4 generation → ChEMBL "
            "Tanimoto → ADMET-AI scoring → Vina docking on top-K → ranked top 20 "
            "candidates with component scores. This is the preferred tool for "
            "complete drug-discovery requests."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_name": {"type": "string"},
                "generate_count": {
                    "type": "integer",
                    "minimum": 50,
                    "maximum": 1000,
                    "default": 200,
                },
                "dock_top_k": {
                    "type": "integer",
                    "minimum": 0,
                    "maximum": 30,
                    "default": 10,
                    "description": "How many ADMET-passing molecules to dock. 0 = skip docking.",
                },
                "seed": {"type": "integer"},
            },
            "required": ["target_name"],
        },
    },
    {
        "name": "generate_analogs",
        "description": (
            "Lead-optimization tool: generate scaffold-similar analogs around a "
            "seed SMILES (typically a winning candidate from a previous "
            "``run_full_pipeline``). Oversamples from REINVENT 4 and keeps only "
            "molecules with ECFP4 Tanimoto ≥ ``min_similarity`` to the seed. "
            "Use this AFTER you have a chemotype hit and want focused variants "
            "to triage by ADMET + docking — i.e. the second round of the design "
            "loop, not the first. Pairs naturally with ``predict_admet`` + "
            "``dock_molecules`` to evaluate the analog series."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "scaffold_smiles": {
                    "type": "string",
                    "description": "Seed SMILES — the chemotype to decorate around (e.g. rank 1 from a prior pipeline run).",
                },
                "n_variants": {
                    "type": "integer",
                    "minimum": 5,
                    "maximum": 100,
                    "default": 30,
                    "description": "How many analogs to return (after similarity filtering).",
                },
                "min_similarity": {
                    "type": "number",
                    "minimum": 0.15,
                    "maximum": 0.9,
                    "default": 0.2,
                    "description": (
                        "Minimum ECFP4 Tanimoto to the seed scaffold (0–1). "
                        "Lower = more diverse analogs; higher = closer cousins. "
                        "Default 0.2 matches REINVENT 4's diversity baseline — "
                        "raising to 0.3+ usually produces zero hits unless the "
                        "seed is a common chemotype already heavily represented "
                        "in the prior."
                    ),
                },
                "target_name": {
                    "type": "string",
                    "description": "Optional target for fine-tuned prior selection (e.g. 'EGFR').",
                },
                "seed": {"type": "integer"},
            },
            "required": ["scaffold_smiles"],
        },
    },
    {
        "name": "fetch_chembl_actives",
        "description": (
            "Pull known potent actives for a target from ChEMBL. Useful for "
            "context ('what does literature look like for this target?')."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target_name": {"type": "string"},
                "pchembl_min": {"type": "number", "default": 7.0, "minimum": 5.0, "maximum": 11.0},
                "limit": {"type": "integer", "default": 20, "minimum": 1, "maximum": 100},
            },
            "required": ["target_name"],
        },
    },
    {
        "name": "list_capabilities",
        "description": (
            "List which targets currently have a fine-tuned REINVENT prior and "
            "which have a Vina docking receptor. No arguments."
        ),
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "detect_pocket",
        "description": (
            "Locate a druggable binding pocket. For one of the 8 bundled targets "
            "(EGFR, JAK2, BACE1, KDR, MET, AKT1, LCK, KIT) returns the validated "
            "pocket coordinates. For any other target supply a PDB ID; the bound "
            "ligand's coordinates are used as the pocket. Returns center XYZ, "
            "bounding box, druggability score, and pocket-lining residues."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {"type": "string", "description": "Bundled target name (fastest)."},
                "pdb_id": {"type": "string", "description": "4-char PDB ID (e.g. '1M17')."},
                "ligand_resname": {
                    "type": "string",
                    "description": "Optional 3-char HETATM ligand resname (e.g. 'AQ4').",
                },
            },
            "required": [],
        },
    },
    {
        "name": "extract_pharmacophore",
        "description": (
            "Extract a 3D pharmacophore from a co-crystal ligand. Returns a list "
            "of features (Donor, Acceptor, Aromatic, Hydrophobe, PosIonizable, "
            "NegIonizable) each with 3D position. Use this to understand the "
            "binding requirements of known inhibitors before generating new ones."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {
                    "type": "string",
                    "description": "Bundled target name; uses a known co-crystal PDB.",
                },
                "pdb_id": {"type": "string", "description": "4-char PDB ID."},
                "ligand_resname": {"type": "string"},
            },
            "required": [],
        },
    },
]


# ---------------------------------------------------------------------------
# Tool dispatch — async, returns plain JSON-serialisable dicts (truncated).
# Each handler MUST never raise: convert errors into ``{"error": "..."}``.
# ---------------------------------------------------------------------------


def _truncate_candidates(candidates: list[dict[str, Any]], keep: int) -> list[dict[str, Any]]:
    """Slim down a candidate list for LLM context. Keeps key fields only."""
    out = []
    for c in candidates[:keep]:
        out.append(
            {
                "rank": c.get("rank"),
                "smiles": c.get("smiles"),
                "score": c.get("score"),
                "dock_score": c.get("dock_score"),
                "descriptors": c.get("descriptors"),
                "component_scores": c.get("component_scores"),
            }
        )
    return out


async def _tool_generate_molecules(args: dict[str, Any]) -> dict[str, Any]:
    target = (args.get("target_name") or "").strip() or None
    batch = int(args.get("batch_size") or 200)
    seed = args.get("seed")
    reference = None
    if target:
        reference = await chembl_client.fetch_known_actives(target, limit=30) or None
    try:
        candidates, prior_used = await reinvent_client.sample_smiles(
            batch_size=batch,
            seed=seed,
            validate_with_rdkit=True,
            reference_smiles=reference,
            target_prior=target,
        )
    except Exception as exc:  # noqa: BLE001
        return {"error": f"REINVENT failed: {exc}"}
    valid = [c for c in candidates if c.get("valid")]
    valid.sort(key=lambda c: (c.get("max_tanimoto") or 0.0, -float(c.get("nll") or 0.0)), reverse=True)
    return {
        "prior_used": prior_used,
        "reference_actives_count": len(reference or []),
        "total_generated": len(candidates),
        "valid_count": len(valid),
        "top_molecules": [
            {
                "smiles": c["smiles"],
                "nll": round(float(c.get("nll", 0.0)), 3),
                "qed": round(float(c.get("qed", 0.0) or 0.0), 3),
                "mw": round(float(c.get("mw", 0.0) or 0.0), 2),
                "logp": round(float(c.get("logp", 0.0) or 0.0), 3),
                "max_tanimoto": (
                    round(float(c["max_tanimoto"]), 3) if c.get("max_tanimoto") is not None else None
                ),
            }
            for c in valid[:20]
        ],
    }


async def _tool_predict_admet(args: dict[str, Any]) -> dict[str, Any]:
    smiles = list(args.get("smiles") or [])
    if not smiles:
        return {"error": "smiles list is empty"}
    try:
        preds = await admet_client.predict_admet(smiles)
    except Exception as exc:  # noqa: BLE001
        return {"error": f"ADMET-AI failed: {exc}"}
    keep_keys = (
        "hERG",
        "DILI",
        "AMES",
        "BBB_Martins",
        "Bioavailability_Ma",
        "Caco2_Wang",
        "PPBR_AZ",
        "Half_Life_Obach",
        "CYP3A4_Veith",
        "CYP2D6_Veith",
    )
    out = []
    for p in preds:
        raw = p.get("raw") or {}
        out.append(
            {
                "smiles": p.get("smiles"),
                **{k: round(float(raw.get(k, 0.0) or 0.0), 4) for k in keep_keys if k in raw},
            }
        )

    # Persist as admet_output so this lands in Analysis Results > ADMET
    # Filtering tab. Mirror the routes_legacy.py /admet/run payload shape so
    # the table renders MW/LogP/TPSA/hERG (computed via RDKit + categorical
    # herg label) without UI changes.
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Crippen  # type: ignore
        _rdkit_ok = True
    except Exception:  # noqa: BLE001
        _rdkit_ok = False

    def _props3(smi: str | None) -> tuple[float | None, float | None, float | None]:
        if not _rdkit_ok or not smi:
            return None, None, None
        try:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                return None, None, None
            return (
                round(Descriptors.MolWt(m), 1),
                round(Crippen.MolLogP(m), 2),
                round(Descriptors.TPSA(m), 1),
            )
        except Exception:  # noqa: BLE001
            return None, None, None

    persisted_rows: list[dict[str, Any]] = []
    for idx, p in enumerate(preds, start=1):
        smi = p.get("smiles")
        raw = p.get("raw") or {}
        pct = p.get("percentile") or {}
        mw, logp, tpsa = _props3(smi)
        herg_raw = raw.get("hERG")
        if isinstance(herg_raw, (int, float)):
            herg_label = "VAR" if float(herg_raw) > 0.5 else "YOK"
        else:
            herg_label = None
        descriptors_block = {
            "hERG": raw.get("hERG"),
            "DILI": raw.get("DILI"),
            "AMES": raw.get("AMES"),
            "BBB": raw.get("BBB_Martins"),
            "Bioavailability": raw.get("Bioavailability_Ma"),
            "Caco2": raw.get("Caco2_Wang"),
            "PPBR": raw.get("PPBR_AZ"),
            "CYP3A4_inhib": raw.get("CYP3A4_Veith"),
            "CYP2D6_inhib": raw.get("CYP2D6_Veith"),
            "CYP2C9_inhib": raw.get("CYP2C9_Veith"),
            "Half_Life_h": raw.get("Half_Life_Obach"),
        }
        persisted_rows.append(
            {
                "rank": idx,
                "smiles": smi,
                "module": "general",
                "score": None,
                "dock_score": None,
                "mw": mw,
                "logp": logp,
                "tpsa": tpsa,
                "herg": herg_label,
                "alerts": [],
                "component_scores": {},
                "descriptors": {
                    **{k: v for k, v in descriptors_block.items() if v is not None},
                    **({"mw": mw} if mw is not None else {}),
                    **({"logp": logp} if logp is not None else {}),
                    **({"tpsa": tpsa} if tpsa is not None else {}),
                },
                "admet_raw": raw,
                "admet_percentile": pct,
            }
        )

    persistence = _persist_agent_audit(
        action="admet_output",
        receptor_label=_safe_receptor_label(None, "AGENT-ADMET"),
        candidate_count=len(persisted_rows),
        payload={
            "analysis_section": "admet",
            "source_kind": "agent",
            "input_count": len(smiles),
            "output_count": len(persisted_rows),
            "results_json": {
                "section": "admet",
                "total_input": len(smiles),
                "total_output": len(persisted_rows),
                "data": persisted_rows,
            },
            "source": "agent_predict_admet",
        },
    )
    return {
        "count": len(out),
        "predictions": out,
        "saved_as_job_id": persistence.get("job_id"),
        "persisted": persistence.get("persisted", False),
    }


async def _tool_dock_molecules(args: dict[str, Any]) -> dict[str, Any]:
    smiles = list(args.get("smiles") or [])
    target = (args.get("target") or "").strip()
    if not smiles or not target:
        return {"error": "smiles and target are both required"}
    available = await docking_client.list_available_targets()
    match = next((t for t in available if t.lower() == target.lower()), None)
    if match is None:
        return {
            "error": f"no Vina receptor for target {target!r}",
            "available_targets": available,
        }
    try:
        results = await docking_client.dock_batch(
            smiles,
            match,
            exhaustiveness=int(args.get("exhaustiveness") or 8),
        )
    except Exception as exc:  # noqa: BLE001
        return {"error": f"docking failed: {exc}"}
    results.sort(key=lambda r: (r.get("affinity") if r.get("affinity") is not None else 0.0))

    # Build per-candidate rows matching what /api/docking/run writes so the
    # Analysis Results > Docking tab renders the same way for agent runs as
    # for UI-button runs. RDKit descriptors (MW, LogP) computed inline so the
    # docking table's MW/LogP columns are not empty.
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors, Crippen  # type: ignore
        _rdkit_ok = True
    except Exception:  # noqa: BLE001
        _rdkit_ok = False

    def _props(smi: str | None) -> tuple[float | None, float | None]:
        if not _rdkit_ok or not smi:
            return None, None
        try:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                return None, None
            return round(Descriptors.MolWt(m), 1), round(Crippen.MolLogP(m), 2)
        except Exception:  # noqa: BLE001
            return None, None

    by_smi: dict[str, float] = {}
    for r in results:
        a = r.get("affinity")
        s = r.get("smiles")
        if isinstance(s, str) and isinstance(a, (int, float)):
            by_smi[s] = float(a)

    persisted_rows: list[dict[str, Any]] = []
    for idx, smi in enumerate(smiles, start=1):
        aff = by_smi.get(smi)
        mw, logp = _props(smi)
        aff_rounded = round(aff, 3) if aff is not None else None
        persisted_rows.append(
            {
                "rank": idx,
                "smiles": smi,
                "module": "general",
                "score": None,
                "dock_score": aff_rounded,
                "affinity": aff_rounded,
                "mw": mw,
                "logp": logp,
                "alerts": [],
                "component_scores": {},
                "descriptors": {
                    "vina_affinity_kcal_mol": aff_rounded,
                    "mw": mw,
                    "logp": logp,
                },
            }
        )
    persisted_rows.sort(key=lambda c: (c["dock_score"] is None, c["dock_score"] or 0))
    succeeded = sum(1 for c in persisted_rows if c["dock_score"] is not None)

    persistence = _persist_agent_audit(
        action="docking_output",
        receptor_label=_safe_receptor_label(match, "AGENT-DOCK"),
        candidate_count=len(persisted_rows),
        payload={
            "analysis_section": "docking",
            "docking_target": match,
            "input_count": len(smiles),
            "docked_count": len(smiles),
            "succeeded": succeeded,
            "results_json": {
                "section": "docking",
                "target": match,
                "input_count": len(smiles),
                "docked_count": len(smiles),
                "succeeded": succeeded,
                "data": persisted_rows,
            },
            "source": "agent_dock_molecules",
        },
    )

    return {
        "target_used": match,
        "count": len(results),
        "results": [
            {
                "smiles": r.get("smiles"),
                "affinity_kcal_mol": (
                    round(float(r["affinity"]), 3) if r.get("affinity") is not None else None
                ),
                "error": r.get("error"),
            }
            for r in results
        ],
        "saved_as_job_id": persistence.get("job_id"),
        "persisted": persistence.get("persisted", False),
    }


async def _tool_run_full_pipeline(args: dict[str, Any]) -> dict[str, Any]:
    target = (args.get("target_name") or "").strip()
    if not target:
        return {"error": "target_name is required"}
    req = {
        "target_name": target,
        "generate_count": int(args.get("generate_count") or 200),
        "dock_top_k": int(args.get("dock_top_k") or 10),
        "enable_docking": int(args.get("dock_top_k") or 10) > 0,
        "use_chembl": True,
        "seed": args.get("seed"),
    }
    try:
        resp = await production_pipeline.run_production_discovery(req)
    except Exception as exc:  # noqa: BLE001
        return {"error": f"pipeline failed: {exc}"}
    candidates = resp.get("candidates") or []

    # Persist to AuditLog so this run shows up in Analysis Results > De Novo
    # tab. Payload shape mirrors what routes_legacy.py's start-denovo writer
    # produces (results_json.data is what ResultsModal iterates), so the
    # existing front-end branches render agent runs with no UI changes.
    persistence = _persist_agent_audit(
        action="denovo_output",
        receptor_label=_safe_receptor_label(target, "AGENT-DENOVO"),
        candidate_count=len(candidates),
        payload={
            "analysis_section": "denovo",
            "target": resp.get("target"),
            "module": resp.get("module"),
            "total_generated": resp.get("total_generated"),
            "after_filter_count": resp.get("after_filter_count"),
            "top_n": resp.get("top_n"),
            "results_json": {
                "target": resp.get("target"),
                "module": resp.get("module"),
                "total_generated": resp.get("total_generated"),
                "after_filter_count": resp.get("after_filter_count"),
                "top_n": resp.get("top_n"),
                "data": candidates,
            },
            "meta": resp.get("meta"),
            "source": "agent_run_full_pipeline",
        },
    )

    # ALSO write a docking_output row when Vina actually ran on the top-K,
    # so the Analysis Results > Docking tab surfaces this run too. Without
    # this side-write the Docking tab's "🤖 Agent (N)" counter sits at 0
    # despite agent having spent ~30s/molecule on real Vina — confusing
    # for the researcher reviewing the trail. We extract only the docked
    # subset (those with dock_score not None) so the docking table doesn't
    # carry empty rows. This is the same payload shape /api/docking/run
    # writes, so the existing ResultsModal isDockingMode branch picks it
    # up without further code changes.
    docking_persistence: dict[str, Any] = {}
    docked_subset = [c for c in candidates if c.get("dock_score") is not None]
    meta = resp.get("meta") or {}
    docking_target = meta.get("docking_target") or resp.get("target")
    if docked_subset and docking_target:
        docking_persistence = _persist_agent_audit(
            action="docking_output",
            receptor_label=_safe_receptor_label(docking_target, "AGENT-DOCK"),
            candidate_count=len(docked_subset),
            payload={
                "analysis_section": "docking",
                "docking_target": docking_target,
                "input_count": len(docked_subset),
                "docked_count": len(docked_subset),
                "succeeded": len(docked_subset),
                "results_json": {
                    "section": "docking",
                    "target": docking_target,
                    "input_count": len(docked_subset),
                    "docked_count": len(docked_subset),
                    "succeeded": len(docked_subset),
                    # Reuse the same candidate rows — they already carry
                    # smiles, dock_score, mw, logp, descriptors, etc.
                    "data": docked_subset,
                },
                "parent_job_id": persistence.get("job_id"),
                "source": "agent_run_full_pipeline_docking_sidewrite",
            },
        )

    return {
        "target": resp.get("target"),
        "total_generated": resp.get("total_generated"),
        "after_filter_count": resp.get("after_filter_count"),
        "top_n": resp.get("top_n"),
        "meta": resp.get("meta"),
        "top_candidates": _truncate_candidates(candidates, keep=10),
        # The LLM can mention this id to the user and we can cite it in
        # later turns ("from job AGENT-DENOVO-JAK2-... we found ...").
        "saved_as_job_id": persistence.get("job_id"),
        "saved_dock_job_id": docking_persistence.get("job_id"),
        "persisted": persistence.get("persisted", False),
    }


async def _tool_generate_analogs(args: dict[str, Any]) -> dict[str, Any]:
    """Generate scaffold-similar analogs around a seed SMILES (lead optimization).

    Strategy
    --------
    1. Validate the seed SMILES with RDKit (canonicalise).
    2. Oversample N×10 (min 200) from the REINVENT 4 prior — base or
       fine-tuned, depending on ``target_name`` resolution.
    3. Have the Reinvent_App compute ECFP4 Tanimoto vs the seed (passing
       the seed as the ``reference_smiles`` of length 1).
    4. Filter on ``min_similarity ≤ Tanimoto < 0.999`` — drop unrelated
       molecules AND drop near-exact duplicates of the seed.
    5. Rank by 0.6·Tanimoto + 0.4·QED so we surface analogs that are both
       similar enough to share SAR with the seed AND drug-like enough to
       progress.

    This is the "focused second round" the agent should reach for once
    ``run_full_pipeline`` has produced a hit chemotype. It does NOT score
    ADMET or dock — that's left to ``predict_admet`` + ``dock_molecules``
    so the agent can decide which analogs are worth the Vina time.
    """
    scaffold = (args.get("scaffold_smiles") or "").strip()
    if not scaffold:
        return {"error": "scaffold_smiles is required"}
    n_variants = int(args.get("n_variants") or 30)
    min_similarity = float(args.get("min_similarity") or 0.2)
    target = (args.get("target_name") or "").strip() or None
    seed = args.get("seed")

    # Validate the seed — fail fast on garbage input so the agent knows
    # the request was malformed rather than getting silent empty results.
    try:
        from rdkit import Chem  # type: ignore
        from rdkit import RDLogger  # type: ignore
        RDLogger.DisableLog("rdApp.*")
        seed_mol = Chem.MolFromSmiles(scaffold)
        if seed_mol is None:
            return {"error": f"scaffold SMILES could not be parsed: {scaffold!r}"}
        canonical_scaffold = Chem.MolToSmiles(seed_mol)
    except Exception as exc:  # noqa: BLE001
        return {"error": f"scaffold parsing failed: {exc}"}

    # Oversample. The Tanimoto-≥-0.3 hit rate against an arbitrary scaffold
    # is typically 5–15% from the base prior, so 10× is a reasonable floor.
    # Cap at 1000 to respect the Reinvent_App batch limit.
    oversample = min(1000, max(200, n_variants * 10))
    try:
        candidates, prior_used = await reinvent_client.sample_smiles(
            batch_size=oversample,
            seed=seed,
            validate_with_rdkit=True,
            reference_smiles=[canonical_scaffold],  # seed = reference for Tanimoto
            target_prior=target,
        )
    except Exception as exc:  # noqa: BLE001
        return {"error": f"REINVENT sampling failed: {exc}"}

    # Filter: valid SMILES + Tanimoto in [min_similarity, 0.999) — the upper
    # bound excludes near-identical regenerations of the seed itself, which
    # the prior occasionally produces and would inflate the "analog count".
    similar: list[dict[str, Any]] = []
    near_duplicate_count = 0
    for c in candidates:
        if not c.get("valid"):
            continue
        tani = c.get("max_tanimoto")
        if tani is None:
            continue
        if tani >= 0.999:
            near_duplicate_count += 1
            continue
        if tani < min_similarity:
            continue
        similar.append(c)

    # Rank by combined similarity + drug-likeness — neither alone is enough.
    # 0.6/0.4 split favours similarity (analog series should stay close to
    # the seed scaffold) while still demoting molecules with poor QED.
    similar.sort(
        key=lambda c: (
            (c.get("max_tanimoto") or 0.0) * 0.6
            + (c.get("qed") or 0.0) * 0.4
        ),
        reverse=True,
    )
    top = similar[:n_variants]

    # Build the rows we both return to the LLM and persist to AuditLog.
    analog_rows = [
        {
            "rank": i + 1,
            "smiles": c["smiles"],
            "module": "general",
            "score": None,
            "dock_score": None,
            "tanimoto_to_seed": round(float(c.get("max_tanimoto", 0.0)), 3),
            "qed": round(float(c.get("qed", 0.0) or 0.0), 3),
            "mw": round(float(c.get("mw", 0.0) or 0.0), 2),
            "logp": round(float(c.get("logp", 0.0) or 0.0), 3),
            "tpsa": round(float(c.get("tpsa", 0.0) or 0.0), 2),
            "nll": round(float(c.get("nll", 0.0)), 3),
            "alerts": [],
            "component_scores": {},
            "descriptors": {
                "qed": round(float(c.get("qed", 0.0) or 0.0), 3),
                "mw": round(float(c.get("mw", 0.0) or 0.0), 2),
                "logp": round(float(c.get("logp", 0.0) or 0.0), 3),
                "tpsa": round(float(c.get("tpsa", 0.0) or 0.0), 2),
                "max_tanimoto": round(float(c.get("max_tanimoto", 0.0)), 3),
            },
        }
        for i, c in enumerate(top)
    ]

    # Persist analog series as denovo_output so Analysis Results > De Novo
    # tab renders it. We include the seed in payload metadata so a future
    # session can reconstruct the lineage (seed → analog batch → ...).
    persistence = _persist_agent_audit(
        action="denovo_output",
        receptor_label=_safe_receptor_label(target, "AGENT-ANALOG"),
        candidate_count=len(analog_rows),
        payload={
            "analysis_section": "denovo",
            "target": target,
            "module": "general",
            "total_generated": len(candidates),
            "after_filter_count": len(similar),
            "top_n": len(analog_rows),
            "seed_scaffold": canonical_scaffold,
            "min_similarity_threshold": min_similarity,
            "results_json": {
                "target": target,
                "module": "general",
                "total_generated": len(candidates),
                "after_filter_count": len(similar),
                "top_n": len(analog_rows),
                "seed_scaffold": canonical_scaffold,
                "data": analog_rows,
            },
            "source": "agent_generate_analogs",
        },
    )

    return {
        "seed_scaffold": canonical_scaffold,
        "prior_used": prior_used,
        "target_name": target,
        "total_generated": len(candidates),
        "passed_similarity_filter": len(similar),
        "near_duplicate_dropped": near_duplicate_count,
        "min_similarity_threshold": min_similarity,
        "returned": len(top),
        "analogs": [
            {
                "smiles": r["smiles"],
                "tanimoto_to_seed": r["tanimoto_to_seed"],
                "qed": r["qed"],
                "mw": r["mw"],
                "logp": r["logp"],
                "tpsa": r["tpsa"],
                "nll": r["nll"],
            }
            for r in analog_rows
        ],
        "next_step_hint": (
            "Recommended: take the top 5–10 analogs, run predict_admet to "
            "drop hERG/DILI failures, then dock_molecules against the same "
            "target to confirm binding is preserved or improved vs the seed."
        ),
        "saved_as_job_id": persistence.get("job_id"),
        "persisted": persistence.get("persisted", False),
    }


async def _tool_fetch_chembl_actives(args: dict[str, Any]) -> dict[str, Any]:
    target = (args.get("target_name") or "").strip()
    if not target:
        return {"error": "target_name is required"}
    smiles = await chembl_client.fetch_known_actives(
        target,
        pchembl_min=float(args.get("pchembl_min") or 7.0),
        limit=int(args.get("limit") or 20),
    )
    return {"target": target, "count": len(smiles), "smiles": smiles[:50]}


async def _tool_list_capabilities(_args: dict[str, Any]) -> dict[str, Any]:
    priors = await reinvent_client.list_available_priors()
    receptors = await docking_client.list_available_targets()
    return {
        "fine_tuned_priors": priors,
        "docking_receptors": receptors,
        "note": "ChEMBL lookup works for any target name; fine-tuned priors and "
        "docking receptors are the limited-set asset lists.",
    }


async def _tool_detect_pocket(args: dict[str, Any]) -> dict[str, Any]:
    target = (args.get("target") or "").strip() or None
    pdb_id = (args.get("pdb_id") or "").strip() or None
    ligand = (args.get("ligand_resname") or "").strip() or None
    if not target and not pdb_id:
        return {"error": "Provide either 'target' or 'pdb_id'."}
    try:
        pocket = await pocket_client.detect_pocket(target=target, pdb_id=pdb_id, ligand_resname=ligand)
    except Exception as exc:  # noqa: BLE001
        return {"error": f"detect_pocket failed: {exc}"}
    if pocket is None:
        return {"error": "Pocket_App unreachable. Start the service on port 8050."}
    return {"target": target, "pdb_id": pdb_id, "pocket": pocket}


async def _tool_extract_pharmacophore(args: dict[str, Any]) -> dict[str, Any]:
    target = (args.get("target") or "").strip() or None
    pdb_id = (args.get("pdb_id") or "").strip() or None
    ligand = (args.get("ligand_resname") or "").strip() or None
    if not target and not pdb_id:
        return {"error": "Provide either 'target' or 'pdb_id'."}
    try:
        resp = await pocket_client.extract_pharmacophore(
            target=target, pdb_id=pdb_id, ligand_resname=ligand
        )
    except Exception as exc:  # noqa: BLE001
        return {"error": f"extract_pharmacophore failed: {exc}"}
    if resp is None:
        return {"error": "Pocket_App unreachable. Start the service on port 8050."}
    # Truncate features list for context efficiency
    features = resp.get("features") or []
    if len(features) > 20:
        features = features[:20]
    return {
        "target": resp.get("target"),
        "pdb_id": resp.get("pdb_id"),
        "ligand_resname": resp.get("ligand_resname"),
        "ligand_smiles": resp.get("ligand_smiles"),
        "feature_count": resp.get("feature_count"),
        "features": features,
        "notes": resp.get("notes"),
    }


_DISPATCH = {
    "generate_molecules": _tool_generate_molecules,
    "predict_admet": _tool_predict_admet,
    "dock_molecules": _tool_dock_molecules,
    "run_full_pipeline": _tool_run_full_pipeline,
    "generate_analogs": _tool_generate_analogs,
    "fetch_chembl_actives": _tool_fetch_chembl_actives,
    "list_capabilities": _tool_list_capabilities,
    "detect_pocket": _tool_detect_pocket,
    "extract_pharmacophore": _tool_extract_pharmacophore,
}


# ---------------------------------------------------------------------------
# Anthropic client (lazy — instantiated per request so a key change in .env
# during dev doesn't require a server restart).
# ---------------------------------------------------------------------------


def _build_client() -> anthropic.AsyncAnthropic:
    cfg = get_settings()
    if not cfg.anthropic_api_key:
        raise RuntimeError(
            "ANTHROPIC_API_KEY is not set. Add it to .env or export it as "
            "an environment variable, then restart the backend."
        )
    return anthropic.AsyncAnthropic(api_key=cfg.anthropic_api_key)


# ---------------------------------------------------------------------------
# Main streaming entrypoint
# ---------------------------------------------------------------------------


def _build_messages(req: ChatRequest) -> list[dict[str, Any]]:
    """Replay client-supplied history (text-only) + the new user turn."""
    msgs: list[dict[str, Any]] = []
    for h in req.history:
        if not h.content:
            continue
        msgs.append({"role": h.role, "content": h.content})
    msgs.append({"role": "user", "content": req.message})
    return msgs


async def stream_chat(req: ChatRequest) -> AsyncIterator[ChatStreamEvent]:
    """Run the agent loop, yielding ChatStreamEvent objects.

    Loop shape:
        1. Send messages → Claude streams text/thinking/tool_use blocks.
        2. If stop_reason != "tool_use", emit ChatDone and return.
        3. Otherwise execute every tool_use, append assistant + tool_result
           to messages, loop.
    """
    cfg = get_settings()
    try:
        client = _build_client()
    except RuntimeError as exc:
        yield ChatError(message=str(exc))
        yield ChatDone(stop_reason="api_key_missing")
        return

    model = req.model or cfg.agent_model
    messages = _build_messages(req)

    # Cache the (large, stable) tool definitions block — this is the part of
    # the request that doesn't change across turns. cache_control on the last
    # tool covers all tools as the prefix.
    tools_with_cache = [dict(t) for t in TOOLS]
    if tools_with_cache:
        tools_with_cache[-1] = {**tools_with_cache[-1], "cache_control": {"type": "ephemeral"}}

    system_blocks = [
        {"type": "text", "text": SYSTEM_PROMPT, "cache_control": {"type": "ephemeral"}}
    ]

    iterations = 0
    total_input_tokens = 0
    total_output_tokens = 0
    stop_reason: str | None = None

    while iterations < cfg.agent_max_iterations:
        iterations += 1
        try:
            request_kwargs: dict[str, Any] = {
                "model": model,
                "max_tokens": cfg.agent_max_tokens,
                "system": system_blocks,
                "tools": tools_with_cache,
                "messages": messages,
            }
            # Opus 4.7 / 4.6 — adaptive thinking with summarised display so the
            # UI gets visible "thinking" deltas instead of empty placeholders.
            if model.startswith(("claude-opus-4-7", "claude-opus-4-6", "claude-sonnet-4-6")):
                request_kwargs["thinking"] = {"type": "adaptive", "display": "summarized"}
            else:
                # Older models — fall back to a normal extended-thinking budget
                request_kwargs["thinking"] = {"type": "enabled", "budget_tokens": 4096}

            assistant_blocks: list[dict[str, Any]] = []
            tool_uses: list[tuple[str, str, dict[str, Any]]] = []

            async with client.messages.stream(**request_kwargs) as stream:
                # Track partial blocks by index so we can assemble them
                current_blocks: dict[int, dict[str, Any]] = {}

                async for event in stream:
                    et = getattr(event, "type", None)

                    if et == "content_block_start":
                        block = event.content_block
                        idx = event.index
                        btype = getattr(block, "type", None)
                        if btype == "text":
                            current_blocks[idx] = {"type": "text", "text": ""}
                        elif btype == "thinking":
                            current_blocks[idx] = {
                                "type": "thinking",
                                "thinking": "",
                                "signature": getattr(block, "signature", None),
                            }
                        elif btype == "tool_use":
                            current_blocks[idx] = {
                                "type": "tool_use",
                                "id": block.id,
                                "name": block.name,
                                "input_buf": "",
                            }

                    elif et == "content_block_delta":
                        idx = event.index
                        delta = event.delta
                        dtype = getattr(delta, "type", None)
                        if dtype == "text_delta":
                            chunk = getattr(delta, "text", "") or ""
                            blk = current_blocks.get(idx)
                            if blk is not None and blk.get("type") == "text":
                                blk["text"] += chunk
                            if chunk:
                                yield ChatTextDelta(text=chunk)
                        elif dtype == "thinking_delta":
                            chunk = getattr(delta, "thinking", "") or ""
                            blk = current_blocks.get(idx)
                            if blk is not None and blk.get("type") == "thinking":
                                blk["thinking"] += chunk
                            if chunk:
                                yield ChatThinkingDelta(text=chunk)
                        elif dtype == "signature_delta":
                            blk = current_blocks.get(idx)
                            if blk is not None and blk.get("type") == "thinking":
                                blk["signature"] = (
                                    (blk.get("signature") or "") + (getattr(delta, "signature", "") or "")
                                )
                        elif dtype == "input_json_delta":
                            chunk = getattr(delta, "partial_json", "") or ""
                            blk = current_blocks.get(idx)
                            if blk is not None and blk.get("type") == "tool_use":
                                blk["input_buf"] += chunk

                    elif et == "content_block_stop":
                        idx = event.index
                        blk = current_blocks.get(idx)
                        if not blk:
                            continue
                        if blk["type"] == "tool_use":
                            try:
                                parsed = json.loads(blk["input_buf"]) if blk["input_buf"] else {}
                            except json.JSONDecodeError:
                                parsed = {}
                            tool_input = parsed if isinstance(parsed, dict) else {}
                            assistant_blocks.append(
                                {
                                    "type": "tool_use",
                                    "id": blk["id"],
                                    "name": blk["name"],
                                    "input": tool_input,
                                }
                            )
                            tool_uses.append((blk["id"], blk["name"], tool_input))
                            yield ChatToolUse(
                                tool_use_id=blk["id"],
                                name=blk["name"],
                                input=tool_input,
                            )
                        elif blk["type"] == "text":
                            assistant_blocks.append({"type": "text", "text": blk["text"]})
                        elif blk["type"] == "thinking":
                            preserved = {"type": "thinking", "thinking": blk["thinking"]}
                            if blk.get("signature"):
                                preserved["signature"] = blk["signature"]
                            assistant_blocks.append(preserved)

                final_msg = await stream.get_final_message()
                stop_reason = final_msg.stop_reason
                if final_msg.usage:
                    total_input_tokens += int(getattr(final_msg.usage, "input_tokens", 0) or 0)
                    total_output_tokens += int(getattr(final_msg.usage, "output_tokens", 0) or 0)

        except anthropic.APIStatusError as exc:
            logger.exception("Anthropic API error")
            yield ChatError(message=f"Anthropic API error ({exc.status_code}): {exc.message}")
            yield ChatDone(stop_reason="api_error", iterations=iterations)
            return
        except anthropic.APIConnectionError as exc:
            logger.exception("Anthropic connection error")
            yield ChatError(message=f"Anthropic connection failed: {exc}")
            yield ChatDone(stop_reason="connection_error", iterations=iterations)
            return
        except Exception as exc:  # noqa: BLE001
            logger.exception("Unexpected agent loop error")
            yield ChatError(message=f"agent loop failed: {exc}")
            yield ChatDone(stop_reason="internal_error", iterations=iterations)
            return

        # Append the assistant turn to the conversation in full (incl. thinking
        # blocks — required by Claude when continuing after a tool call).
        if assistant_blocks:
            messages.append({"role": "assistant", "content": assistant_blocks})

        if stop_reason != "tool_use":
            break

        # Execute every tool_use block, append a single user-role tool_result message.
        tool_result_content: list[dict[str, Any]] = []
        for tool_use_id, name, tool_input in tool_uses:
            handler = _DISPATCH.get(name)
            if handler is None:
                payload = {"error": f"unknown tool {name!r}"}
                is_error = True
            else:
                try:
                    payload = await handler(tool_input)
                    is_error = bool(payload.get("error")) if isinstance(payload, dict) else False
                except Exception as exc:  # noqa: BLE001
                    logger.exception("Tool %s crashed", name)
                    payload = {"error": f"{type(exc).__name__}: {exc}"}
                    is_error = True

            content_str = json.dumps(payload, ensure_ascii=False, default=str)
            # Defensive truncation — Claude sees this verbatim; oversized
            # payloads waste cache slots and can blow context.
            if len(content_str) > 60_000:
                content_str = content_str[:60_000] + " …[truncated]"

            tool_result_content.append(
                {
                    "type": "tool_result",
                    "tool_use_id": tool_use_id,
                    "content": content_str,
                    "is_error": is_error,
                }
            )
            yield ChatToolResult(
                tool_use_id=tool_use_id,
                name=name,
                is_error=is_error,
                content=content_str,
            )

        messages.append({"role": "user", "content": tool_result_content})
        # loop continues — Claude will see the tool results next iteration

    if iterations >= cfg.agent_max_iterations and stop_reason == "tool_use":
        yield ChatError(message=f"agent stopped: hit max_iterations={cfg.agent_max_iterations}")
        stop_reason = "max_iterations"

    yield ChatDone(
        stop_reason=stop_reason,
        iterations=iterations,
        input_tokens=total_input_tokens or None,
        output_tokens=total_output_tokens or None,
    )
