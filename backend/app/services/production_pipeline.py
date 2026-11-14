"""Production discovery pipeline: REINVENT 4 generator + ADMET-AI predictor.

Replaces the synthetic ``universal_discovery_engine.cjs`` motor with real
chemistry. Output is shape-compatible with ``DiscoveryResponse`` so the
existing UI consumes it without changes.

Flow:
    1. Reinvent_App → N raw SMILES + NLL + RDKit basic descriptors (mw, logp, qed, tpsa)
    2. Filter invalid / non-RDKit-parsable SMILES
    3. Admet_App → 52 ADMET endpoints + 52 DrugBank percentiles per molecule
    4. Apply user filter_overrides (base / safety / pk) using mapped ADMET fields
    5. Heuristic score from available signals; never interpret as activity probability
    6. Sort, slice top-N, format as DiscoveryResponse-compatible candidates
"""

from __future__ import annotations

import logging
import math
import hashlib
import json
from uuid import uuid4
from datetime import datetime, timezone
from threading import Lock
from typing import Any, Iterable, Optional

from .admet_client import predict_admet
from .chembl_client import fetch_known_actives
from .docking_client import dock_batch as run_docking, list_available_targets
from .reinvent_client import list_available_priors, sample_smiles

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# RDKit substructure-alert catalogs (PAINS / BRENK) and on-demand descriptors
# ---------------------------------------------------------------------------
# These are built lazily on first use to avoid paying the RDKit import cost at
# server startup (RDKit is ~100 MB of native deps). Both catalogs are
# thread-safe to query but we still guard initialization with a Lock to keep
# the first-call race window tight.
_FILTER_CATALOG_LOCK = Lock()
_PAINS_CATALOG: Optional[Any] = None
_BRENK_CATALOG: Optional[Any] = None


def _get_pains_catalog() -> Any:
    """Lazy-build the RDKit PAINS FilterCatalog (singleton).

    PAINS = Pan-Assay Interference Compounds. Catalog ships with RDKit and
    matches the original Baell & Holloway 2010 substructure set (~480
    SMARTS). A positive match means the molecule is *likely* an assay
    artifact (aggregator, redox-cycler, reactive, fluorescer, etc.) and
    pharma medicinal chemistry teams routinely reject such hits before
    advancing them. We expose this as a hard filter so the De Novo top-N
    can be trusted as "real" drug candidates.
    """
    global _PAINS_CATALOG
    if _PAINS_CATALOG is not None:
        return _PAINS_CATALOG
    with _FILTER_CATALOG_LOCK:
        if _PAINS_CATALOG is None:
            from rdkit.Chem import FilterCatalog
            params = FilterCatalog.FilterCatalogParams()
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
            _PAINS_CATALOG = FilterCatalog.FilterCatalog(params)
    return _PAINS_CATALOG


def _get_brenk_catalog() -> Any:
    """Lazy-build the RDKit BRENK FilterCatalog (singleton).

    BRENK = Brenk et al. 2008 substructure alert set (~150 patterns). Flags
    chemotypes with known liabilities: nitro groups (mutagenic), aromatic
    amines, vinyl chlorides, peroxides, hydrazines, etc. Together with PAINS
    this is the de-facto "drug-like-ish floor" used in academic + early
    industrial screening cascades.
    """
    global _BRENK_CATALOG
    if _BRENK_CATALOG is not None:
        return _BRENK_CATALOG
    with _FILTER_CATALOG_LOCK:
        if _BRENK_CATALOG is None:
            from rdkit.Chem import FilterCatalog
            params = FilterCatalog.FilterCatalogParams()
            params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
            _BRENK_CATALOG = FilterCatalog.FilterCatalog(params)
    return _BRENK_CATALOG


def _structural_alerts(smiles: str) -> tuple[list[str], list[str]]:
    """Return (hard_failures, soft_warnings) for substructure-alert catalogs.

    - **PAINS hits** go to ``hard_failures`` and the upstream caller rejects
      the molecule. PAINS captures assay-interference / nuisance chemotypes
      that pharma DD always reject (aggregators, redox-cyclers, reactives).

    - **BRENK hits** go to ``soft_warnings`` only. They are surfaced in the
      candidate's output for medicinal-chemistry review but do **not** filter
      the molecule out, because BRENK is far too aggressive to use as a hard
      gate: it rejects aspirin, metformin, osimertinib (acrylamide warhead),
      remibrutinib (cyanoacrylamide), and most other covalent inhibitors.
      Industry convention treats BRENK as an annotation, not a rejection.

    Returns ``([], [])`` on any RDKit parse failure (defensive no-op; the
    upstream pipeline already filters parse failures).
    """
    try:
        from rdkit import Chem
    except Exception:  # noqa: BLE001
        return [], []
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return [], []
    hard: list[str] = []
    soft: list[str] = []
    try:
        if _get_pains_catalog().HasMatch(mol):
            hard.append("PAINS")
    except Exception:  # noqa: BLE001
        pass
    try:
        if _get_brenk_catalog().HasMatch(mol):
            soft.append("BRENK")
    except Exception:  # noqa: BLE001
        pass
    return hard, soft


def _extra_rdkit_descriptors(smiles: str) -> dict[str, float]:
    """Compute on-demand RDKit descriptors not already in the Reinvent_App payload.

    Reinvent_App returns ``mw``, ``logp``, ``qed``, ``tpsa``. The UI filter
    set also covers HBD/HBA, rotatable bonds, formal charge, and synthetic
    accessibility (SA). We compute those here so the filter engine can
    enforce every UI input — previously several UI fields were silently
    ignored. Returns an empty dict on any failure so callers can keep
    checking the filters that don't need these extras.
    """
    out: dict[str, float] = {}
    try:
        from rdkit import Chem
        from rdkit.Chem import Lipinski, Descriptors, rdMolDescriptors
    except Exception:  # noqa: BLE001
        return out
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return out
    try:
        out["hbd"] = float(Lipinski.NumHDonors(mol))
        out["hba"] = float(Lipinski.NumHAcceptors(mol))
        out["rotatable_bonds"] = float(Lipinski.NumRotatableBonds(mol))
        out["formal_charge_abs"] = float(abs(Chem.GetFormalCharge(mol)))
    except Exception:  # noqa: BLE001
        pass
    # SA score is in RDKit Contrib (sascorer.py). Not all installs ship it.
    try:
        from rdkit.Chem import RDConfig  # type: ignore
        import os, sys
        sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_path not in sys.path:
            sys.path.append(sa_path)
        import sascorer  # type: ignore
        out["sa_score"] = float(sascorer.calculateScore(mol))
    except Exception:  # noqa: BLE001
        # SA filter just won't activate for this molecule — graceful skip.
        pass
    return out

# ---------------------------------------------------------------------------
# Filter spec → ADMET-AI field mapping
# ---------------------------------------------------------------------------
# Map the user-facing filter keys (from the existing UI) to ADMET-AI raw keys.
# Comparator: "max" means value <= threshold; "min" means value >= threshold;
#             "range" means [low, high] both inclusive.

SAFETY_MAX_MAP = {
    "herg_risk_max":     "hERG",
    "dili_risk_max":     "DILI",
    "bsep_risk_max":     None,  # ADMET-AI has no BSEP endpoint; ignored gracefully
    "cyp3a4_inhib_max":  "CYP3A4_Veith",
    "cyp2d6_inhib_max":  "CYP2D6_Veith",
    "cyp2c9_inhib_max":  "CYP2C9_Veith",
    "ames_max":          "AMES",
}

PK_MIN_MAP = {
    # NOTE: ADMET-AI's Caco2_Wang is log10(cm/s) ~ [-7..-4]; UI eşiği 0..1
    # aralığında olduğu için doğrudan karşılaştırma yapamayız. Filter'ı
    # normalize'lı olarak _check_filters içinde uyguluyoruz (aşağıda).
    "oral_bioavailability_min": "Bioavailability_Ma",  # 0..1 olasılık → eşik aynı birimde
}
PK_MAX_MAP: dict[str, str] = {
    # PPB_AZ ADMET-AI'da yüzde olarak (0..100), UI eşiği 0..1 oranı —
    # _check_filters içinde 100'e bölerek normalize ediyoruz.
}
# Caco2_Wang ve PPBR_AZ için manuel normalize edilmiş özel kontroller
# _check_filters içinde; eski "doğrudan karşılaştırma" yaklaşımı birim
# uyumsuzluğu nedeniyle her molekülü eliyordu.

# Clearance ve half-life için: ADMET-AI değerleri farklı birimde, varsayılan
# aralıklar (clearance 0.2-0.8) gerçek değerleri elediği için bu range
# kontrollerini şimdilik DEVRE DIŞI bırakıyoruz. UI'da kalsa bile pipeline
# uygulamaz. Doğru birim normalizasyonu Phase 2.
PK_RANGE_MAP: dict[str, tuple[str, str, str]] = {}


class UnsupportedFilterError(ValueError):
    pass


def validate_supported_filters(overrides):
    """Fail closed on thresholds that cannot be enforced or are not numeric.

    The legacy demo engine has a different schema. This check is only for the
    production scientific workflow, so callers must explicitly choose demo.
    """
    if not overrides:
        return
    if not isinstance(overrides, dict):
        raise UnsupportedFilterError("filter_overrides must be a JSON object")
    supported = {
        "base_filters": {"qed_min", "sa_score_max", "mw_min", "mw_max", "logp_min", "logp_max", "tpsa_min", "tpsa_max", "rotatable_bonds_max", "hbd_max", "hba_max", "formal_charge_abs_max"},
        "safety_filters": set(SAFETY_MAX_MAP) - {"bsep_risk_max"},
        "pk_filters": set(PK_MIN_MAP) | set(PK_MAX_MAP) | {"caco2_min", "ppb_max"},
        "docking": {"enabled", "score_threshold"},
    }
    unsupported = []
    import math
    for category, payload in overrides.items():
        if category not in supported or not isinstance(payload, dict):
            unsupported.append(str(category))
            continue
        for key, value in payload.items():
            if key not in supported[category]:
                unsupported.append(f"{category}.{key}")
                continue
            if category == "docking" and key == "enabled":
                if not isinstance(value, bool):
                    unsupported.append(f"{category}.{key}: must be boolean")
                continue
            try:
                numeric = float(value)
                if isinstance(value, bool) or not math.isfinite(numeric):
                    raise ValueError("non-finite/non-numeric value")
            except (TypeError, ValueError, OverflowError):
                unsupported.append(f"{category}.{key}: must be finite numeric")
                continue
            if category == "docking" and key == "score_threshold" and numeric != 0:
                unsupported.append("docking.score_threshold: nonzero threshold is not implemented")
    if unsupported:
        raise UnsupportedFilterError("Unsupported model filter(s): " + ", ".join(sorted(unsupported)))


def _normalize_nll(nll: float) -> float:
    """Map REINVENT NLL → [0,1] where higher means more prior-likely (= more drug-like skeleton).

    Empirically, REINVENT prior NLLs span ~15..60 for valid SMILES. Linear clamp to [0,1].
    """
    lo, hi = 15.0, 60.0
    value = float(nll)
    if not math.isfinite(value):
        raise ValueError("REINVENT likelihood must be finite")
    x = max(lo, min(hi, value))
    return 1.0 - (x - lo) / (hi - lo)


def _check_filters(
    rdkit_desc: dict[str, float],
    admet_raw: dict[str, float],
    filter_overrides: dict[str, dict[str, Any]] | None,
    smiles: str | None = None,
) -> tuple[list[str], list[str]]:
    """Return ``(failed, warnings)`` — hard rejections + soft annotations.

    ``failed`` is the list of filter names the molecule violated; a non-empty
    failed list means the caller should reject the candidate. ``warnings`` is
    a parallel list of annotations that do NOT block (e.g. BRENK substructure
    alerts — informational, not rejecting). Survivors store the warnings on
    their candidate output for downstream medicinal-chemistry review.

    ``smiles`` is optional but required to enforce the structural-alert
    filters (PAINS hard / BRENK soft) and the on-demand RDKit descriptors
    (HBD, HBA, rotatable bonds, formal charge, SA score). When provided,
    every UI filter input becomes enforceable; when omitted, those checks
    are silently skipped (legacy behavior preserved).
    """
    if not filter_overrides:
        return [], []
    bf = filter_overrides.get("base_filters") or {}
    sf = filter_overrides.get("safety_filters") or {}
    pf = filter_overrides.get("pk_filters") or {}
    failed: list[str] = []
    warnings: list[str] = []

    # A specified safety/PK threshold MUST NOT silently pass when a model did
    # not predict the endpoint. It is rejected with an explicit missing label.
    required_admet = {
        **{key: field for key, field in SAFETY_MAX_MAP.items() if field is not None},
        **PK_MIN_MAP, **PK_MAX_MAP,
        "caco2_min": "Caco2_Wang", "ppb_max": "PPBR_AZ",
    }
    def has_finite_value(mapping: dict, name: str) -> bool:
        """A NaN/Inf must never satisfy a numeric scientific threshold."""
        value = mapping.get(name)
        try:
            return not isinstance(value, bool) and math.isfinite(float(value))
        except (TypeError, ValueError, OverflowError):
            return False

    for configured, field in required_admet.items():
        section = sf if configured in sf else pf
        if configured in section and not has_finite_value(admet_raw, field):
            failed.append(f"{field}_missing")

    # Apply the same rule to mandatory RDKit descriptors; None is not a pass.
    base_fields = {
        "qed_min": "qed", "mw_min": "mw", "mw_max": "mw",
        "logp_min": "logp", "logp_max": "logp",
        "tpsa_min": "tpsa", "tpsa_max": "tpsa",
    }
    for key, desc in base_fields.items():
        if key in bf and not has_finite_value(rdkit_desc, desc):
            failed.append(f"{desc}_missing")

    # Base filters — operate on RDKit descriptors (mw, logp, qed, tpsa)
    mw = rdkit_desc.get("mw")
    logp = rdkit_desc.get("logp")
    qed = rdkit_desc.get("qed")
    tpsa = rdkit_desc.get("tpsa")
    if "qed_min" in bf and qed is not None and qed < float(bf["qed_min"]):
        failed.append("qed")
    if "mw_min" in bf and mw is not None and mw < float(bf["mw_min"]):
        failed.append("mw_low")
    if "mw_max" in bf and mw is not None and mw > float(bf["mw_max"]):
        failed.append("mw_high")
    if "logp_min" in bf and logp is not None and logp < float(bf["logp_min"]):
        failed.append("logp_low")
    if "logp_max" in bf and logp is not None and logp > float(bf["logp_max"]):
        failed.append("logp_high")
    if "tpsa_min" in bf and tpsa is not None and tpsa < float(bf["tpsa_min"]):
        failed.append("tpsa_low")
    if "tpsa_max" in bf and tpsa is not None and tpsa > float(bf["tpsa_max"]):
        failed.append("tpsa_high")

    # Extra RDKit descriptors — HBD, HBA, rotatable bonds, formal charge, SA.
    # These were UI fields silently ignored by the original implementation;
    # we now compute them on-demand from SMILES so every input does something.
    extra_filters = ("hbd_max", "hba_max", "rotatable_bonds_max", "formal_charge_abs_max", "sa_score_max")
    if smiles is None:
        # An API caller lacking molecular structures cannot enforce these
        # controls. Explicitly reject the candidate rather than skipping them.
        for key, extra_name in {
            "hbd_max": "hbd", "hba_max": "hba", "rotatable_bonds_max": "rotatable_bonds",
            "formal_charge_abs_max": "formal_charge_abs", "sa_score_max": "sa_score"
        }.items():
            if key in bf:
                failed.append(f"{extra_name}_missing")
    if smiles is not None and any(k in bf for k in extra_filters):
        extras = _extra_rdkit_descriptors(smiles)
        hbd = extras.get("hbd")
        hba = extras.get("hba")
        rb = extras.get("rotatable_bonds")
        chg = extras.get("formal_charge_abs")
        sa = extras.get("sa_score")
        for key, extra_key in {
            "hbd_max": "hbd", "hba_max": "hba",
            "rotatable_bonds_max": "rotatable_bonds",
            "formal_charge_abs_max": "formal_charge_abs",
            "sa_score_max": "sa_score",
        }.items():
            if key in bf and extras.get(extra_key) is None:
                failed.append(f"{extra_key}_missing")
        if "hbd_max" in bf and hbd is not None and hbd > float(bf["hbd_max"]):
            failed.append("hbd_high")
        if "hba_max" in bf and hba is not None and hba > float(bf["hba_max"]):
            failed.append("hba_high")
        if "rotatable_bonds_max" in bf and rb is not None and rb > float(bf["rotatable_bonds_max"]):
            failed.append("rotatable_bonds_high")
        if "formal_charge_abs_max" in bf and chg is not None and chg > float(bf["formal_charge_abs_max"]):
            failed.append("formal_charge_high")
        if "sa_score_max" in bf and sa is not None and sa > float(bf["sa_score_max"]):
            failed.append("sa_score_high")

    # Structural-alert catalogs — PAINS (hard reject) + BRENK (soft warning).
    # See _structural_alerts() docstring for why BRENK is informational only:
    # it rejects aspirin / metformin / osimertinib / remibrutinib and most
    # covalent inhibitors, so industry treats it as annotation, not gate.
    if smiles is not None:
        hard_alerts, soft_alerts = _structural_alerts(smiles)
        failed.extend(hard_alerts)
        warnings.extend(soft_alerts)

    # Safety filters — ADMET-AI raw fields, "max" semantics
    for key, admet_field in SAFETY_MAX_MAP.items():
        if admet_field is None or key not in sf:
            continue
        threshold = float(sf[key])
        v = admet_raw.get(admet_field)
        if v is not None and v > threshold:
            failed.append(admet_field)

    # PK filters — same-unit fields handled generically
    for key, admet_field in PK_MIN_MAP.items():
        if key not in pf:
            continue
        threshold = float(pf[key])
        v = admet_raw.get(admet_field)
        if v is not None and v < threshold:
            failed.append(f"{admet_field}_low")
    for key, admet_field in PK_MAX_MAP.items():
        if key not in pf:
            continue
        threshold = float(pf[key])
        v = admet_raw.get(admet_field)
        if v is not None and v > threshold:
            failed.append(f"{admet_field}_high")

    # Caco2_Wang: ADMET-AI returns log10(cm/s) ~ [-7..-4]. Map UI's 0..1
    # threshold to the log scale: 0.0 → -7.0 (no permeability), 1.0 → -4.0
    # (high permeability). So caco2_min=0.5 means "log10 >= -5.5".
    if "caco2_min" in pf:
        threshold_unit = float(pf["caco2_min"])  # 0..1
        threshold_log = -7.0 + 3.0 * max(0.0, min(1.0, threshold_unit))
        v = admet_raw.get("Caco2_Wang")
        if v is not None and float(v) < threshold_log:
            failed.append("Caco2_Wang_low")

    # PPBR_AZ: ADMET-AI returns percentage (0..100). UI threshold is 0..1
    # ratio. Convert before compare: ppb_max=0.95 → percent 95.
    if "ppb_max" in pf:
        threshold_pct = 100.0 * float(pf["ppb_max"])
        v = admet_raw.get("PPBR_AZ")
        if v is not None and float(v) > threshold_pct:
            failed.append("PPBR_AZ_high")

    # Clearance / half-life range checks intentionally disabled — see comment
    # above the (empty) PK_RANGE_MAP. The ADMET-AI units (mL/min/kg, hours)
    # don't align cleanly with the existing UI ranges and naive comparison
    # eliminated every candidate.

    return failed, warnings


def _normalize_dock_affinity(affinity_kcal_mol: float | None) -> float | None:
    """Map Vina scoring output to [0,1]; this is not measured affinity.

    Linear clamp: -12 → 1.0 (very strong), -4 → 0.0 (weak/none).
    Returns None if no affinity (docking failed).
    """
    if affinity_kcal_mol is None:
        return None
    lo, hi = -12.0, -4.0  # most negative = strongest binding
    a = max(lo, min(hi, float(affinity_kcal_mol)))
    return (hi - a) / (hi - lo)  # -12 -> 1.0, -4 -> 0.0


# ---------------------------------------------------------------------------
# Novelty / IP-defensibility thresholds.
# Tanimoto similarity to closest known active (ECFP4) is the proxy for how
# "new" a generated molecule is relative to ChEMBL prior art.
#
#   tanimoto >= NOVELTY_HARD_DROP  →  drop near-identical fingerprints relative
#                                     to the *available reference set* only
#   NOVELTY_SOFT_LO <= t < hard    →  keep but penalise score linearly toward
#                                     0.5× at the hard threshold ("high-similarity"
#                                     region of the reference set)
#   tanimoto < NOVELTY_SOFT_LO     →  no penalty (real novel chemistry)
# ---------------------------------------------------------------------------
NOVELTY_HARD_DROP = 0.99   # essentially identical fingerprint
NOVELTY_SOFT_LO = 0.85     # heuristic similarity cutoff; NOT a legal/IP test


def _novelty_penalty(max_tanimoto: float | None) -> float:
    """Return a [0.5, 1.0] multiplier applied to the final score.

    1.0 = no penalty relative to this reference set (NOT proof of global novelty).
    Linearly ramps to 0.5 between SOFT_LO and HARD_DROP.
    Molecules at or above HARD_DROP are filtered upstream and never reach here.
    """
    if max_tanimoto is None:
        return 1.0
    t = float(max_tanimoto)
    if t < NOVELTY_SOFT_LO:
        return 1.0
    if t >= NOVELTY_HARD_DROP:
        return 0.5  # caller should have already dropped — defensive
    # linear from 1.0 at NOVELTY_SOFT_LO to 0.5 at NOVELTY_HARD_DROP
    span = NOVELTY_HARD_DROP - NOVELTY_SOFT_LO
    return 1.0 - 0.5 * (t - NOVELTY_SOFT_LO) / span


def _component_scores(
    rdkit_desc: dict[str, float],
    admet_raw: dict[str, float],
    nll: float,
    max_tanimoto: float | None = None,
    dock_affinity: float | None = None,
) -> dict[str, float]:
    """Map per-molecule signals into [0,1] components for the weighted final score.

    When ``max_tanimoto`` (similarity to known target actives) is available,
    we expose it as a separate ``target_relevance`` component so the final
    score reflects target-awareness. We also derive ``novelty = 1 -
    max_tanimoto`` so the model is rewarded for finding genuinely new chemistry
    rather than reproducing literature.

    Docking contributes only when supplied by the ranking policy. Generator
    NLL never substitutes for missing binding evidence.
    """
    qed = rdkit_desc.get("qed") or 0.0
    real_binding = _normalize_dock_affinity(dock_affinity)
    # Generator likelihood is not target binding evidence. Missing docking
    # stays absent rather than being replaced by a generation score.
    # Safety: invert risk probabilities (lower risk → higher safety score)
    risks = []
    for f in ("hERG", "DILI", "AMES", "CYP3A4_Veith"):
        v = admet_raw.get(f)
        if v is not None:
            risks.append(float(v))
    safety = 1.0 - (sum(risks) / len(risks)) if risks else 0.5
    # ADMET aggregate: oral bioavailability + Caco2 (normalized) + (1 - PPB normalized)
    bioavail = admet_raw.get("Bioavailability_Ma")
    bioavail = float(bioavail) if bioavail is not None else 0.5
    caco2 = admet_raw.get("Caco2_Wang")
    # Caco2_Wang is log10(cm/s), typical range ~ -7 to -4, normalize
    caco2_norm = max(0.0, min(1.0, (float(caco2) + 7.0) / 3.0)) if caco2 is not None else 0.5
    ppb = admet_raw.get("PPBR_AZ")
    ppb_norm = max(0.0, min(1.0, 1.0 - float(ppb) / 100.0)) if ppb is not None else 0.5
    admet_score = (bioavail + caco2_norm + ppb_norm) / 3.0
    target_relevance = float(max_tanimoto) if max_tanimoto is not None else 0.5
    # Real novelty score (was previously a 0.5 placeholder). Higher Tanimoto
    # to known actives = lower novelty. When no reference set, default to 0.5.
    novelty = (1.0 - float(max_tanimoto)) if max_tanimoto is not None else 0.5
    return {
        **({"binding": round(real_binding, 4)} if real_binding is not None else {}),
        "admet": round(max(0.0, min(1.0, admet_score)), 4),
        "safety": round(max(0.0, min(1.0, safety)), 4),
        "druglikeness": round(max(0.0, min(1.0, qed)), 4),
        "target_relevance": round(max(0.0, min(1.0, target_relevance)), 4),
        "novelty": round(max(0.0, min(1.0, novelty)), 4),
    }


def _final_score(
    components: dict[str, float],
    target_aware: bool,
    max_tanimoto: float | None = None,
) -> float:
    """Weighted sum × novelty penalty.

    Target-aware weights rebalance toward novelty (was 0.03, now 0.15) so
    reference-set dissimilarity matters more than "found a known one again".
    The novelty penalty multiplier additionally squeezes scores in the
    near-exact-copy zone (Tanimoto > 0.85), so even if other signals are
    great a literature-derivative still ranks below a true novel hit.
    """
    if target_aware:
        weights = {
            "target_relevance": 0.20,    # was 0.30 — still rewarded but less dominant
            "binding": 0.20,
            "admet": 0.20,
            "safety": 0.20,
            "druglikeness": 0.05,
            "novelty": 0.15,             # was 0.03 — now meaningful
        }
    else:
        weights = {
            "binding": 0.30, "admet": 0.25, "safety": 0.25,
            "druglikeness": 0.15, "novelty": 0.05,
        }
    available_weight = sum(w for k, w in weights.items() if k in components)
    base = sum(components.get(k, 0.0) * w for k, w in weights.items()) / available_weight if available_weight else 0.0
    penalty = _novelty_penalty(max_tanimoto)
    return round(base * penalty, 4)


def _format_candidate(
    rank: int,
    smi: str,
    nll: float,
    rdkit_desc: dict,
    admet_raw: dict,
    alerts: list[str],
    module: str,
    max_tanimoto: float | None = None,
    target_aware: bool = False,
    dock_affinity: float | None = None,
    dock_pose_pdb: str | None = None,
    dock_pose_pdb_all: str | None = None,
    docking_for_ranking: bool = True,
) -> dict:
    components = _component_scores(
        rdkit_desc, admet_raw, nll,
        max_tanimoto=max_tanimoto,
        dock_affinity=dock_affinity if docking_for_ranking else None,
    )
    score = _final_score(components, target_aware=target_aware, max_tanimoto=max_tanimoto)
    # dock_score field on DiscoveryResponse: ONLY populated when Vina actually
    # ran against a bundled receptor for this target. When docking is
    # unavailable (target not bundled, Docking_App down, individual candidate
    # failed to dock) we set it to None rather than fabricating a placeholder.
    # The previous "-10.0 * components['binding']" formula was an NLL-derived
    # surrogate that looked like real kcal/mol in the UI and was at risk of
    # being misread as genuine binding affinity in pharma DD — that synthetic
    # holdover is removed here. Downstream consumers (UI table, ranking) must
    # treat None as "no docking; rank from ADMET + similarity + drug-likeness".
    if dock_affinity is not None:
        dock_score_field = round(float(dock_affinity), 3)
    else:
        dock_score_field = None
    def observed(source, key):
        value = source.get(key)
        if value is None:
            return None
        try:
            return round(float(value), 4)
        except (ValueError, TypeError):
            return None

    # Missing ADMET fields are genuinely unavailable, not a measured zero.
    descriptors = {
        "qed": observed(rdkit_desc, "qed"),
        "mw": observed(rdkit_desc, "mw"),
        "logp": observed(rdkit_desc, "logp"),
        "tpsa": observed(rdkit_desc, "tpsa"),
        "hERG": observed(admet_raw, "hERG"),
        "DILI": observed(admet_raw, "DILI"),
        "AMES": observed(admet_raw, "AMES"),
        "BBB": observed(admet_raw, "BBB_Martins"),
        "Bioavailability": observed(admet_raw, "Bioavailability_Ma"),
        "Half_Life_h": observed(admet_raw, "Half_Life_Obach"),
        "nll": round(float(nll), 4),
    }
    if max_tanimoto is not None:
        descriptors["max_tanimoto"] = round(float(max_tanimoto), 4)
    if dock_affinity is not None:
        descriptors["vina_affinity_kcal_mol"] = round(float(dock_affinity), 3)
    out = {
        "rank": rank,
        "smiles": smi,
        "module": module,
        "score": score,
        "dock_score": dock_score_field,
        "alerts": alerts,
        "component_scores": components,
        "descriptors": descriptors,
    }
    # Persist the docked-pose PDB block alongside the row when Vina
    # returned one. The NGL viewer fetches this via
    # /api/jobs/{id}/poses/{rank}.pdb to render the real binding geometry.
    # Stored as a string here, not a file — the audit log payload is JSON
    # so a 3-15 KB string per pose is the simplest path; for a 30-row
    # top-k that's ~150 KB of payload bloat, still well within SQLite's
    # comfort zone (default page size 4 KB, no hard row-size cap on
    # SQLite's TEXT column).
    if dock_pose_pdb:
        out["pose_pdb"] = dock_pose_pdb
    if dock_pose_pdb_all:
        out["pose_pdb_all"] = dock_pose_pdb_all
    return out


async def _resolve_reference_set(req: dict, target: str) -> tuple[list[str], str]:
    """Build the target-relevance reference set for this run.

    Priority:
      1. If caller passed ``known_actives_smiles`` explicitly, use that list verbatim.
      2. Else, if ``use_chembl`` is truthy (default), fetch potent actives from ChEMBL.
      3. Else, return empty list (no target-aware scoring).

    Returns (smiles_list, source_label).
    """
    explicit = req.get("known_actives_smiles") or []
    if isinstance(explicit, list) and explicit:
        cleaned = [s for s in (str(x).strip() for x in explicit) if s]
        if cleaned:
            return cleaned, "user_provided"

    use_chembl = req.get("use_chembl", True)
    if not use_chembl or not target or target == "TARGET":
        return [], "none"

    pchembl_min = float(req.get("chembl_pchembl_min", 7.0))
    limit = int(req.get("chembl_limit", 50))
    smiles = await fetch_known_actives(target, pchembl_min=pchembl_min, limit=limit)
    return smiles, ("chembl" if smiles else "chembl_no_match")


class TargetValidationError(Exception):
    """Raised when target_name has zero meaningful signal across all three
    target-aware axes (ChEMBL actives, fine-tuned REINVENT prior, bundled
    Vina receptor). The caller can choose to translate this into a 422
    HTTPException, or pass ``force_generic=True`` in the request to bypass
    the gate and run a target-agnostic drug-like generation instead.

    Carries the per-axis findings so the API can surface "ChEMBL: 0 hits,
    no prior, no receptor" in the error detail — a silent fallback to
    generic drug-likeness was the bug that prompted this gate (a user
    typed ``target_name="KLM17XY"`` and got back 15 "candidates" that
    were just random drug-like molecules with no target relevance).
    """

    def __init__(self, target: str, findings: dict):
        self.target = target
        self.findings = findings
        super().__init__(
            f"target_name {target!r} produced no target-aware signal: "
            f"chembl_hits={findings.get('chembl_hits', 0)}, "
            f"fine_tuned_prior={findings.get('has_fine_tuned_prior', False)}, "
            f"bundled_receptor={findings.get('has_bundled_receptor', False)}"
        )


async def run_production_discovery(req: dict) -> dict:
    """End-to-end pipeline. Returns a dict matching DiscoveryResponse schema."""
    target = (req.get("target_name") or "").strip() or "TARGET"
    # The source-block control study did not validate ALK5 docking ranking.
    # Keep poses/scores for inspection but exclude them from prioritization.
    alk5_target = target.upper().replace("-", "").replace(" ", "") in {
        "ALK5", "TGFBR1", "CHEMBL4439", "TGFBETARECEPTORTYPE1",
    }
    module = req.get("disease_module") or "general"
    n = int(req.get("generate_count") or 400)
    n = max(50, min(3000, n))
    seed = req.get("seed")  # optional integer; none means non-deterministic REINVENT sampling
    overrides = req.get("filter_overrides")
    validate_supported_filters(overrides)
    return_top_n = 20
    # When True the caller knowingly accepts a target-agnostic run; the
    # pipeline still produces molecules but the engine label and meta
    # block flag the run as "drug-like only" so downstream UI shows the
    # warning banner instead of pretending Vina + ChEMBL Tanimoto ran.
    force_generic = bool(req.get("force_generic", False))

    # 0. Resolve target-relevance reference set (ChEMBL actives or user list)
    reference_smiles, ref_source = await _resolve_reference_set(req, target)
    target_aware = bool(reference_smiles)
    logger.info(
        "Production pipeline: target=%r reference_source=%s reference_count=%d",
        target, ref_source, len(reference_smiles),
    )

    # 0b. Validation gate — refuse to silently degrade into a target-agnostic
    # run. We require at least ONE of the three target-aware signals before
    # the pipeline will proceed:
    #   (a) ChEMBL has potent actives for the target name
    #   (b) A fine-tuned REINVENT prior exists for the target
    #   (c) Vina has a bundled receptor for the target
    # Missing all three almost always means the target_name is a typo or
    # made-up string. Without ``force_generic=True`` we surface this as a
    # validation error so the API can return 422 and the UI can prompt
    # for correction — much better than silently producing 15 generic
    # drug-like molecules labeled "reinvent4_admet_vina" that have zero
    # connection to the user's stated target.
    if not force_generic and target != "TARGET":
        try:
            available_priors = list(await list_available_priors())
        except Exception:  # noqa: BLE001
            available_priors = []
        has_fine_tuned_prior = any(
            p.lower() == target.lower() for p in available_priors if p.lower() != "base"
        )
        try:
            available_receptors = list(await list_available_targets())
        except Exception:  # noqa: BLE001
            available_receptors = []
        has_bundled_receptor = any(
            r.lower() == target.lower() for r in available_receptors
        )
        chembl_hits = len(reference_smiles)
        if chembl_hits == 0 and not has_fine_tuned_prior and not has_bundled_receptor:
            raise TargetValidationError(
                target=target,
                findings={
                    "chembl_hits": chembl_hits,
                    "has_fine_tuned_prior": has_fine_tuned_prior,
                    "has_bundled_receptor": has_bundled_receptor,
                    "available_fine_tuned_priors": [p for p in available_priors if p.lower() != "base"],
                    "available_bundled_receptors": available_receptors,
                },
            )

    # 1. Generate (passing reference_smiles enables Tanimoto scoring at the sampler)
    # Use a target-specific fine-tuned prior if one exists; fall back to base.
    # Caller can override the auto-resolution by passing ``target_prior``
    # explicitly in the request body — used by the UI dropdown that lets
    # the chemist pick between "auto-by-target-name", "base", or a
    # specifically-named promoted RL checkpoint (e.g. "EGFR_rl").
    requested_prior = (req.get("target_prior") or "").strip() or target
    sampled, prior_used = await sample_smiles(
        batch_size=n,
        seed=seed,
        validate_with_rdkit=True,
        reference_smiles=reference_smiles or None,
        target_prior=requested_prior,  # Reinvent_App matches by name, falls back gracefully
    )
    valid_only = [s for s in sampled if s.get("valid")]
    logger.info(
        "REINVENT generated %d / %d valid (prior=%s)", len(valid_only), len(sampled), prior_used
    )

    # 1b. Novelty hard drop — fingerprint-identical (or near-identical) copies
    #     of known ChEMBL actives are not distinct in the provided fingerprint
    #     comparison. This is NOT a patentability or novelty assessment. Only applies in target-aware
    #     mode where we have a reference set to compare against.
    novelty_drop_count = 0
    if target_aware:
        kept = []
        for s in valid_only:
            t = s.get("max_tanimoto")
            if t is not None and float(t) >= NOVELTY_HARD_DROP:
                novelty_drop_count += 1
                continue
            kept.append(s)
        if novelty_drop_count:
            logger.info(
                "Novelty filter: dropped %d literal/near-identical copies (Tanimoto >= %.2f)",
                novelty_drop_count, NOVELTY_HARD_DROP,
            )
        valid_only = kept

    if not valid_only:
        return _empty_response(target, module, total=len(sampled), note="no valid SMILES")

    # 2. ADMET predict
    smiles_list = [s["smiles"] for s in valid_only]
    admet_preds = await predict_admet(smiles_list)
    # Preserve row identity (including duplicate generated molecules). A dict
    # keyed by SMILES loses multiplicity and can misattribute predictions.
    if len(admet_preds) != len(valid_only) or any(
        a.get("smiles") != s["smiles"] for s, a in zip(valid_only, admet_preds)
    ):
        raise ValueError("ADMET response is not in one-to-one generated molecule order")
    logger.info("ADMET predicted for %d", len(admet_preds))

    # 3. Filter + score
    survivors: list[dict] = []
    for s, a in zip(valid_only, admet_preds):
        smi = s["smiles"]
        rdkit_desc = {"mw": s.get("mw"), "logp": s.get("logp"), "qed": s.get("qed"), "tpsa": s.get("tpsa")}
        admet_raw = a.get("raw") or {}
        # Pass SMILES so _check_filters can enforce the on-demand RDKit
        # descriptors (HBD, HBA, rotatable bonds, formal charge, SA) and the
        # PAINS / BRENK structural-alert filters that the original implementation
        # was silently skipping. Returns (hard_failures, soft_warnings) — the
        # warnings include BRENK substructure alerts that pharma reviewers
        # want to see but that we don't auto-reject (BRENK is too aggressive
        # to use as a hard gate; it nukes aspirin / osimertinib / etc.).
        failed_filters, soft_warnings = _check_filters(rdkit_desc, admet_raw, overrides, smiles=smi)
        if failed_filters:
            continue
        survivors.append(
            {
                "smiles": smi,
                "nll": float(s.get("nll", 0.0)),
                "rdkit_desc": rdkit_desc,
                "admet_raw": admet_raw,
                "admet_pct": a.get("percentile") or {},
                "alerts": soft_warnings,  # informational only — BRENK etc.
                "max_tanimoto": s.get("max_tanimoto"),
            }
        )

    # 4a. Pre-rank by ADMET+similarity score, take top dock_top_k for Vina docking.
    #     Docking is expensive (~5-30 sec per ligand), so we never dock the full
    #     surviving population — only the top candidates that already passed
    #     ADMET + filters and look promising on cheap signals.
    pre_formatted = [
        _format_candidate(
            0, x["smiles"], x["nll"], x["rdkit_desc"], x["admet_raw"], x["alerts"], module,
            max_tanimoto=x.get("max_tanimoto"),
            target_aware=target_aware,
            dock_affinity=None,  # filled in step 4b
        )
        for x in survivors
    ]
    pre_formatted.sort(key=lambda c: c["score"], reverse=True)

    # 4b. Run Vina on the top-K (configurable) of pre-ranked survivors
    # The UI sends docking.enabled nested inside filter_overrides. Respect
    # it unless the API caller explicitly provided enable_docking top-level.
    docking_filters = (overrides or {}).get("docking") or {}
    enable_docking = bool(req.get("enable_docking", docking_filters.get("enabled", True)))
    dock_top_k = int(req.get("dock_top_k", 30))
    docking_status = "skipped"
    docking_target_used: str | None = None
    docked_affinities: dict[str, float] = {}
    # Parallel pose map (smiles -> PDB block) populated by the docking
    # call below. Persisted in _format_candidate so the Results modal can
    # render the real Vina pose via /api/jobs/{id}/poses/{rank}.pdb.
    docked_poses: dict[str, str] = {}
    # Same idea but for the multi-MODEL block (all Vina poses concatenated).
    # NGL's `?mode=all` viewer cycles through these. Separate map so we
    # don't bloat the per-row payload for callers who only need top-1.
    docked_poses_all: dict[str, str] = {}

    if enable_docking and pre_formatted:
        # Resolve which Vina target receptor to use. Allow explicit override; else
        # default to the user's target_name if it matches a bundled receptor.
        explicit_dock_target = (req.get("dock_target") or "").strip()
        if explicit_dock_target:
            dock_target_candidate = explicit_dock_target
        else:
            dock_target_candidate = target

        try:
            available = await list_available_targets()
        except Exception:  # noqa: BLE001
            available = []
        # case-insensitive resolve
        match = None
        for t in available:
            if t.lower() == dock_target_candidate.lower():
                match = t
                break

        if match is None:
            docking_status = (
                f"no_receptor_for_{dock_target_candidate!r}; available={available[:8]}"
            )
            logger.info("Docking skipped: %s", docking_status)
        else:
            docking_target_used = match
            top_for_docking = pre_formatted[:dock_top_k]
            top_smiles = [c["smiles"] for c in top_for_docking]
            try:
                dock_results = await run_docking(top_smiles, match)
                for r in dock_results:
                    if r.get("affinity") is not None:
                        docked_affinities[r["smiles"]] = float(r["affinity"])
                    # Parallel pose map keyed on SMILES so _format_candidate
                    # can stamp the right pose onto the right row regardless
                    # of post-sort reordering. Vina may return a pose without
                    # an affinity (rare) — keep it anyway so the NGL viewer
                    # can still render it.
                    pose = r.get("pose_pdb")
                    if isinstance(pose, str) and pose.strip():
                        docked_poses[r["smiles"]] = pose
                    pose_all = r.get("pose_pdb_all")
                    if isinstance(pose_all, str) and pose_all.strip():
                        docked_poses_all[r["smiles"]] = pose_all
                docking_status = (
                    f"ok: docked={len(top_smiles)}, succeeded={len(docked_affinities)},"
                    f" poses={len(docked_poses)}, multi_poses={len(docked_poses_all)}"
                )
                logger.info("Docking done: %s", docking_status)
            except Exception as exc:  # noqa: BLE001
                docking_status = f"failed: {exc}"
                logger.warning("Docking failed: %s", exc)

    # 4c. Re-format with real Vina affinities and re-rank.
    #
    # Two-bucket ranking: docked candidates ALWAYS come first, then the rest.
    # Why: the previous flat re-rank could push docked candidates below
    # un-docked ones whose ChEMBL Tanimoto score was high enough — meaning
    # the user's returned top-N would contain `dock_score: null` rows
    # despite Vina having actually run for that pipeline. The Phase 2
    # agent caught this and surfaced it as "scores might be stored
    # elsewhere or only calculated for a subset of candidates" — correct
    # diagnosis. Fix: bucket-then-sort so the contract "docking ran on
    # the top N you see" actually holds.
    #
    # Within each bucket we still sort by the combined score, so
    # target-relevance (Tanimoto) + ADMET + safety still drive the
    # within-bucket order. Real Vina kcal/mol is the bucket selector,
    # not the rank-key.
    formatted_all = [
        _format_candidate(
            0, x["smiles"], x["nll"], x["rdkit_desc"], x["admet_raw"], x["alerts"], module,
            max_tanimoto=x.get("max_tanimoto"),
            target_aware=target_aware,
            dock_affinity=docked_affinities.get(x["smiles"]),
            dock_pose_pdb=docked_poses.get(x["smiles"]),
            dock_pose_pdb_all=docked_poses_all.get(x["smiles"]),
            docking_for_ranking=not alk5_target,
        )
        for x in survivors
    ]
    docked_bucket = [c for c in formatted_all if c.get("dock_score") is not None]
    undocked_bucket = [c for c in formatted_all if c.get("dock_score") is None]
    docked_bucket.sort(key=lambda c: c["score"], reverse=True)
    undocked_bucket.sort(key=lambda c: c["score"], reverse=True)
    # Guarantee every docked candidate makes the returned set when docking
    # ran — otherwise the agent has no way to see the Vina kcal/mol it just
    # spent ~10s/molecule computing.
    effective_top_n = max(return_top_n, len(docked_bucket)) if docked_bucket else return_top_n
    ranked_all = sorted(formatted_all, key=lambda c: c["score"], reverse=True) if alk5_target else docked_bucket + undocked_bucket
    top = ranked_all[:effective_top_n]
    for i, c in enumerate(top, 1):
        c["rank"] = i

    # Engine label reflects actual capability used. The historical label
    # "reinvent4_admet_vina" implied Vina ran and ChEMBL Tanimoto re-ranked
    # — calling that on a target_aware=False / docked_count=0 run is the
    # exact misrepresentation that allowed garbage inputs ("KLM17XY") to
    # produce confident-looking output. Downgrade the label when neither
    # ChEMBL nor Vina contributed, so audit logs + UI banners can tell
    # the truth about what actually happened.
    vina_ran = len(docked_affinities) > 0
    if target_aware and vina_ran:
        engine_label = "reinvent4_admet_vina"
    elif target_aware:
        engine_label = "reinvent4_admet_chembl"
    elif vina_ran:
        engine_label = "reinvent4_admet_vina_generic"
    else:
        engine_label = "reinvent4_drug_like_only"
    degraded = not (target_aware or vina_ran)
    return {
        "target": target,
        "module": module,
        "total_generated": len(sampled),
        "after_filter_count": len(survivors),
        "top_n": len(top),
        "candidates": top,
        "meta": {
            "engine": "production",
            "discovery_engine": engine_label,
            "run_id": str(uuid4()),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "engine_type": "computational_prediction",
            "is_synthetic": False,
            "research_validated": False,
            "seed": seed,
            "reference_hash": (hashlib.sha256(json.dumps(reference_smiles, sort_keys=True,
                                     separators=(",", ":")).encode("utf-8")).hexdigest()
                               if reference_smiles else None),
            "source_data": ref_source,
            "model_versions": {"generator_prior_identifier": prior_used,
                               "reinvent_version": None, "admet_version": None, "vina_version": None},
            "generator": "REINVENT 4 prior (Apache 2.0, AstraZeneca)",
            "admet": "ADMET-AI (MIT, Stanford / Swanson et al.)",
            "docking": "AutoDock Vina 1.2.7 (Apache 2.0)",
            "prior_used": prior_used,
            "target_relevance_source": ref_source,
            "target_relevance_reference_count": len(reference_smiles),
            "target_aware_scoring": target_aware,
            "ranking_interpretation": "heuristic composite; not an activity probability or safety validation",
            "activity_model_status": "not_deployed",
            "activity_validation_status": "source_generalization_not_demonstrated" if alk5_target else "not_assessed",
            "activity_validation_report": "alk5_source_stress_20261019" if alk5_target else None,
            "docking_used_for_ranking": not alk5_target,
            "alk5_validation_report": "alk5_five_block_validation_20261018" if alk5_target else None,
            "novelty_dropped_count": novelty_drop_count,
            "novelty_threshold": NOVELTY_HARD_DROP,
            "docking_status": docking_status,
            "docking_target": docking_target_used,
            "docking_top_k": dock_top_k if enable_docking else 0,
            "docked_count": len(docked_affinities),
            "filter_overrides": overrides or {},
            "valid_smiles": len(valid_only),
            # Validation summary — frontend reads this to render a big
            # yellow banner when degraded=True. Populated whether or not
            # the user passed force_generic=true (the gate runs first for
            # validation, then we honour the bypass and continue).
            "validation": {
                "degraded": degraded,
                "force_generic": force_generic,
                "chembl_hits": len(reference_smiles),
                "prior_is_fine_tuned": (
                    (prior_used or "").lower() not in {"", "base", "default", "generic"}
                ),
                "vina_ran": vina_ran,
            },
        },
    }


def _empty_response(target: str, module: str, total: int, note: str) -> dict:
    return {
        "target": target,
        "module": module,
        "total_generated": total,
        "after_filter_count": 0,
        "top_n": 0,
        "candidates": [],
        "meta": {
            "engine": "production",
            "discovery_engine": "reinvent4_admet",
            "note": note,
        },
    }
