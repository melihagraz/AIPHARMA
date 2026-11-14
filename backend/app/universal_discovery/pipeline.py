from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import random
from copy import deepcopy

from .config import AppConfig
from .models import Candidate, DiscoveryRequest, DiscoveryResponse


@dataclass
class Molecule:
    smiles: str
    descriptors: dict[str, float]
    admet: dict[str, float]
    safety: dict[str, float]
    dock_score: float


def _unit_float(key: str) -> float:
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return int(digest[:12], 16) / float(16**12 - 1)


def _make_smiles(seed: str | None, idx: int, context_salt: str = "") -> str:
    atoms = ["C", "N", "O", "S", "F", "Cl"]
    stable_seed = int.from_bytes(hashlib.sha256(
        json.dumps([seed or "seed", idx, context_salt], separators=(",", ":")).encode("utf-8")
    ).digest()[:8], "big")
    rng = random.Random(stable_seed)
    length = rng.randint(12, 32)
    chunks: list[str] = []
    for _ in range(length):
        token = atoms[rng.randint(0, len(atoms) - 1)]
        if rng.random() < 0.1:
            token += str(rng.randint(1, 9))
        chunks.append(token)
    return "".join(chunks)


def _predict_descriptors(smiles: str) -> dict[str, float]:
    u = _unit_float(smiles)
    return {
        "qed": 0.35 + 0.6 * u,
        "sa_score": 2.0 + 6.0 * _unit_float(smiles + "|sa"),
        "mw": 180 + 420 * _unit_float(smiles + "|mw"),
        "logp": 0.8 + 5.0 * _unit_float(smiles + "|logp"),
        "tpsa": 20 + 180 * _unit_float(smiles + "|tpsa"),
        "rotatable_bonds": int(2 + 12 * _unit_float(smiles + "|rot")),
        "hbd": int(0 + 7 * _unit_float(smiles + "|hbd")),
        "hba": int(1 + 12 * _unit_float(smiles + "|hba")),
        "formal_charge_abs": int(0 + 3 * _unit_float(smiles + "|chg")),
        "cns_mpo": 2.0 + 3.5 * _unit_float(smiles + "|cns_mpo"),
        "logd74": 0.5 + 4.5 * _unit_float(smiles + "|logd74"),
    }


def _predict_safety(smiles: str) -> dict[str, float]:
    return {
        "herg_risk": _unit_float(smiles + "|herg"),
        "dili_risk": _unit_float(smiles + "|dili"),
        "bsep_risk": _unit_float(smiles + "|bsep"),
        "cyp3a4_inhib": _unit_float(smiles + "|cyp3a4"),
        "cyp2d6_inhib": _unit_float(smiles + "|cyp2d6"),
        "cyp2c9_inhib": _unit_float(smiles + "|cyp2c9"),
        "mitotox": _unit_float(smiles + "|mito"),
        "qt_risk": _unit_float(smiles + "|qt"),
    }


def _predict_admet(smiles: str) -> dict[str, float]:
    return {
        "caco2": _unit_float(smiles + "|caco2"),
        "oral_bioavailability": _unit_float(smiles + "|f"),
        "clearance": _unit_float(smiles + "|cl"),
        "ppb": _unit_float(smiles + "|ppb"),
        "half_life_h": 0.5 + 30.0 * _unit_float(smiles + "|t12"),
        "solubility": _unit_float(smiles + "|sol"),
        "cell_permeability": _unit_float(smiles + "|perm"),
        "lung_affinity": _unit_float(smiles + "|lung"),
        "bbb_score": _unit_float(smiles + "|bbb"),
    }


def _dock_score(target_name: str, pdb_id: str | None, smiles: str, context_salt: str = "") -> float:
    key = f"{target_name}|{pdb_id or 'none'}|{smiles}|dock|{context_salt}"
    return -4.5 - 6.5 * _unit_float(key)


def _pharmacophore_context_salt(ctx: dict | None) -> str:
    if not ctx:
        return ""
    try:
        blob = json.dumps(ctx, sort_keys=True, ensure_ascii=True)
    except (TypeError, ValueError):
        blob = str(ctx)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:24]


def _build_runtime_config(cfg: AppConfig, request: DiscoveryRequest) -> dict:
    runtime = {
        "base_filters": deepcopy(cfg.base_filters),
        "safety_filters": deepcopy(cfg.safety_filters),
        "pk_filters": deepcopy(cfg.pk_filters),
        "docking": deepcopy(cfg.docking),
        "scoring": deepcopy(cfg.scoring),
        "modules": deepcopy(cfg.modules),
    }
    overrides = request.filter_overrides or {}
    allowed_sections = {"base_filters", "safety_filters", "pk_filters", "docking"}
    for section, section_overrides in overrides.items():
        if section not in allowed_sections or not isinstance(section_overrides, dict):
            continue
        target = runtime[section]
        for key, value in section_overrides.items():
            if key in target:
                target[key] = value
    return runtime


def _passes_filters(m: Molecule, runtime_cfg: dict, module_name: str) -> tuple[bool, list[str]]:
    bf = runtime_cfg["base_filters"]
    sf = runtime_cfg["safety_filters"]
    pf = runtime_cfg["pk_filters"]
    module = runtime_cfg["modules"].get(module_name, runtime_cfg["modules"]["general"])
    alerts: list[str] = []

    d, s, a = m.descriptors, m.safety, m.admet
    checks = [
        ("qed", d["qed"] >= bf["qed_min"]),
        ("sa_score", d["sa_score"] <= bf["sa_score_max"]),
        ("mw", bf["mw_min"] <= d["mw"] <= bf["mw_max"]),
        ("logp", bf["logp_min"] <= d["logp"] <= bf["logp_max"]),
        ("tpsa", bf["tpsa_min"] <= d["tpsa"] <= bf["tpsa_max"]),
        ("rotatable_bonds", d["rotatable_bonds"] <= bf["rotatable_bonds_max"]),
        ("hbd", d["hbd"] <= bf["hbd_max"]),
        ("hba", d["hba"] <= bf["hba_max"]),
        ("formal_charge_abs", d["formal_charge_abs"] <= bf["formal_charge_abs_max"]),
        ("herg_risk", s["herg_risk"] <= sf["herg_risk_max"]),
        ("dili_risk", s["dili_risk"] <= sf["dili_risk_max"]),
        ("bsep_risk", s["bsep_risk"] <= sf["bsep_risk_max"]),
        ("cyp3a4", s["cyp3a4_inhib"] <= sf["cyp3a4_inhib_max"]),
        ("cyp2d6", s["cyp2d6_inhib"] <= sf["cyp2d6_inhib_max"]),
        ("cyp2c9", s["cyp2c9_inhib"] <= sf["cyp2c9_inhib_max"]),
        ("caco2", a["caco2"] >= pf["caco2_min"]),
        ("oral_bioavailability", a["oral_bioavailability"] >= pf["oral_bioavailability_min"]),
        ("clearance", pf["clearance_min"] <= a["clearance"] <= pf["clearance_max"]),
        ("ppb", a["ppb"] <= pf["ppb_max"]),
        ("half_life_h", pf["half_life_min_hours"] <= a["half_life_h"] <= pf["half_life_max_hours"]),
    ]

    if runtime_cfg["docking"]["enabled"]:
        checks.append(("dock_score", m.dock_score <= runtime_cfg["docking"]["score_threshold"]))

    if module_name == "pah":
        checks += [
            ("pah_mw", module["mw_min"] <= d["mw"] <= module["mw_max"]),
            ("pah_logd74", module["logd74_min"] <= d["logd74"] <= module["logd74_max"]),
            ("pah_lung_affinity", a["lung_affinity"] >= module["lung_affinity_min"]),
            ("pah_mitotox", s["mitotox"] <= module["mitotox_max"]),
            ("pah_qt_risk", s["qt_risk"] <= module["qt_risk_max"]),
            ("pah_bbb_negative", a["bbb_score"] <= 0.45 if module["bbb_negative"] else True),
        ]
    elif module_name == "cns":
        checks += [
            ("cns_tpsa", d["tpsa"] <= module["tpsa_max"]),
            ("cns_mw", d["mw"] <= module["mw_max"]),
            ("cns_mpo", d["cns_mpo"] >= module["cns_mpo_min"]),
            ("cns_bbb_positive", a["bbb_score"] >= 0.55 if module["bbb_positive"] else True),
        ]
    elif module_name == "oncology":
        checks.append(("oncology_perm", a["cell_permeability"] >= module["cell_permeability_min"]))
    elif module_name == "antiviral":
        checks += [
            ("antiviral_solubility", a["solubility"] >= module["solubility_min"]),
            ("antiviral_clearance", a["clearance"] <= module["clearance_max"]),
        ]

    failed = [name for name, ok in checks if not ok]
    alerts.extend(failed)
    return len(failed) == 0, alerts


def _score(m: Molecule, runtime_cfg: dict, module_name: str) -> tuple[float, dict[str, float]]:
    sc = runtime_cfg["scoring"]
    module_extra = float(runtime_cfg["modules"].get(module_name, runtime_cfg["modules"]["general"]).get("extra_weight", 1.0))

    binding = max(0.0, min(1.0, (-m.dock_score - 4.0) / 6.0))
    admet = (m.admet["caco2"] + m.admet["oral_bioavailability"] + (1.0 - m.admet["ppb"])) / 3.0
    safety = 1.0 - (
        m.safety["herg_risk"] + m.safety["dili_risk"] + m.safety["bsep_risk"] + m.safety["cyp3a4_inhib"]
    ) / 4.0
    druglikeness = (m.descriptors["qed"] + max(0.0, 1.0 - m.descriptors["sa_score"] / 10.0)) / 2.0
    novelty = _unit_float(m.smiles + "|novel")

    final_score = (
        sc["binding_weight"] * binding
        + sc["admet_weight"] * admet
        + sc["safety_weight"] * safety
        + sc["druglikeness_weight"] * druglikeness
        + sc["novelty_weight"] * novelty
    ) * module_extra

    components = {
        "binding": round(binding, 4),
        "admet": round(admet, 4),
        "safety": round(safety, 4),
        "druglikeness": round(druglikeness, 4),
        "novelty": round(novelty, 4),
    }
    return round(final_score, 4), components


def run_discovery(request: DiscoveryRequest, cfg: AppConfig) -> DiscoveryResponse:
    module_name = request.disease_module if request.disease_module in cfg.modules else "general"
    runtime_cfg = _build_runtime_config(cfg, request)
    ctx_salt = _pharmacophore_context_salt(request.pharmacophore_context)

    generated: list[Molecule] = []
    for i in range(request.generate_count):
        smiles = _make_smiles(request.smiles_seed, i, ctx_salt)
        generated.append(
            Molecule(
                smiles=smiles,
                descriptors=_predict_descriptors(smiles),
                admet=_predict_admet(smiles),
                safety=_predict_safety(smiles),
                dock_score=round(_dock_score(request.target_name, request.pdb_id, smiles, ctx_salt), 3),
            )
        )

    passed: list[Candidate] = []
    for m in generated:
        ok, alerts = _passes_filters(m, runtime_cfg, module_name)
        if not ok:
            continue
        score, component_scores = _score(m, runtime_cfg, module_name)
        passed.append(
            Candidate(
                rank=0,
                smiles=m.smiles,
                module=module_name,
                score=score,
                dock_score=m.dock_score,
                alerts=alerts,
                component_scores=component_scores,
                descriptors={k: round(float(v), 4) for k, v in m.descriptors.items()},
            )
        )

    passed.sort(key=lambda x: x.score, reverse=True)
    top_n = min(cfg.output_top_n, len(passed))
    top = passed[:top_n]
    for idx, c in enumerate(top, start=1):
        c.rank = idx

    return DiscoveryResponse(
        target=request.target_name,
        module=module_name,
        total_generated=len(generated),
        after_filter_count=len(passed),
        top_n=top_n,
        candidates=top,
        meta={
            "pdb_id": request.pdb_id,
            "desired_route": request.desired_route,
            "pipeline_steps": cfg.raw["pipeline"]["steps"],
            "filter_overrides": request.filter_overrides or {},
            "pharmacophore_context": request.pharmacophore_context or {},
        },
    )
