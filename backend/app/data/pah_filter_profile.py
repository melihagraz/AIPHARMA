"""Canonical ADMET / drug-likeness filter profile for the PAH pipeline run.

This is the single source of truth for the thresholds applied in the manuscript
end-to-end run (ALK5 / BMPR2). It is consumed by the discovery pipeline as
``filter_overrides`` and it is the exact specification the paper's Methods 2.4
must state — code and manuscript stay in lock-step by both pointing here.

Threshold provenance (every number is literature-grounded, not ad hoc):

  base_filters
    mw_max               500    Lipinski et al. 1997 (Rule of Five)
    logp_max             5      Lipinski et al. 1997
    hbd_max              5      Lipinski et al. 1997
    hba_max              10     Lipinski et al. 1997
    rotatable_bonds_max  10     Veber et al. 2002  (NOT 12 — see note)
    tpsa_max             140    Veber et al. 2002  (NOT 150 — see note)
    qed_min              0.50   Bickerton et al. 2012 (quantitative drug-likeness)
    sa_score_max         6.0    Ertl & Schuffenhauer 2009 (synthetic accessibility)

  safety_filters (ADMET-AI predicted probabilities, 0..1; reject if value > max)
    ames_max             0.50   PRIMARY toxicity exclusion — predicted mutagens
                                 dropped (genotoxic liability). Ames is the hard
                                 gate; the rest are secondary.
    herg_risk_max        0.50   cardiotoxicity (hERG channel block)
    dili_risk_max        0.50   drug-induced liver injury

  structural alerts (enforced in production_pipeline._check_filters, not here)
    PAINS                hard reject   (Baell & Holloway 2010)
    BRENK                soft warning  (Brenk et al. 2008) — annotation only

  pk_filters
    (empty) — ADMET-AI clearance / half-life range checks stay OFF pending unit
    normalization, matching the current pipeline behaviour. Do not silently
    enable them; that would change what the paper reports.

NOTE on Veber (important correctness fix): the manuscript draft stated
"rotatable bonds <= 12" and "TPSA <= 150". Those numbers are wrong for the
Veber 2002 reference it cites — the actual Veber criteria are <=10 rotatable
bonds and TPSA <= 140 A^2 (the "12" is Veber's *H-bond count* limit, a
different criterion). This profile uses the correct 10 / 140, which also
matches the legacy server/ implementation.

ADMET engine: predictions come from ADMET-AI (Swanson et al. 2024), a
graph-neural-network model *trained on* the Therapeutics Data Commons (TDC)
benchmark datasets (Huang et al. 2022). The manuscript should cite ADMET-AI as
the engine and TDC as the training data — not TDC single-task models directly.
"""

from __future__ import annotations

from typing import Any

PAH_FILTER_PROFILE: dict[str, dict[str, Any]] = {
    "base_filters": {
        "mw_max": 500.0,
        "logp_max": 5.0,
        "hbd_max": 5,
        "hba_max": 10,
        "rotatable_bonds_max": 10,
        "tpsa_max": 140.0,
        "qed_min": 0.50,
        "sa_score_max": 6.0,
    },
    "safety_filters": {
        "ames_max": 0.50,
        "herg_risk_max": 0.50,
        "dili_risk_max": 0.50,
    },
    "pk_filters": {},
}
"""Reject a molecule if any base/safety threshold is violated; PAINS is an
additional hard structural gate and BRENK a soft annotation."""
