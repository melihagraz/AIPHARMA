"""Unit tests for the reward / utility helpers in
Reinvent_App/scripts/rl_finetune.py.

We import the module directly via importlib so we test the actual reward
logic the RL loop uses, not a stub. Full training-loop integration is
covered by the smoke-test invocation that ran 2 steps × batch=4 during
development — that ~12s round-trip is too slow for unit-test cadence.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_RL_PATH = (
    Path(__file__).resolve().parents[2] / "Reinvent_App" / "scripts" / "rl_finetune.py"
)


@pytest.fixture(scope="module")
def rl_module():
    """Load rl_finetune.py as a module. Pytest collects this once per
    session — torch import is ~2 s so we don't repeat it per test."""
    # The script lives in scripts/ which is not on sys.path; we need
    # Reinvent_App on the path so the dynamic `from reinvent_core ...`
    # imports inside the module work.
    rein_root = _RL_PATH.parent.parent
    if str(rein_root) not in sys.path:
        sys.path.insert(0, str(rein_root))
    spec = importlib.util.spec_from_file_location("rl_finetune_under_test", _RL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_qed_score_valid_smiles(rl_module):
    """QED returns a float in [0, 1] for a parseable molecule."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCO")
    score = rl_module._qed_score(mol)
    assert 0.0 <= score <= 1.0
    # Ethanol QED is well-known ≈ 0.41; loose bound guards against
    # a future RDKit version that materially changes the formula.
    assert 0.2 < score < 0.7


def test_pains_pass_returns_1_for_clean_molecule(rl_module):
    """A small drug-like molecule (caffeine) should pass PAINS."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")  # caffeine
    assert rl_module._pains_pass(mol) == 1.0


def test_pains_pass_returns_0_for_pains_alert(rl_module):
    """A SMILES with a known PAINS alert should score 0. Curcumin's
    polyhydroxyphenol motif hits the PAINS_A catalog."""
    from rdkit import Chem

    # Catechol motif — PAINS A hit (catechol_A 1) per the RDKit
    # FilterCatalog reference set. If a future RDKit catalog update
    # removes this hit, this test will need a new sentinel.
    mol = Chem.MolFromSmiles("Oc1ccccc1O")  # catechol
    assert rl_module._pains_pass(mol) == 0.0


def test_max_tanimoto_zero_when_no_references(rl_module):
    """Empty ref list returns the neutral 0.5 — we don't want a
    zero-Tanimoto signal in the reward when there's nothing to compare
    against (would unfairly punish every molecule)."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCO")
    assert rl_module._max_tanimoto(mol, []) == 0.5


def test_max_tanimoto_self_match_is_one(rl_module):
    """Comparing a molecule against its own fingerprint must give 1.0."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCO")
    ref_fps = rl_module._build_reference_fingerprints(["CCO"])
    assert rl_module._max_tanimoto(mol, ref_fps) == pytest.approx(1.0)


def test_compose_reward_arithmetic_mean(rl_module):
    """The composer must use weighted ARITHMETIC mean (not geometric),
    so a single zero component shouldn't collapse the reward to zero.
    Tested by passing a known PAINS-hitter — pains contributes 0 but
    qed + tani should still produce a positive reward."""
    weights = {"qed": 1.0, "pains": 1.0, "tani": 0.0}
    rewards, stats = rl_module._compose_reward(
        ["Oc1ccccc1O"],  # catechol — PAINS hit, but QED is non-zero
        weights,
        ref_fps=[],
    )
    assert len(rewards) == 1
    # PAINS=0 weighted 1/2, QED ≈ 0.5 weighted 1/2 — arithmetic mean
    # produces something strictly > 0, which is the whole point.
    assert rewards[0].item() > 0.0
    assert stats["pains_pass_rate"] == 0.0


def test_compose_reward_invalid_smiles_returns_zero(rl_module):
    """Garbage SMILES contributes 0 reward without raising. The RL loop
    feeds the model's raw sample output here — invalid samples are
    expected (especially early in training) and must not crash."""
    weights = {"qed": 1.0, "pains": 1.0, "tani": 1.0}
    rewards, stats = rl_module._compose_reward(
        ["", "$%^not_a_smiles", "CCO"],
        weights,
        ref_fps=[],
    )
    assert rewards[0].item() == 0.0  # empty
    assert rewards[1].item() == 0.0  # un-parseable
    assert rewards[2].item() > 0.0   # valid ethanol
    # validity counts only the parseable ones.
    assert stats["validity"] == pytest.approx(1 / 3)


def test_compose_reward_validity_rate(rl_module):
    """Stats dict's validity should track parseable / total."""
    rewards, stats = rl_module._compose_reward(
        ["CCO", "CCN", "CCC", "%not%"],
        {"qed": 1.0, "pains": 1.0, "tani": 0.0},
        ref_fps=[],
    )
    assert stats["validity"] == pytest.approx(3 / 4)


def test_pdbqt_to_pdb_helper_strips_autodock_columns():
    """The Docking_App pose stripper lives in dock.py — verify it
    truncates the AutoDock charge + atom-type extension cleanly so the
    output is parseable by standard PDB readers (NGL.js, RDKit).

    We load the helper via importlib instead of `from app.dock import …`
    because the backend already owns the `app` package name and a
    sys.path-based import would collide (Docking_App also has an `app/`
    package). Direct file load → no shared module-namespace at all.
    """
    docking_dock_py = (
        Path(__file__).resolve().parents[2] / "Docking_App" / "app" / "dock.py"
    )
    mod_name = "docking_dock_under_test"
    spec = importlib.util.spec_from_file_location(mod_name, docking_dock_py)
    docking_mod = importlib.util.module_from_spec(spec)
    # Must register in sys.modules BEFORE exec_module so DockResult's
    # @dataclass machinery can resolve its own `__module__` lookup —
    # otherwise Python's dataclass typing code raises AttributeError
    # on `sys.modules.get(cls.__module__).__dict__`.
    sys.modules[mod_name] = docking_mod
    try:
        spec.loader.exec_module(docking_mod)
        _pdbqt_block_to_pdb = docking_mod._pdbqt_block_to_pdb
    finally:
        sys.modules.pop(mod_name, None)

    src = (
        "MODEL 1\n"
        "HETATM    1  C   LIG A   1      12.345  56.789   1.234  1.00  0.00     0.123 C \n"
        "ENDMDL\n"
        "MODEL 2\n"
        "HETATM    1  C   LIG A   1      13.000  57.000   1.500  1.00  0.00     0.123 C \n"
        "ENDMDL\n"
    )
    # only_first_model=True: must drop MODEL 2 entirely.
    out_first = _pdbqt_block_to_pdb(src, only_first_model=True)
    assert "MODEL 2" not in out_first
    assert "12.345" in out_first
    assert "END" in out_first
    # AutoDock extension dropped: the "0.123 C" trailing columns are
    # past column 66 so the truncation strips them.
    for line in out_first.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            assert len(line) <= 66

    # only_first_model=False: both MODELs preserved.
    out_all = _pdbqt_block_to_pdb(src, only_first_model=False)
    assert "MODEL 1" in out_all
    assert "MODEL 2" in out_all
