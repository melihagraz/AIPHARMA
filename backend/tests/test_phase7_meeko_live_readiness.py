"""Phase 7: deterministic ligand preparation and Meeko public API contract.

These are isolated adapter tests using a stub Meeko interface. They are NOT
substitutes for an actual Meeko + Vina integration run.
"""
import sys
from types import ModuleType

import numpy as np
from rdkit import Chem

from Docking_App.app.dock import smiles_to_pdbqt
from Docking_App.app.schemas import DockResponse


def _stub_meeko(monkeypatch, num_setups=1, ok=True, pdbqt="ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00\n"):
    fake = ModuleType('meeko')
    observed = []

    class MoleculePreparation:
        def prepare(self, mol):
            observed.append(np.array(mol.GetConformer().GetPositions(), copy=True))
            return [object() for _ in range(num_setups)]

    class PDBQTWriterLegacy:
        @staticmethod
        def write_string(setup):
            return pdbqt, ok, "failed" if not ok else None

    fake.MoleculePreparation = MoleculePreparation
    fake.PDBQTWriterLegacy = PDBQTWriterLegacy
    monkeypatch.setitem(sys.modules, "meeko", fake)
    return observed


def test_ligand_prep_uses_public_meeko_setup_and_stable_rdkit_embedding(monkeypatch):
    poses = _stub_meeko(monkeypatch)
    a = smiles_to_pdbqt("CCO")
    b = smiles_to_pdbqt("CCO")
    assert a == b and a.startswith('ATOM')
    assert len(poses) == 2
    assert np.array_equal(poses[0], poses[1]), "same canonical SMILES should yield same starting conformer"


def test_ligand_prep_refuses_multiple_molecule_setups(monkeypatch):
    _stub_meeko(monkeypatch, num_setups=2)
    assert smiles_to_pdbqt("CCO") is None


def test_ligand_prep_refuses_empty_or_invalid_pdbqt(monkeypatch):
    _stub_meeko(monkeypatch, pdbqt="REMARK not an atom\n")
    assert smiles_to_pdbqt("CCO") is None


def test_ligand_prep_refuses_meeko_error(monkeypatch):
    _stub_meeko(monkeypatch, ok=False)
    assert smiles_to_pdbqt("CCO") is None


def test_ligand_prep_refuses_invalid_smiles(monkeypatch):
    _stub_meeko(monkeypatch)
    assert smiles_to_pdbqt("not_a_smiles") is None


def test_vina_api_does_not_falsely_claim_unmeasured_binary_version():
    assert "1.2.7" not in DockResponse.model_fields['engine'].default


def test_local_e2e_production_service_credential(monkeypatch):
    from scripts.run_true_three_stage_e2e_20261109 import _local_service_headers
    monkeypatch.setenv('APP_ENV', 'production')
    monkeypatch.setenv('INTERNAL_SERVICE_KEY', 'k' * 40)
    assert _local_service_headers() == {'X-Internal-Service-Key': 'k' * 40}


def test_local_e2e_fails_closed_without_production_credentials(monkeypatch):
    import pytest
    from scripts.run_true_three_stage_e2e_20261109 import _local_service_headers
    monkeypatch.setenv('APP_ENV', 'production')
    monkeypatch.delenv('INTERNAL_SERVICE_KEY', raising=False)
    with pytest.raises(RuntimeError, match='INTERNAL_SERVICE_KEY'):
        _local_service_headers()
    monkeypatch.setenv('APP_ENV', 'unknown')
    with pytest.raises(RuntimeError, match='APP_ENV'):
        _local_service_headers()
