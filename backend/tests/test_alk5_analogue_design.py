"""RDKit design branch and real FastAPI router tests; no external services."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fastapi import FastAPI
from fastapi.testclient import TestClient
from rdkit import Chem, DataStructs
from rdkit.Chem import BRICS, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from app.services.alk5_analogue_design import ASSET, design_analogues
from app.api.routes_alk5_design import router


class DesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = design_analogues(30)
        cls.asset = json.loads(ASSET.read_text())

    def test_structure_provenance_and_bounds(self):
        rows = self.result['candidates']
        self.assertGreater(len(rows), 0)
        self.assertGreater(len({r['provenance']['replacement_fragment'] for r in rows}), 1)
        self.assertEqual(len({r['normalized_parent'] for r in rows}), len(rows))
        known = set(self.asset['excluded_normalized_parents'])
        refs = {r['id']: r for r in self.asset['references']}
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        for row in rows:
            m = Chem.MolFromSmiles(row['smiles'])
            p = row['provenance']
            parent = Chem.MolFromSmiles(p['parent_smiles'])
            self.assertFalse(any(a.GetAtomicNum() == 0 for a in m.GetAtoms()))
            self.assertEqual(Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)), Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(parent)))
            self.assertNotIn(row['normalized_parent'], known)
            self.assertGreaterEqual(row['parent_similarity'], .5)
            self.assertLessEqual(row['parent_similarity'], .9)
            self.assertAlmostEqual(row['parent_similarity'], DataStructs.TanimotoSimilarity(gen.GetFingerprint(m), gen.GetFingerprint(parent)))
            self.assertEqual(p['parent_activity_ids'], refs[p['parent_id']]['activity_ids'])
            self.assertEqual(p['parent_documents'], refs[p['parent_id']]['documents'])
            self.assertTrue(set(p['donor_ids']) <= refs.keys())
            for donor_id in p['donor_ids']:
                donor = Chem.MolFromSmiles(refs[donor_id]['smiles'])
                donor_parts = {Chem.MolToSmiles(part) for bond in BRICS.FindBRICSBonds(donor) for part in Chem.GetMolFrags(BRICS.BreakBRICSBonds(donor, bonds=[bond]), asMols=True)}
                self.assertIn(p['replacement_fragment'], donor_parts)
            self.assertEqual(row['activity_status'], 'unmeasured')
            self.assertEqual(row['safety_status'], 'not_evaluated')
            self.assertNotIn('nll', row)
            self.assertNotIn('score', row)
            self.assertNotIn('activity_probability', row)

    def test_determinism_and_prefix(self):
        smaller = design_analogues(5)
        self.assertEqual(smaller['candidates'], self.result['candidates'][:5])

    def test_invalid_counts_and_corrupt_reference_fail_explicitly(self):
        for count in [0, 101, True, 2.5]:
            with self.assertRaises(ValueError):
                design_analogues(count)
        bad = copy.deepcopy(self.asset)
        bad['references'][0]['smiles'] = 'not_a_structure'
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'invalid.json'
            path.write_text(json.dumps(bad))
            with self.assertRaises(ValueError):
                design_analogues(1, asset_path=path)

    def test_router_serializes_real_design_and_rejects_bad_input(self):
        app = FastAPI()
        app.include_router(router)
        with TestClient(app) as client:
            response = client.post('/api/motor/alk5/analogues', json={'count': 5})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response.json()['candidates'], self.result['candidates'][:5])
            self.assertFalse(response.json()['meta']['activity_validated'])
            for payload in [{'count': 0}, {'count': 101}, {'count': True}, {'count': '5'}, {'target': 'EGFR'}]:
                self.assertEqual(client.post('/api/motor/alk5/analogues', json=payload).status_code, 422)


if __name__ == '__main__':
    unittest.main()
