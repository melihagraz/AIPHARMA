"""Paired ADMET comparison: alignment, missingness and real router contract."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import AsyncMock, patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fastapi import FastAPI
from fastapi.testclient import TestClient
from app.services.alk5_design_admet import compare_profiles, profile_inputs
from app.api.routes_alk5_design import router


class ProfileTests(unittest.TestCase):
    def setUp(self):
        self.design = dict(requested=1, candidates=[dict(candidate_id='design_1', smiles='CCO', provenance=dict(parent_id='parent_1', parent_smiles='CCN'))], meta=dict(reference_asset_sha256='test-asset'))

    def test_alignment_and_missing_values_never_become_zero(self):
        predictions = [dict(smiles='CCN', raw=dict(DILI=0, hERG=.7, Caco2_Wang=-6)), dict(smiles='CCO', raw=dict(hERG=.2, Caco2_Wang=-5))]
        result = compare_profiles(self.design, predictions)
        e = result['pairs'][0]['endpoints']
        self.assertEqual(e['DILI'], dict(parent=0., candidate=None, delta=None, direction='missing'))
        self.assertAlmostEqual(e['hERG']['delta'], -.5)
        self.assertEqual(e['Caco2_Wang']['delta'], 1.)
        self.assertEqual(result['summary']['DILI']['missing'], 1)
        self.assertFalse(result['meta']['safety_validated'])

    def test_invalid_classification_outputs_and_duplicate_rows(self):
        for v in [float('nan'), float('inf'), -1, 2, True, 'not-a-number']:
            r = compare_profiles(self.design, [dict(smiles='CCO', raw=dict(DILI=v)), dict(smiles='CCN', raw=dict(DILI=.4))])
            self.assertIsNone(r['pairs'][0]['endpoints']['DILI']['candidate'])
        with self.assertRaises(ValueError):
            compare_profiles(self.design, [dict(smiles='CCO'), dict(smiles='CCO')])

    def test_missing_parent_and_shared_parent_deduplication(self):
        design = copy.deepcopy(self.design)
        design['candidates'].append(dict(candidate_id='design_2', smiles='CCC', provenance=dict(parent_id='parent_1', parent_smiles='CCN')))
        self.assertEqual(profile_inputs(design), ['CCO', 'CCN', 'CCC'])
        r = compare_profiles(design, [dict(smiles='CCO', raw=dict(DILI=.2)), dict(smiles='CCC', raw=dict(DILI=.3))])
        self.assertEqual(r['meta']['missing_prediction_structures'], 1)
        self.assertEqual(r['summary']['DILI']['missing'], 2)

    def test_real_router_alignment_and_service_error(self):
        app = FastAPI()
        app.include_router(router)
        predictions = [dict(smiles='CCO', raw=dict(DILI=.2)), dict(smiles='CCN', raw=dict(DILI=.5))]
        with TestClient(app) as client, patch('app.api.routes_alk5_design.alk5_analogues', new=AsyncMock(return_value=self.design)), patch('app.api.routes_alk5_design.predict_profiles', new=AsyncMock(return_value=predictions)) as predictor:
            r = client.post('/api/motor/alk5/analogues/admet', json=dict(count=1))
            self.assertEqual(r.status_code, 200)
            self.assertAlmostEqual(r.json()['pairs'][0]['endpoints']['DILI']['delta'], -.3)
            predictor.assert_awaited_once_with(['CCO', 'CCN'])
            predictor.return_value = predictions + [predictions[0]]
            self.assertEqual(client.post('/api/motor/alk5/analogues/admet', json=dict(count=1)).status_code, 502)


if __name__ == '__main__':
    unittest.main()
