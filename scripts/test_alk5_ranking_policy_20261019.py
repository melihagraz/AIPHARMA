"""Dependency-isolated regression tests against actual production function bodies.

External HTTP services and filters are stubbed; this is not an integration test.
"""
import ast
import asyncio
import logging
from pathlib import Path
import unittest
from unittest.mock import AsyncMock


class RankingPolicyTest(unittest.TestCase):
    def setUp(self):
        source = Path(__file__).resolve().parents[1]/'backend/app/services/production_pipeline.py'
        tree = ast.parse(source.read_text(encoding='utf-8-sig'))
        names = {'_normalize_dock_affinity', '_novelty_penalty', '_component_scores', '_final_score', '_format_candidate', 'run_production_discovery'}
        bodies = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in names]
        self.ns = dict(NOVELTY_SOFT_LO=.85, NOVELTY_HARD_DROP=.99, logger=logging.getLogger(__name__))
        exec(compile(ast.Module(body=bodies, type_ignores=[]), str(source), 'exec'), self.ns)

    def test_likelihood_is_not_binding(self):
        score = self.ns['_component_scores']
        self.assertEqual(score({'qed': .6}, {}, 1), score({'qed': .6}, {}, 100))
        self.assertNotIn('binding', score({}, {}, 1))

    def test_alk5_docking_retained_without_score_influence(self):
        fmt = self.ns['_format_candidate']
        a = fmt(1, 'CCO', 1, {'qed': .6}, {}, [], 'general', dock_affinity=-5, docking_for_ranking=False)
        b = fmt(1, 'CCO', 100, {'qed': .6}, {}, [], 'general', dock_affinity=-12, docking_for_ranking=False)
        self.assertEqual(a['score'], b['score'])
        self.assertEqual(b['dock_score'], -12)
        self.assertNotIn('binding', b['component_scores'])

    def test_partial_docking_does_not_promote_docked_alk5_candidate(self):
        ns = self.ns
        ns['_resolve_reference_set'] = AsyncMock(return_value=(['CCN'], 'test'))
        ns['list_available_priors'] = AsyncMock(return_value=['base'])
        ns['list_available_targets'] = AsyncMock(return_value=['ALK5'])
        ns['sample_smiles'] = AsyncMock(return_value=([
            dict(smiles='CCO', valid=True, nll=1, qed=.9, max_tanimoto=.5),
            dict(smiles='CCC', valid=True, nll=100, qed=.2, max_tanimoto=.5)], 'base'))
        ns['predict_admet'] = AsyncMock(return_value=[dict(smiles=s, raw={}) for s in ['CCO', 'CCC']])
        ns['_check_filters'] = lambda *a, **kw: ([], [])
        ns['run_docking'] = AsyncMock(return_value=[dict(smiles='CCC', affinity=-12, pose_pdb='POSE')])
        result = asyncio.run(ns['run_production_discovery'](dict(target_name='ALK5')))
        self.assertEqual(result['candidates'][0]['smiles'], 'CCO')
        self.assertEqual(result['candidates'][1]['dock_score'], -12)
        self.assertEqual(result['candidates'][1]['pose_pdb'], 'POSE')
        self.assertFalse(result['meta']['docking_used_for_ranking'])
        self.assertEqual(result['meta']['activity_model_status'], 'not_deployed')
        self.assertEqual(result['meta']['activity_validation_status'], 'source_generalization_not_demonstrated')


if __name__ == '__main__':
    unittest.main()
