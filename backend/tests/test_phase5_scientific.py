"""Check Phase 5 offline audit arithmetic. No claims of new measurements."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from alk5_negative_influence_phase5 import analyze, ranking_share_of_one_negative


def test_pairwise_negative_contribution_counts_ties_as_half():
    assert ranking_share_of_one_negative(np.asarray([0.,1.,2.]),1.)==pytest.approx(0.5)


def test_fixed_archived_splits_reconcile_prior_reports():
    result=analyze()
    assert result['counts']['total']==356
    assert result['counts']['weak_or_inactive']==11
    assert result['counts']['distinct_scaffolds']==98
    assert result['results']['document_and_scaffold/logistic']['roc_auc']==pytest.approx(0.9175230566534915)
    assert result['results']['document_and_scaffold/similarity']['roc_auc']==pytest.approx(0.7761528326745717)
    for metrics in result['results'].values():
        assert len(metrics['weak_compound_contributions'])==11
        assert metrics['negative_leave_one_out_auc_range'][0]<=metrics['roc_auc']+0.1
        assert metrics['average_precision_above_active_only_baseline']<=1-result['counts']['active_fraction']+1e-12
