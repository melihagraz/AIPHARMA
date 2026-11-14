"""Archival scientific sensitivity checks; these do not run live ADMET models."""
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path: sys.path.insert(0, str(ROOT))
from scripts.alk5_parent_sensitivity_phase4 import analyze


def test_parent_sensitivity_counts_and_provenance():
    result=analyze()
    assert result['not_experimental_or_new_model_predictions'] is True
    assert result['parent_cluster_count'] == 16
    assert result['candidate_count'] == 30
    assert result['comparisons'] == 180
    assert result['source_sha256'] == '3d9c7e2b90fb94baedcf1247dd3face254303b766b0bdfa889cea1ae26a95665'
    assert result['endpoints']['DILI']['improved_count'] == 13
    assert result['endpoints']['Caco2_Wang']['improved_count'] == 23


def test_parent_deletion_is_descriptive_not_wetlab():
    result=analyze()
    for ep, info in result['endpoints'].items():
        assert len(info['parent_deletions']) == 16
        assert info['leave_one_parent_out_mean_min'] <= info['overall_mean_predicted_delta'] <= info['leave_one_parent_out_mean_max']
        assert 0 <= info['improved_count'] <= 30
        assert all(row['remaining_designs'] in (28,29) for row in info['parent_deletions'])


def test_malformed_source_is_rejected(tmp_path):
    bad=tmp_path/'bad.csv'
    bad.write_text('candidate_id,parent_id,endpoint,parent,candidate,delta\na,b,DILI,0.2,0.6,0.2\n')
    with pytest.raises(ValueError):
        analyze(bad)


def test_paired_auc_archive_reproducible():
    from scripts.alk5_paired_auc_phase4 import analyze
    a=analyze(seed=2, n_boot=75)
    b=analyze(seed=2, n_boot=75)
    assert a==b
    assert a['weak_inactive_count']==11
    assert abs(a['splits']['document_and_scaffold']['logistic_auc'] - 0.9175230566534915) < 1e-10
    assert abs(a['splits']['document_and_scaffold']['similarity_auc'] - 0.7761528326745717) < 1e-10
    assert a['not_retrained_not_external_validation'] is True
