#!/usr/bin/env python3
"""Leave-one-parent-out sensitivity of archived ALK5 ADMET predictions.

No activity measurements, prospective validation, or new ADMET inference.
This descriptive perturbation probes whether mean predicted shifts rely on a
single parent scaffold/family; it is not a confidence or p-value calculation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'reports/alk5_design_admet_20261019/paired_outputs.csv'
RESULT = ROOT / 'reports/phase4/alk5_parent_sensitivity.json'
FIGURE = ROOT / 'reports/phase4/alk5_direction_summary.png'
ENDPOINTS = {
    'DILI': -1,
    'hERG': -1,
    'AMES': -1,
    'Bioavailability_Ma': 1,
    'Caco2_Wang': 1,
    'Solubility_AqSolDB': 1,
}

def analyze(path: Path = SOURCE) -> dict:
    source = path.read_bytes()
    df = pd.read_csv(path)
    required = {'candidate_id','parent_id','endpoint','parent','candidate','delta'}
    if not required.issubset(df.columns):
        raise ValueError(f'Missing columns: {sorted(required - set(df.columns))}')
    if df.duplicated(['candidate_id', 'endpoint']).any():
        raise ValueError('Duplicate candidate-endpoint pair')
    if not np.isfinite(df[['parent', 'candidate', 'delta']].to_numpy(dtype=float)).all():
        raise ValueError('Nonfinite archive predictions')
    if not np.allclose(df.candidate - df.parent, df.delta, atol=1e-6):
        raise ValueError('Archived differences not internally consistent')
    if set(df.endpoint) != set(ENDPOINTS):
        raise ValueError('Unexpected endpoint identities')
    if df.candidate_id.nunique() != 30 or df.parent_id.nunique() != 16 or len(df) != 180:
        raise ValueError('Archive does not match 30 candidates, 16 parents and 180 comparisons')

    result = {
        'analysis':'descriptive_leave_one_parent_out_archived_predicted_admet',
        'not_experimental_or_new_model_predictions':True,
        'source_relative_path':str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
        'source_sha256':hashlib.sha256(source).hexdigest(),
        'candidate_count':30,'parent_cluster_count':16, 'comparisons':180,
        'interpretation':'Mean shifts and leave-one-parent-out ranges are descriptive sensitivity, not inferential CIs. Endpoints use distinct, non-comparable numerical scales.',
        'endpoints':{},
    }
    for ep, desired in ENDPOINTS.items():
        subset = df[df.endpoint.eq(ep)]
        parent_ids = sorted(subset.parent_id.unique())
        loo = []
        for parent_id in parent_ids:
            remaining = subset[subset.parent_id.ne(parent_id)]
            loo.append({'removed_parent_id':parent_id,
                        'remaining_designs':len(remaining),
                        'mean_delta':float(remaining.delta.mean())})
        overall=float(subset.delta.mean())
        fraction=int((desired*subset.delta > 0).sum())
        result['endpoints'][ep] = {
            'desired_direction':'lower' if desired < 0 else 'higher',
            'overall_mean_predicted_delta':overall,
            'improved_count':fraction,
            'total_count':len(subset),
            'leave_one_parent_out_mean_min':min(x['mean_delta'] for x in loo),
            'leave_one_parent_out_mean_max':max(x['mean_delta'] for x in loo),
            'leave_one_parent_out_sign_consistent':all(
                x['mean_delta']*overall > 0 for x in loo
            ) if overall != 0 else False,
            'parent_deletions':loo,
        }
    return result


def make_plot(data: dict, path: Path = FIGURE):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    keys=list(ENDPOINTS)
    vals=[data['endpoints'][ep]['improved_count'] for ep in keys]
    labels=['DILI risk', 'hERG risk', 'AMES risk', 'Bioavailability', 'Caco-2', 'Solubility']
    fig, ax=plt.subplots(figsize=(7.4,3.6))
    bars=ax.barh(labels, [i/30*100 for i in vals])
    ax.axvline(50, linestyle='--', linewidth=1, color='0.45', label='50% of designs')
    for bar, n in zip(bars,vals):
        ax.text(bar.get_width()+1,bar.get_y()+bar.get_height()/2, f'{n}/30',va='center',fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0,107)
    ax.set_xlabel('Designs with a numerically favorable predicted shift (%)')
    ax.set_title('Archived ALK5 parent-to-design ADMET changes')
    ax.legend(loc='lower right', fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=180)
    plt.close(fig)


def main():
    data=analyze()
    RESULT.parent.mkdir(parents=True,exist_ok=True)
    RESULT.write_text(json.dumps(data,indent=2)+'\n', encoding='utf-8')
    make_plot(data)
    print(json.dumps({'result':str(RESULT),'figure':str(FIGURE),
                      'endpoint_summary': {k:{x:v for x,v in info.items() if x not in {'parent_deletions'}} for k,info in data['endpoints'].items()}},indent=2))

if __name__ == '__main__':main()
