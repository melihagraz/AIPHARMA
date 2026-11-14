#!/usr/bin/env python3
"""Parent-cluster bootstrap of ARCHIVED model-predicted ALK5 ADMET deltas.

These are deterministic post hoc summaries of existing predictions, not new
ADMET model runs or experimental molecule measurements. A parent-level cluster
bootstrap partly accounts for multiple designs derived from the same parent;
it does NOT capture model or assay uncertainty, nor prove safety improvements.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "reports/alk5_design_admet_20261019/paired_outputs.csv"
OUTPUT = ROOT / "reports/validation_20261109_extended/alk5_cluster_bootstrap_20261109.json"
DIRECTION = {"DILI": -1, "hERG": -1, "AMES": -1,
             "Bioavailability_Ma": 1, "Caco2_Wang": 1,
             "Solubility_AqSolDB": 1}


def analyze(input_path: Path = INPUT, n_boot: int = 5000, seed: int = 20261109):
    raw_bytes = input_path.read_bytes()
    df = pd.read_csv(input_path)
    required = {"candidate_id", "parent_id", "endpoint", "parent", "candidate", "delta"}
    if not required.issubset(df.columns):
        raise ValueError(f"missing columns {required - set(df.columns)}")
    if df.duplicated(["candidate_id", "endpoint"]).any():
        raise ValueError("duplicate design-endpoint pairs")
    if not np.isfinite(df[["parent", "candidate", "delta"]].to_numpy(float)).all():
        raise ValueError("non-finite prediction values")
    if not np.allclose(df["candidate"] - df["parent"], df["delta"], atol=1e-6):
        raise ValueError("archived delta does not match candidate minus parent")
    if set(df.endpoint) != set(DIRECTION):
        raise ValueError(f"unexpected endpoints {set(df.endpoint) ^ set(DIRECTION)}")
    designs = df.candidate_id.nunique()
    parents = df.parent_id.nunique()
    if designs != 30 or parents != 16 or len(df) != 180:
        raise ValueError(f"unexpected archive size: {designs} designs, {parents} parents, {len(df)} comparisons")
    rng = np.random.default_rng(seed)
    output = {
        "analysis": "parent_cluster_resampling_of_archived_predicted_alk5_admet_deltas",
        "input": str(input_path.relative_to(ROOT)),
        "input_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "n_designs": designs, "n_distinct_parent_compounds": parents,
        "n_endpoints": len(DIRECTION), "n_comparisons": len(df),
        "bootstrap": {"clusters": "parent_id", "n_resamples": n_boot, "seed": seed,
                      "interval": "percentile_2.5_to_97.5_percent", "unit": "model_output_scale"},
        "interpretation_limitations": [
            "Predictions are retrospective ADMET-AI outputs, not laboratory measurements.",
            "Resampling parent compounds does not quantify model training, assay, or domain-shift uncertainty.",
            "All models and endpoints are exploratory; no prespecified confirmatory inference.",
            "Predicted ADMET changes do not demonstrate compound safety, activity, or clinical efficacy.",
        ], "endpoints": {},
    }
    for endpoint, sign in DIRECTION.items():
        x = df.loc[df.endpoint == endpoint].copy()
        grouped = {str(parent): chunk.delta.to_numpy(float)
                   for parent, chunk in x.groupby("parent_id", sort=True)}
        labels = sorted(grouped)
        boots = np.empty(n_boot, dtype=float)
        for i in range(n_boot):
            sampled = rng.choice(labels, size=len(labels), replace=True)
            values = np.concatenate([grouped[label] for label in sampled])
            boots[i] = np.mean(values)
        ci = np.quantile(boots, [0.025, 0.975])
        output["endpoints"][endpoint] = {
            "direction_of_desired_change": "decrease" if sign < 0 else "increase",
            "mean_predicted_delta": float(x.delta.mean()),
            "median_predicted_delta": float(x.delta.median()),
            "improved_direction_count": int((sign*x.delta>0).sum()),
            "worsened_direction_count": int((sign*x.delta<0).sum()),
            "zero_delta_count": int((x.delta==0).sum()),
            "cluster_bootstrap_mean_delta_95pct_interval": [float(ci[0]), float(ci[1])],
            "fraction_designs_in_desired_direction": float((sign*x.delta>0).mean()),
        }
    return output


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,default=INPUT)
    p.add_argument("--output",type=Path,default=OUTPUT)
    p.add_argument("--resamples",type=int,default=5000)
    args=p.parse_args()
    result=analyze(args.input,n_boot=args.resamples)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    for name, v in result["endpoints"].items():
        print(name, 'mean_delta',f'{v["mean_predicted_delta"]:+.4f}',
              'CI', [round(a,4) for a in v['cluster_bootstrap_mean_delta_95pct_interval']],
              'improved',v['improved_direction_count'],'/ 30')
    print('saved:',args.output)

if __name__ == '__main__':
    main()
