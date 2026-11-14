"""Source-document deletion diagnostic for archived ALK5 fixed OOF scores.

This is a *descriptive* sensitivity check: it neither retrains a model nor
simulates a new held-out assay. Documents may overlap; deleting a document
removes every compound citing that document, once per compound.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATAPATH = ROOT / "reports/alk5_endpoint_archive_20261016/dataset.json"
PREDPATH = ROOT / "reports/alk5_endpoint_archive_20261016/predictions.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def evaluate_source_deletions(dataset: list[dict], predictions: list[dict]) -> dict:
    ds = {r["smiles"]: r for r in dataset}
    if len(ds) != len(dataset):
        raise ValueError("Duplicate reference SMILES; source deletion is ambiguous")
    for row in dataset:
        if row["label"] not in (0, 1) or not isinstance(row.get("documents"), list):
            raise ValueError("Invalid labels or document provenance")
    groups = Counter(doc for row in dataset if row["label"] == 0 for doc in set(row["documents"]))
    if not groups:
        raise ValueError("No negative-containing source documents")
    out = {"dataset_n": len(dataset), "negative_n": sum(r["label"] == 0 for r in dataset),
           "negative_containing_documents": len(groups), "splits": {}}
    for split in sorted({r["split"] for r in predictions}):
        methods = {}
        for method in sorted({r["method"] for r in predictions if r["split"] == split}):
            scored = {r["smiles"]: r for r in predictions if r["split"] == split and r["method"] == method}
            if len(scored) != len(dataset) or set(scored) != set(ds):
                raise ValueError(f"Missing/duplicate score rows in {split}/{method}")
            for smiles, r in scored.items():
                if r["label"] != ds[smiles]["label"]:
                    raise ValueError(f"Label mismatch for {smiles}")
            original = list(scored.values())
            original_auc = roc_auc_score([r["label"] for r in original], [r["score"] for r in original])
            original_ap = average_precision_score([r["label"] for r in original], [r["score"] for r in original])
            deletions = []
            for doc, negatives_in_doc in sorted(groups.items()):
                retained = [r for r in original if doc not in ds[r["smiles"]]["documents"]]
                labels = [r["label"] for r in retained]
                if len(set(labels)) < 2:
                    continue
                score = [r["score"] for r in retained]
                auc = float(roc_auc_score(labels, score))
                ap = float(average_precision_score(labels, score))
                deletions.append({"excluded_document": doc, "negative_rows_at_source": negatives_in_doc,
                                  "removed_compounds": len(dataset) - len(retained),
                                  "remaining_active": sum(labels), "remaining_negative": len(labels)-sum(labels),
                                  "roc_auc": auc, "delta_auc": auc-float(original_auc),
                                  "average_precision": ap})
            if not deletions:
                raise ValueError(f"No class-balanced source deletions for {split}/{method}")
            methods[method] = {"original_roc_auc": float(original_auc),
                               "original_average_precision": float(original_ap),
                               "deleted_source_count": len(deletions),
                               "deletion_auc_min": min(row["roc_auc"] for row in deletions),
                               "deletion_auc_max": max(row["roc_auc"] for row in deletions),
                               "deletions": deletions}
        out["splits"][split] = methods
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "reports/phase6/alk5_source_deletion.json")
    args = ap.parse_args()
    data = json.loads(DATAPATH.read_text())
    predictions = json.loads(PREDPATH.read_text())
    report = evaluate_source_deletions(data, predictions)
    report["source_sha256"] = {"dataset": sha256(DATAPATH), "predictions": sha256(PREDPATH)}
    report["interpretation"] = ("Fixed historical score sensitivity only. Deleting all rows linked "
      "to each negative-containing document does not constitute re-training, independent validation, "
      "or confidence interval construction. Documents can overlap.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False)+"\n")
    for split, methods in report["splits"].items():
        for method, item in methods.items():
            print(f"{split} {method}: fixed AUC={item['original_roc_auc']:.4f}, "
                  f"leave-one-negative-source-out [{item['deletion_auc_min']:.4f}, "
                  f"{item['deletion_auc_max']:.4f}]")

if __name__ == "__main__":
    main()
