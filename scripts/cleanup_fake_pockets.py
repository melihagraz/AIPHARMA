"""Remove legacy hardcoded fake pocket_output AuditLog rows.

Background
----------
Until 2026-05-11, ``JobService._ensure_pocket_output`` injected a hardcoded
2-pocket placeholder (``druggability_score`` 0.84 / 0.66, fixed residues
``ASP25 / GLU54 / TYR82`` and ``SER98 / LYS102 / PHE110``) into the
``AuditLog`` table whenever a Job auto-completed. Those rows poisoned the
UI's "Completed Analyses" dropdown with synthetic data.

That generator has now been deleted. Real pocket data is written only when
the synchronous Pocket_App pipeline completes (via
``POST /api/jobs/log-pocket-analysis``). This script removes the orphaned
fake rows so the dropdown shows a clean slate.

Run
---
    cd D:\\ilac-master\\ilac-master
    .venv\\Scripts\\python.exe scripts\\cleanup_fake_pockets.py            # dry-run
    .venv\\Scripts\\python.exe scripts\\cleanup_fake_pockets.py --apply    # actually delete

The dry-run prints exactly which rows would be deleted, with their job id,
PDB id, and pocket counts, so you can review before committing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Make the backend package importable when running this script from the repo root.
ROOT = Path(__file__).resolve().parents[1]
BACKEND_DIR = ROOT / "backend"
sys.path.insert(0, str(BACKEND_DIR))

from app.db.session import SessionLocal  # noqa: E402
from app.models.audit_log import AuditLog  # noqa: E402

# The exact fingerprint of the hardcoded generator's output.
FAKE_DRUGGABILITY = {0.84, 0.66}
FAKE_RESIDUE_SUFFIXES = {"-ASP25", "-GLU54", "-TYR82", "-SER98", "-LYS102", "-PHE110"}


def _looks_fake(payload: dict) -> bool:
    """True when the payload matches the legacy hardcoded fingerprint."""
    pockets = payload.get("pockets") or []
    if not isinstance(pockets, list) or not pockets:
        return False
    seen_scores = set()
    seen_suffixes = set()
    for p in pockets:
        if not isinstance(p, dict):
            continue
        score = p.get("druggability_score")
        if isinstance(score, (int, float)):
            seen_scores.add(round(float(score), 2))
        for r in p.get("residues") or []:
            if isinstance(r, str):
                for suffix in FAKE_RESIDUE_SUFFIXES:
                    if r.endswith(suffix):
                        seen_suffixes.add(suffix)
    # Be conservative: require BOTH the score fingerprint AND at least one
    # residue suffix from the hardcoded list. That way a real analysis that
    # happens to score 0.84 won't be incorrectly flagged.
    return bool(seen_scores & FAKE_DRUGGABILITY) and bool(seen_suffixes)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually delete the rows. Without this flag, only a dry-run is performed.",
    )
    args = parser.parse_args()

    db = SessionLocal()
    try:
        rows = (
            db.query(AuditLog)
            .filter(AuditLog.resource_type == "job_output", AuditLog.action == "pocket_output")
            .all()
        )
        if not rows:
            print("No pocket_output AuditLog rows present — nothing to clean up.")
            return 0

        fake_rows: list[AuditLog] = []
        for row in rows:
            try:
                payload = json.loads(row.payload or "{}")
            except Exception:
                payload = {}
            if _looks_fake(payload):
                fake_rows.append(row)

        print(f"Scanned {len(rows)} pocket_output rows.")
        print(f"Found {len(fake_rows)} fake hardcoded rows.")
        if not fake_rows:
            return 0

        for row in fake_rows:
            try:
                payload = json.loads(row.payload or "{}")
                pdb = payload.get("pdb_id", "—")
                n = len(payload.get("pockets") or [])
            except Exception:
                pdb, n = "?", "?"
            print(f"  - audit_id={row.id}  job_id={row.resource_id}  pdb={pdb}  pockets={n}")

        if not args.apply:
            print()
            print("Dry-run only. Re-run with --apply to delete these rows.")
            return 0

        for row in fake_rows:
            db.delete(row)
        db.commit()
        print(f"\nDeleted {len(fake_rows)} fake pocket_output rows.")
        return 0
    finally:
        db.close()


if __name__ == "__main__":
    raise SystemExit(main())
