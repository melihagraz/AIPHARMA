"""Remove legacy hardcoded fake pharmacophore_output AuditLog rows.

Background
----------
Until 2026-05-11 the ``/jobs/start-pharmacophore`` endpoint generated
pharmacophore feature counts via the placeholder formula

    acceptor    = max(1, len(residues) // 3)
    donor       = max(1, len(residues) // 3)
    hydrophobic = max(1, len(residues) // 2)
    positive    = 1   # hardcoded
    negative    = 1   # hardcoded

— i.e. acceptor count ALWAYS equals donor count, and positive/negative
are always exactly 1. No actual RDKit pharmacophore extraction was
performed. The Analysis Results table showed those numbers, looking
deceptively scientific.

The endpoint has now been rewritten to (a) tag manual-residue runs with
``"source": "manual_residue_hint"`` and an explicit note, and (b) accept
real Pocket_App results via the new ``/jobs/log-pharmacophore-analysis``
endpoint. This script ages out the old placeholder rows so the Analysis
Results table reflects only honest data.

Detection fingerprint (conservative — ALL must match before deletion):
  - acceptor == donor (the // 3 ratio gives identical numbers)
  - positive == 1 AND negative == 1 (hardcoded)
  - source field is missing (new rows include it)

Run
---
    cd D:\\ilac-master\\ilac-master
    .venv\\Scripts\\python.exe scripts\\cleanup_fake_pharmacophores.py            # dry-run
    .venv\\Scripts\\python.exe scripts\\cleanup_fake_pharmacophores.py --apply    # actually delete
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BACKEND_DIR = ROOT / "backend"
sys.path.insert(0, str(BACKEND_DIR))

from app.db.session import SessionLocal  # noqa: E402
from app.models.audit_log import AuditLog  # noqa: E402


def _looks_fake(payload: dict) -> bool:
    """True when the payload matches the legacy placeholder fingerprint."""
    # New rows include a "source" tag; legacy ones don't.
    if "source" in payload:
        return False
    # Top-level numeric fields the placeholder always wrote.
    acceptor = payload.get("acceptor")
    donor = payload.get("donor")
    positive = payload.get("positive")
    negative = payload.get("negative")
    if not all(isinstance(v, (int, float)) for v in (acceptor, donor, positive, negative)):
        return False
    # The smoking gun: acceptor==donor and positive==negative==1.
    return (
        int(acceptor) == int(donor)
        and int(positive) == 1
        and int(negative) == 1
    )


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
            .filter(AuditLog.resource_type == "job_output", AuditLog.action == "pharmacophore_output")
            .all()
        )
        if not rows:
            print("No pharmacophore_output AuditLog rows present — nothing to clean up.")
            return 0

        fake_rows: list[AuditLog] = []
        for row in rows:
            try:
                payload = json.loads(row.payload or "{}")
            except Exception:
                payload = {}
            if _looks_fake(payload):
                fake_rows.append(row)

        print(f"Scanned {len(rows)} pharmacophore_output rows.")
        print(f"Found {len(fake_rows)} legacy placeholder rows.")
        if not fake_rows:
            return 0

        for row in fake_rows:
            try:
                payload = json.loads(row.payload or "{}")
                pdb = payload.get("pdb_id", "?")
                a = payload.get("acceptor", "?")
                d = payload.get("donor", "?")
                h = payload.get("hydrophobic", "?")
            except Exception:
                pdb, a, d, h = "?", "?", "?", "?"
            print(f"  - audit_id={row.id}  job_id={row.resource_id}  pdb={pdb}  A/D/H={a}/{d}/{h}")

        if not args.apply:
            print()
            print("Dry-run only. Re-run with --apply to delete these rows.")
            return 0

        for row in fake_rows:
            db.delete(row)
        db.commit()
        print(f"\nDeleted {len(fake_rows)} legacy placeholder pharmacophore rows.")
        return 0
    finally:
        db.close()


if __name__ == "__main__":
    raise SystemExit(main())
