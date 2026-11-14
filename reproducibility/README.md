# AIPHARMA research-software reproducibility release candidate

**State:** `0.1.0-rc1`, LOCAL RESEARCH RELEASE CANDIDATE, NOT PUBLIC, NOT OPEN SOURCE YET.

The evidence pipeline has three distinct scopes:

1. **Code-level tests** (offline/CI): Python unit tests for EAVE, gate and replay logic. Synthetic fixtures are not docking validation.
2. **Local scientific replay** (WSL/Linux): inspection of frozen EGFR/ALK5 scores and one molecule's true ALK5 Vina replay. Requires installed scientific software and local source/input data.
3. **Independent reproduction from scratch:** requires exact environment lock, native Vina binary, target structure/grid files, full data acquisition and model weights. **Not demonstrated** by this release candidate.

## Prerequisites

- Ubuntu/WSL2 recommended. Work done in Ubuntu 24.04 under Python virtual environment; inspect local `python --version` and `node --version` rather than assuming versions.
- Python 3.11 is the conservative clean-environment baseline; Node.js 22 for frontend. Native AutoDock Vina, Meeko, RDKit, ADMET-AI 2.0.1, and REINVENT checkpoint/configs are separate scientific dependencies.
- Use the actual tested, user-generated pip freeze before making reproducibility claims. The pre-existing source requirements are not one consistent exact scientific lock across services (e.g., `Reinvent_App/requirements.txt` specifies `numpy<2`).
- `npm ci` uses the tracked npm lockfile. Do **not** claim a successful clean frontend installation until tested.

## Reproducibility commands (existing local installation)

```bash
cd ~/AIPHARMA
source .venv/bin/activate
python scripts/aipharma_release.py inspect
python scripts/aipharma_release.py freeze
python scripts/aipharma_release.py test
bash reproducibility/clean_eave_smoke.sh
python scripts/aipharma_release.py bundle
```

`inspect` reports publication blockers honestly. `freeze` writes `reproducibility/local/environment.json` and `pip-freeze-pyXY.txt` for **local review**; do not publish until checking for private URLs and local paths. After review, copy the approved freeze into `reproducibility/locks/` (example below). `bundle` makes a conservative **source-only** ZIP in `dist_release/` plus SHA-256 digest and source file manifest. Dataset files, model weights, `.env`, `.venv`, reports and potentially sensitive files are excluded.

```bash
mkdir -p reproducibility/locks
cp reproducibility/local/pip-freeze-py311.txt reproducibility/locks/
# ONLY if pip-freeze-py311.txt was generated in Python 3.11 and manually reviewed.
# If Python version differs, use actual filename and document the difference.
```

The source candidate is **not** a standalone full scientific replication bundle. The source manifest enumerates included files and SHA-256 hashes. Review its missing assets before public release.

## Tests and worked examples

```bash
# Real-data EAVE audit (no new docking; depends on existing reports and source data):
python scripts/eave_audit.py --all
python scripts/phase14_eave_gate.py verify --all

# Fast code tests (synthetic fixtures):
python -m unittest discover -s tests -p 'test_eave_engine.py' -q
python -m unittest discover -s tests -p 'test_phase14_gate.py' -q
python -m unittest discover -s tests -p 'test_phase15_live_replay.py' -q
python -m unittest discover -s tests -p 'test_release_packaging.py' -q

# Real docking replay — opt in only; uses prior Phase11 data:
python scripts/phase15_live_replay.py check
python scripts/phase15_live_replay.py run

# Crystal redocking, where external 3D structure preparation is configured:
PYTHONPATH="$PWD" python scripts/phase9_redock.py --target ALK5 --exhaustiveness 16
```

To reproduce **the benchmark AUCs from frozen scores**, retain source `.ism` and panel file hashes, reactor/grid parameters and the original `analysis.json`, not just a summary value. For a truly independent run, repeat frozen panel generation and all 100 dockings per target/receptor on the same pinned environment and log the result.

## Release procedure

1. Confirm distribution rights and license; inspect third-party bundled components and file notices. Root `LICENSE` is intentionally absent until approved.
2. Verify no passwords, API tokens, credentials, copyrighted/unlicensed datasets, patient data, vendor checkpoints, and private URLs enter the repository.
3. Capture actual WSL environment and native binary/target hashes; resolve dependency conflicts and exercise a fresh virtual environment.
4. Execute unit tests, scientific audits, single real docking replay, and at least one clean install in a separate machine/container.
5. Create a public Git repository only after clearance. Commit reviewed files; tag `v0.1.0-rc1` only on the corresponding commit. Push with your own credentials. Verify hash of published release ZIP.
6. Create a GitHub release and optionally Zenodo archive; fill the **actual**, not anticipated, URL/DOI into `CITATION.cff`. Deposit input data or external download instructions according to dataset licensing and attribution requirements.
7. Report commit SHA, release SHA-256, runtime versions, exact data/weights digests and failures in manuscript. Never replace pending work with `verified`.

See `DATA_AVAILABILITY.md`, `THIRD_PARTY_NOTICES.md`, and `docs/MANUSCRIPT_REPRODUCIBILITY_EN.md`.
