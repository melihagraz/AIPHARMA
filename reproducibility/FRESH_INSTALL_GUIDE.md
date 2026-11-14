# Clean-install tiers (Linux / WSL2)

## A. Evidence-audit software tests — stdlib only

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv
cd /path/to/AIPHARMA
bash reproducibility/clean_eave_smoke.sh
```

This reproduces **unit-test behavior**, not molecular model outputs. No model weights or data are needed.

## B. Frontend package installation (to test independently)

```bash
cd /path/to/AIPHARMA
node --version  # intended test baseline: Node 22
npm ci
npm run build
npm test         # package.json uses vitest run
```

Use lockfile `package-lock.json`; record whether each command completed. These commands have NOT been executed as a clean fresh-install test by this release kit.

## C. Python scientific install and replay — pending cross-service lock reconciliation

```bash
cd /path/to/AIPHARMA
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
# After reviewing actual platform-specific lock:
python -m pip install -r reproducibility/locks/pip-freeze-py311.txt
python -m pip check
python scripts/phase15_live_replay.py check
python scripts/phase15_live_replay.py run
```

**Prerequisite:** the approved `pip-freeze-py311.txt` must exist. Current REINVENT requirements include `numpy<2` whereas the recorded user's scientific environment may use NumPy 2.x. Reconcile with isolated service environments or tested compatible constraints **before** claiming a single clean install. Vina is a native binary, not delivered by `pip freeze`. Model weights, datasets and receptor/grid preparation need their own URLs, checksums and licenses. If absent, report a blocker rather than treating a unit test as scientific reproduction.
