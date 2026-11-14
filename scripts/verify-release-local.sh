#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: bash scripts/verify-release-local.sh /path/to/assets.tar.gz"
  exit 1
fi

ARCHIVE="$(realpath "$1")"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

echo "=== PYTHON VERSION ==="
python3 -c 'import sys; assert sys.version_info[:2] == (3,12), "Python 3.12 required"'

echo "=== VERIFY ASSET ARCHIVE ==="
printf '%s  %s\n' \
  '350be962b406d2b03fd5999ae2145f4772d50027ab2c69c8fa5c6353e07a0500' \
  "$ARCHIVE" | sha256sum -c -

echo "=== EXTRACT VERIFIED ASSETS ==="
tar -xzf "$ARCHIVE"

echo "=== VERIFY ASSET FILES ==="
python3 - <<'PY'
import hashlib
import json
from pathlib import Path

manifest = json.loads(Path("assets_manifest.json").read_text())

for item in manifest["files"]:
    p = Path(item["path"])
    assert p.is_file(), f"Missing: {p}"
    assert hashlib.sha256(p.read_bytes()).hexdigest() == item["sha256"], f"Hash mismatch: {p}"

print("Verified:", len(manifest["files"]), "assets")
PY

echo "=== CREATE PYTHON ENVIRONMENT ==="
if [ -e .venv_release ]; then
  echo "Existing .venv_release found; use a fresh clone."
  exit 1
fi
python3 -m venv .venv_release
PYTHON="$ROOT/.venv_release/bin/python"

echo "=== INSTALL CPU PYTORCH ==="
"$PYTHON" -m pip install \
  torch==2.14.0 \
  --index-url https://download.pytorch.org/whl/cpu

echo "=== INSTALL VERIFIED DEPENDENCIES ==="
TEMP_REQ="$(mktemp)"
trap 'rm -f "$TEMP_REQ"' EXIT

grep -v '^torch==' requirements-release-verified.txt > "$TEMP_REQ"
"$PYTHON" -m pip install -r "$TEMP_REQ"

echo "=== DEPENDENCY VALIDATION ==="
"$PYTHON" -m pip check

echo "=== BACKEND TESTS ==="
PYTHONPATH="$ROOT" "$PYTHON" -m pytest \
  backend/tests -q --tb=short -m "not slow"

echo "=== VERIFIED GIT COMMIT ==="
git rev-parse HEAD

echo "=== RELEASE VERIFICATION COMPLETE ==="
