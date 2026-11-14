#!/usr/bin/env bash
# Clone ONLY version-controlled content; verify the actual fresh-clone test path.
# Usage: bash clean_clone_smoke.sh ~/AIPHARMA [~/aipharma_cleanclone_test]
set -uo pipefail
SOURCE="${1:-$HOME/AIPHARMA}"
DEST="${2:-$HOME/aipharma_cleanclone_test}"
if [ -e "$DEST" ]; then echo "REFUSING: destination exists: $DEST" >&2; exit 2; fi
if ! command -v git >/dev/null; then echo "BLOCKED: git missing";exit 2;fi
if ! git -C "$SOURCE" rev-parse --is-inside-work-tree >/dev/null 2>&1;then echo 'BLOCKED: source is not a Git checkout';exit 2;fi
mkdir -p "$(dirname "$DEST")"
git clone --no-local -- "$SOURCE" "$DEST" || exit 2
GIT_SHA="$(git -C "$DEST" rev-parse HEAD)"
mkdir -p "$DEST/release_check_logs"
printf '%s\n' "$GIT_SHA" >"$DEST/release_check_logs/commit.txt"
printf '%s\n' "SOURCE=$SOURCE" "DEST=$DEST" "SHA=$GIT_SHA" >"$DEST/release_check_logs/protocol.txt"
if ! command -v docker >/dev/null; then
  echo 'BLOCKED: docker not available; clean-clone test not executed'
  echo 'BLOCKED: docker unavailable' > "$DEST/release_check_logs/status.txt"
  exit 3
fi
if ! docker info >/dev/null 2>&1; then
  echo 'BLOCKED: Docker daemon unavailable'
  echo 'BLOCKED: Docker daemon unavailable' > "$DEST/release_check_logs/status.txt"
  exit 3
fi
(
 cd "$DEST" || exit 3
 set -o pipefail
 docker build --progress=plain -f backend/Dockerfile.test -t "aipharma-repro:${GIT_SHA:0:12}" . \
  2>&1 | tee release_check_logs/build.log || exit 4
 docker run --rm "aipharma-repro:${GIT_SHA:0:12}" -q --tb=short -m 'not slow' \
  2>&1 | tee release_check_logs/pytest.log || exit 5
)
CODE=$?
if [ "$CODE" -eq 0 ];then
 echo 'PASS: clean local-Git-clone Docker backend test. Not full scientific E2E and not public URL availability.' | tee "$DEST/release_check_logs/status.txt"
else
 echo "FAIL: clean-clone test exit=$CODE; see logs" | tee "$DEST/release_check_logs/status.txt"
fi
exit "$CODE"
