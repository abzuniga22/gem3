#!/usr/bin/env bash
set -euo pipefail

: "${DATA_DIR:=/tmp}"        # quick test: use /tmp; switch to /var/data later for persistence
mkdir -p "$DATA_DIR"
DEST="$DATA_DIR/GEM3.db"

# If we don't already have a DB on this machine/instance, provision one
if [ ! -s "$DEST" ]; then
  echo "Seeding DB into $DEST ..."
  if [ -f "GEM3.db" ] && head -c 16 "GEM3.db" | grep -q "SQLite format 3"; then
    # If a real DB is present (e.g., local dev), just copy it
    cp "GEM3.db" "$DEST"
  elif [ -f "GEM3.db.xz" ]; then
    # Decompress the compressed DB we committed
    DEST="$DEST" python - <<'PY'
import os, lzma, shutil
dst = os.environ["DEST"]
with lzma.open("GEM3.db.xz","rb") as fsrc, open(dst,"wb") as fdst:
    shutil.copyfileobj(fsrc, fdst)
print("Decompressed to:", dst)
PY
  elif [ -n "${DB_URL:-}" ]; then
    # Optional: allow a direct download URL fallback if you later set DB_URL in Render
    curl -L --fail --retry 5 --retry-delay 2 -o "$DEST" "$DB_URL"
  else
    echo "ERROR: No DB source found (no real GEM3.db, no GEM3.db.xz, no DB_URL)"; exit 1
  fi
fi

# sanity check header
if ! head -c 16 "$DEST" | grep -q "SQLite format 3"; then
  echo "ERROR: $DEST is not a valid SQLite file"; exit 1
fi

exec gunicorn --bind "0.0.0.0:${PORT:-10000}" app:app
