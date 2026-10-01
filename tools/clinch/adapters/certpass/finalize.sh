#!/bin/bash
# finalize — mechanical end-of-pass steps: table, manifest, seal.
set -e
CP="$(cd "$(dirname "$0")" && pwd)"
cd "$CP"
python3 summarize_certs.py | tee receipts/CERT_TABLE.md
# sha manifest over receipts + adapter code + data inputs
( find receipts -type f ! -name MANIFEST.sha256 | sort
  ls *.py | sort | sed 's|^|./|'
  ls data/*.json 2>/dev/null | sort ) | while read f; do
    sha256sum "$f"
done > receipts/MANIFEST.sha256
wc -l receipts/MANIFEST.sha256
chmod -R a-w receipts/* || true
chmod u+w receipts 2>/dev/null || true
echo "sealed: $(date -u)"
