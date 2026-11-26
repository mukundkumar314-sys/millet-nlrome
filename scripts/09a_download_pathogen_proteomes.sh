#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/mag
ok=0; fail=0
while read acc; do
  [ -s prot/$acc.faa ] && { ok=$((ok+1)); continue; }
  url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$acc/download?include_annotation_type=PROT_FASTA"
  if curl -s -L --max-time 180 -o /tmp/$acc.zip "$url" && unzip -o -q -j /tmp/$acc.zip "*protein.faa" -d /tmp/x_$acc 2>/dev/null; then
     f=$(ls /tmp/x_$acc/*.faa 2>/dev/null | head -1)
     if [ -n "$f" ] && [ -s "$f" ]; then mv "$f" prot/$acc.faa; ok=$((ok+1)); else fail=$((fail+1)); echo "NOPROT $acc"; fi
  else fail=$((fail+1)); echo "DLFAIL $acc"; fi
  rm -rf /tmp/$acc.zip /tmp/x_$acc
done < acc.txt
echo "DONE ok=$ok fail=$fail"
