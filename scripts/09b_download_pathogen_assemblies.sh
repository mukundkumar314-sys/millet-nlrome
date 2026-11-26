#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/mag
while read acc; do
  [ -s gen/$acc.fna ] && continue
  curl -s -L --max-time 300 -o /tmp/g_$acc.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" \
  && unzip -o -q -j /tmp/g_$acc.zip "*_genomic.fna" -d /tmp/gx_$acc 2>/dev/null \
  && f=$(ls /tmp/gx_$acc/*.fna 2>/dev/null | head -1) && [ -n "$f" ] && mv "$f" gen/$acc.fna
  rm -rf /tmp/g_$acc.zip /tmp/gx_$acc
done < gacc.txt
echo GENOMES_DONE $(ls gen/*.fna 2>/dev/null | wc -l)
