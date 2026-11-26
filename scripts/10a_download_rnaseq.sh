#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/rna
mkdir -p fq
cat urls.txt | xargs -P 6 -I{} sh -c 'f=fq/$(basename "{}"); [ -s "$f" ] || curl -s -L --max-time 600 -o "$f" "{}"'
echo FASTQ_DONE $(ls fq/*.gz 2>/dev/null | wc -l) files $(du -sh fq | cut -f1)
