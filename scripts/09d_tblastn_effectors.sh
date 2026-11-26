#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/mag
: > tblastn_all.tsv
for g in gen/*.fna; do
  acc=$(basename $g .fna)
  tblastn -query fam_members.fa -subject $g -outfmt "6 qseqid sseqid pident length evalue bitscore sstart send" \
          -evalue 1e-5 -max_target_seqs 20 2>/dev/null | awk -v a=$acc '{print a"\t"$0}' >> tblastn_all.tsv
done
echo TBLASTN_DONE $(wc -l < tblastn_all.tsv)
