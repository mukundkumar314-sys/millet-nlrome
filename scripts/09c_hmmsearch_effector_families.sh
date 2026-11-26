#!/usr/bin/env bash
# Stage 09c - search five curated blast effector families against the pathogen proteomes.
#
# max.hmm is the concatenation of five Pfam-A profiles:
#   AVR-Pik_HID (PF22383)  MAX_effector (PF23364)  AvrPiz-t (PF22347)
#   AvrPi9 (PF32858)       Avr_Pia (PF33093)
# Each is applied at its own gathering threshold.
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/mag

# one FASTA of every protein, each header prefixed with its assembly accession
if [ ! -s all_prot.fa ]; then
  : > all_prot.fa
  while read -r acc; do
    sed "s/^>/>${acc}|/" "prot/${acc}.faa" >> all_prot.fa
  done < acc.txt
fi

hmmsearch --cut_ga --tblout max_ga.tbl max.hmm all_prot.fa > max_ga.out

python3 - <<'PY'
import collections, json
counts = collections.defaultdict(dict)
for line in open("max_ga.tbl"):
    if line.startswith("#"):
        continue
    f = line.split()
    acc = f[0].split("|", 1)[0]
    counts[f[2]][acc] = counts[f[2]].get(acc, 0) + 1
json.dump(counts, open("max_counts.json", "w"), indent=1)
print("families:", {k: len(v) for k, v in counts.items()})
PY
echo "HMMSEARCH_DONE"
