#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/rna
python3 -c "
import json
d=json.load(open('design.json'))
print('\n'.join(x['run'] for x in d))
" > runlist.txt
while read r; do
  [ -s quant/$r/quant.sf ] && continue
  salmon quant -i ec_idx -l A -1 fq/${r}_1.fastq.gz -2 fq/${r}_2.fastq.gz \
    -p 2 --validateMappings -o quant/$r --quiet >/dev/null 2>&1
done < runlist.txt
echo QUANT_DONE $(ls quant/*/quant.sf 2>/dev/null | wc -l)
