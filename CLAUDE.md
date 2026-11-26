# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

Analysis code for a comparative NLR genomics paper across eleven grass genomes, written up
for *Physiology and Molecular Biology of Plants*. Every result is in-silico and derived from
public data. The repository holds code only: genomes, read files and intermediates live
under `$PROJECT_ROOT` and are deliberately untracked.

## The rule that matters most

**Never invent a number.** Every figure, table and sentence in the manuscript traces to a
command that was actually run. If a value cannot be recomputed from the data under
`$PROJECT_ROOT`, say so plainly rather than estimating, interpolating or carrying a
plausible number forward from memory. When a stored number and a fresh computation
disagree, investigate the discrepancy and report it; do not quietly adopt whichever looks
tidier. Three real errors were caught exactly this way: a transcript-count mismatch, a
mislabelled read-sum column, and an isolate count that was wrong in two places.

## Conventions this project already follows

* **Paths.** Every script reads `PROJECT_ROOT` from the environment, defaulting to
  `/home/claude`. Keep new scripts to that pattern; never hardcode an absolute path.
* **Counting.** HMA and NLR counts are per gene, not per protein. Proteomes differ in
  isoform annotation, so a per-protein count is not comparable across genomes.
* **Thresholds.** Where a Pfam gathering threshold would exclude a validated positive
  control, a lower bit cutoff is used and the sensitivity of the result to that choice is
  reported. Do not silently restore a gathering threshold.
* **Verification.** A stage that regenerates a published table should be run against the
  original intermediates and diffed against that table before it is considered done.
* **Manuscript prose.** The authors prefer positive phrasing and keep the word "not" to a
  minimum; titles are conventional and descriptive rather than claim-shaped.

## Running a stage

```bash
export PROJECT_ROOT=/path/to/data
python scripts/08_hma_host_targets.py --out /tmp/S21.csv --sweep
```

Then diff against the published table before trusting the output.

## Layout under $PROJECT_ROOT

```
genomes/ gff/ prot/   downloaded plant data
nlr/                  NB-ARC searches, NLR catalogue, orthogroups
hma/                  HMA domain tables, small HMA sets, Clade A trees
mag/                  Pyricularia oryzae proteomes, assemblies, effector searches
rna/                  salmon index, quantifications, design and target sets
paml/                 codeml runs, one directory per orthogroup and model
figures/              rendered figures and the facts ledger
```

`figures/FACTS_LEDGER.txt` records, for each published number, the command that produced it.
Consult it before changing any value that appears in the manuscript, and append to it when a
new number enters the paper.

## External tools

`hmmsearch`, `diamond`, `mafft`, `iqtree2`, `tblastn`, `salmon`, `codeml`. Check with
`which` before assuming one is installed; several stages fail slowly otherwise.

## What to be careful about

* Stages 09a, 09b and 10a download tens of gigabytes. Check for existing files first; both
  download scripts already skip what is present.
* Downloaded FASTQ files have arrived truncated before. `gzip -t` every file before
  quantifying; stage 10a does this.
* `codeml` runs take hours. Confirm `NSsites` is `1 2 7 8` before starting a run; an earlier
  run used the wrong values and had to be discarded.
