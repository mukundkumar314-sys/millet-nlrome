# millet-nlrome

Code for **"Genome-wide identification and comparative analysis of NLR immune receptors and
their integrated domains in eleven grass genomes, with implications for blast resistance in
finger millet (*Eleusine coracana* L. Gaertn.)"**

Mukund Kumar Thakur and Sweta Mishra
Department of Genetics and Plant Breeding, Dr. Rajendra Prasad Central Agricultural
University, Pusa, Samastipur, Bihar 848125, India

---

## What this repository contains

Every analysis in the paper is in-silico and starts from public data, so the whole study can
be rerun from the scripts here plus the accessions listed in Supplementary Table S1 (plant
genomes), Table S25 (*Pyricularia oryzae* isolates) and BioProject PRJNA648385 (RNA-seq).

```
scripts/     numbered pipeline stages, run in order
figures/     one script per figure, plus the shared style module
docs/        a longer walk through the pipeline
```

## Quick start

```bash
export PROJECT_ROOT=/path/to/your/working/directory   # every script reads this
pip install -r requirements.txt
python scripts/08_hma_host_targets.py --sweep
```

`PROJECT_ROOT` holds the downloaded data and intermediates, kept out of the repository
because the genomes and read files run to tens of gigabytes. Its layout is described in
`docs/pipeline.md`.

## Pipeline stages

| Stage | Script | Produces |
|---|---|---|
| 01 | `01_build_nlr_set.py` | NB-ARC search, NLR candidate set |
| 02 | `02_classify.py`, `02b_adddomains.py` | CN/CNL/NL architectures, integrated domains |
| 03 | `03_genecoords.py` | gene coordinates from the GFFs |
| 04 | `04_cluster.py` | physical clusters, head-to-head pairs |
| 05-06 | `05_proso.py`, `06_merge_proso.py` | proso millet added to the set |
| 07 | `07_kaks.py` | Ka/Ks by orthogroup occupancy class |
| 08 | `08_hma_host_targets.py` | **Table S21** - standalone HMA proteins, Clade A host targets |
| 09a-09e | `09a`...`09e` | **Table S22** - effector families by host lineage, annotation gap |
| 10a-10c | `10a`...`10c` | **Table S23** - candidate response to blast inoculation |
| 11a-11b | `11a_run_codeml.sh`, `11b_parse_beb.py` | **Table S24** - codon site models, BEB sites |
| - | `tables.py` | assembles the supplementary workbook |

External tools expected on `PATH`: `hmmsearch` (HMMER 3.3.2), `diamond` (2.1.8),
`mafft` (7.505), `iqtree2` (2.2.2.6), `tblastn` (BLAST+ 2.14.0), `salmon` (1.10.2),
`codeml` (PAML 4.9j).

## Reproducibility

Stages 08, 09c, 09e and 10c were written up from an interactive analysis session into
standalone form. Each was then rerun against the original intermediates and reproduces its
published table exactly:

| Script | Check | Result |
|---|---|---|
| `08_hma_host_targets.py` | all 11 rows of Table S21 | identical |
| `09e_summarise_effectors.py` | all 15 rows of Table S22 | identical |
| `09e_summarise_effectors.py` | AVR-Pik annotation gap | 33/33 in assemblies, 32/33 in proteomes, 95.6-100% identity |
| `10c_candidate_induction.py` | all 24 rows of Table S23 | identical |

Two conventions matter enough to state here, because both change the numbers:

* **HMA hits are kept on a domain bit score of 8.0, below the Pfam PF00403 gathering
  threshold of 24.80.** That threshold excludes both validated AVR-Pik targets: OsHIPP19
  scores 11.9 bits and OsHIPP20 scores 14.2. Gene counts vary by at most one across cutoffs
  from 0 to 10, which `--sweep` reports.
* **HMA proteins are counted per gene, never per protein**, since proteomes differ in how
  many isoforms they annotate. The gene comes from the GFF, except where the protein
  identifier carries its own isoform suffix, in which case stripping that suffix defines the
  gene.

## Data

| Source | Where |
|---|---|
| Plant genomes, annotations, proteomes | accessions in Online Resource 1, Table S1 |
| *Pyricularia oryzae* proteomes and assemblies | accessions in Online Resource 1, Table S25 |
| Blast time-course RNA-seq | ENA, BioProject PRJNA648385 |
| Derived tables S1-S25 | Online Resource 1 of the paper |
| Newick trees behind Figs. 9 and 10 | Online Resource 4 of the paper |

## Citation

See `CITATION.cff`. Please cite the paper rather than this repository alone.

## Licence

MIT, see `LICENSE`. The data analysed here belong to their original depositors and carry
their own terms.
