# Pipeline walkthrough

Every script reads `PROJECT_ROOT` from the environment and defaults to `/home/claude`.
Set it once per shell:

```bash
export PROJECT_ROOT=/path/to/working/directory
```

## Directory layout under `$PROJECT_ROOT`

```
genomes/   plant assemblies
gff/       plant annotations, one gzipped GFF per genome
prot/      plant proteomes
nlr/       NB-ARC searches, the NLR catalogue, orthogroups
hma/       HMA domain tables, small HMA sets, Clade A trees
mag/       Pyricularia oryzae proteomes and assemblies, effector searches
rna/       salmon index, quantifications, experimental design, target sets
paml/      codeml runs, one directory per orthogroup and model
figures/   rendered figures and FACTS_LEDGER.txt
```

## 1. Receptor catalogue (stages 01 to 07)

NB-ARC (PF00931) is searched against each proteome, candidates are classified into CN, CNL,
NL and N architectures, and a complete Pfam-A scan adds integrated domains. Gene coordinates
come from the GFFs, which then support physical clustering and the detection of head-to-head
pairs within 10 kb. Orthogroups are built across the eleven genomes and Ka/Ks is computed
within each.

Proso millet entered the set later, which is why stages 05 and 06 exist separately.

## 2. Host targets of AVR-Pik (stage 08)

```bash
hmmsearch --domtblout hma/<genome>.hma_sens.domtbl -E 10 hma/HMA.hmm prot/<genome>.faa
python scripts/08_hma_host_targets.py --sweep
```

Two choices drive the counts.

The Pfam PF00403 gathering threshold is 24.80 bits. Both validated AVR-Pik targets fall
below it, OsHIPP19 at 11.9 bits and OsHIPP20 at 14.2, so applying the threshold would
discard the very proteins the analysis is about. Hits are kept at a domain bit score of 8.0
instead, and `--sweep` shows gene counts varying by at most one across cutoffs from 0 to 10.

Counting is per gene. Proteomes differ in how many isoforms they annotate, so a per-protein
count is not comparable between genomes. The gene of a protein comes from the GFF, except
where the protein identifier carries its own isoform suffix, in which case stripping that
suffix defines the gene. The script reproduces the published counts for all eleven genomes
under this rule.

Proteins of at most 250 residues are retained as small HMA proteins. Clade A is the smallest
clade containing both reference proteins in the MAFFT and IQ-TREE 2 phylogeny.

## 3. Pathogen effectors (stages 09a to 09e)

```bash
bash scripts/09a_download_pathogen_proteomes.sh
bash scripts/09b_download_pathogen_assemblies.sh
bash scripts/09c_hmmsearch_effector_families.sh
bash scripts/09d_tblastn_effectors.sh
python scripts/09e_summarise_effectors.py
```

Sixty-two annotated proteomes are assigned to a host lineage from the BioSample host field.
The reference isolate 70-15 carries no host field and is assigned to the *Oryza* lineage,
which brings that group to four.

Stage 09d matters for the result. Gene prediction frequently overlooks small effector genes:
in the *Eleusine*-infecting lineage the canonical AVR-Pik family gene appears in one
annotated proteome yet is present in the assembly of all 33, at 95.6 to 100 per cent
amino-acid identity. Presence is therefore scored from the assembly, never from the
annotation alone.

## 4. Expression (stages 10a to 10c)

```bash
bash scripts/10a_download_rnaseq.sh     # verifies every file with gzip -t
bash scripts/10b_salmon_quant.sh
python scripts/10c_candidate_induction.py
```

BioProject PRJNA648385 provides one library per genotype, timepoint and treatment, so
genotype serves as the replicate unit and the contrast rests on twelve paired
inoculated-against-mock comparisons. Testing is confined to 24 pre-specified genes, which is
why no genome-wide differential expression claim appears anywhere in the paper.

Downloaded FASTQ files have arrived truncated before, six of fifty-six in one run, which
drove five samples to near-zero mapping. Stage 10a verifies each file with `gzip -t` and
re-downloads failures.

## 5. Positive selection (stages 11a and 11b)

```bash
bash scripts/11a_run_codeml.sh
python scripts/11b_parse_beb.py
```

Site models are applied only to the four orthogroups holding cloned rice blast resistance
genes, which keeps the multiple-testing burden small. `NSsites` must be `1 2 7 8`, giving
the nested pairs M1a against M2a and M7 against M8; an earlier run used `2 3 10 11` and was
discarded. Sites are taken from the Bayes empirical Bayes posterior under M8 at 0.95 and
mapped onto the finger millet member of each orthogroup.

## 6. Figures and tables

```bash
python figures/fig1.py          # and so on through fig_model.py
python scripts/tables.py        # assembles the supplementary workbook
```

`figures/figstyle.py` holds the shared style: a colourblind-safe palette, single sans-serif
family, and the two standard widths. Figures are written as both 600 dpi PNG and vector PDF.

## Verifying a change

`figures/FACTS_LEDGER.txt` records the command behind every published number. After changing
a stage, rerun it and diff its output against the corresponding supplementary table:

```bash
python scripts/08_hma_host_targets.py --out /tmp/S21.csv
diff <(sort /tmp/S21.csv) <(sort "$PROJECT_ROOT"/hma/TableS21_sHMA_cladeA.csv)
```

A difference is a finding, not a nuisance. Three errors reached the draft and were caught
this way: an isolate count wrong in two places, a mislabelled read-sum column, and a
standalone HMA count carried over from an earlier ten-genome run.
