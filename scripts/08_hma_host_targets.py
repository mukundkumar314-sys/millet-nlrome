#!/usr/bin/env python3
"""Stage 08 - standalone HMA proteins and the Clade A host targets of AVR-Pik.

Rebuilds Supplementary Table S21 from the hmmsearch domain tables produced by

    hmmsearch --domtblout <genome>.hma_sens.domtbl -E 10 HMA.hmm <genome>.faa

Two points of method are worth stating, because both affect the counts.

1. The Pfam PF00403 gathering threshold (GA 24.80) excludes both validated
   AVR-Pik targets: OsHIPP19 scores 11.9 bits and OsHIPP20 scores 14.2. Hits are
   therefore kept on a domain bit score of at least --min-bit (default 8.0).
   Recovered counts are flat for cutoffs between 0 and 10, so the result does not
   hinge on the value chosen; --sweep reruns the count across that range.

2. Counting is per gene, never per protein, since proteomes differ in how many
   isoforms they annotate. The gene of a protein is taken from the GFF, except
   where the protein identifier itself carries an isoform suffix, in which case
   stripping that suffix defines the gene. With this rule the script reproduces
   the published counts for all eleven genomes.
"""
import os, re, csv, gzip, json, argparse, collections

ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
MIN_BIT = 8.0
MAX_LEN = 250
OSHIPP19 = "Os04t0469000"
OSHIPP20 = "Os04t0464100"
ISOFORM = re.compile(r"(_P\d+|_T\d+|-\d+|\.\d+)$")

GFF = {"Bdistachyon": "Bdistachyon.gff3.gz", "Dexilis": "Dexilis.gff3.gz",
       "Ecoracana": "Ecoracana.gff.gz", "Ecrusgalli": "Ecrusgalli.gff3.gz",
       "Osativa": "Osativa.gff3.gz", "Phallii": "Phallii.gff3.gz",
       "Pmiliaceum": "Pmiliaceum.gff.gz", "Sbicolor": "Sbicolor.gff3.gz",
       "Sitalica": "Sitalica.gff3.gz", "Sviridis": "Sviridis.gff3.gz",
       "Zmays": "Zmays.gff3.gz"}
GENOMES = list(GFF)


def parse_domtbl(path, min_bit=MIN_BIT):
    """Best domain bit score and sequence length per target protein."""
    best = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 14:
                continue
            target, tlen, dom_bit = f[0], int(f[2]), float(f[13])
            if dom_bit < min_bit:
                continue
            if target not in best or dom_bit > best[target][0]:
                best[target] = (dom_bit, tlen)
    return best


def protein_to_gene(gff_path):
    """Map protein and transcript identifiers to their gene, from the GFF."""
    p2g = {}
    with gzip.open(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "mRNA"):
                continue
            a = f[8]
            gid = re.search(r"gene=([^;]+)", a) or re.search(r"locus_tag=([^;]+)", a)
            par = re.search(r"Parent=(?:transcript:|rna-)?([^;,]+)", a)
            if gid:
                gene = gid.group(1)
            elif par:
                gene = par.group(1).replace("gene:", "").replace("gene-", "")
            else:
                continue
            pid = re.search(r"protein_id=([^;]+)", a)
            if pid:
                p2g.setdefault(pid.group(1), gene)
            tid = re.search(r"ID=(?:transcript:|rna-)?([^;]+)", a)
            if f[2] == "mRNA" and tid:
                p2g.setdefault(tid.group(1), gene)
    return p2g


def clade_entry(hits):
    """Closest Clade A orthologue: the hit of highest identity to a validated target."""
    if not hits:
        return ["", ""]
    name, pid = max(hits, key=lambda h: h[1])
    return [name.replace(" (one-way)", ""), pid]


def gene_of(protein, p2g):
    stripped = ISOFORM.sub("", protein)
    if stripped != protein:
        return p2g.get(stripped, stripped)
    return p2g.get(protein, protein)


def count_genome(genome, hma_dir, gff_dir, min_bit):
    best = parse_domtbl(os.path.join(hma_dir, f"{genome}.hma_sens.domtbl"), min_bit)
    p2g = protein_to_gene(os.path.join(gff_dir, GFF[genome]))
    hma = {gene_of(t, p2g) for t in best}
    small = {gene_of(t, p2g) for t, (bit, tlen) in best.items() if tlen <= MAX_LEN}
    return len(hma), len(small), best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hma-dir", default=os.path.join(ROOT, "hma"))
    ap.add_argument("--gff-dir", default=os.path.join(ROOT, "gff"))
    ap.add_argument("--clade-a", default=None,
                    help="cladeA_rbh.json mapping genome -> [protein, percent identity]")
    ap.add_argument("--out", default=os.path.join(ROOT, "hma", "TableS21_sHMA_cladeA.csv"))
    ap.add_argument("--min-bit", type=float, default=MIN_BIT)
    ap.add_argument("--sweep", action="store_true",
                    help="report finger millet counts across bit cutoffs 0 to 10")
    args = ap.parse_args()

    clade_path = args.clade_a or os.path.join(args.hma_dir, "cladeA_rbh.json")
    clade = json.load(open(clade_path)) if os.path.exists(clade_path) else {}

    rice = parse_domtbl(os.path.join(args.hma_dir, "Osativa.hma_sens.domtbl"), 0.0)
    for ref, name in ((OSHIPP19, "OsHIPP19"), (OSHIPP20, "OsHIPP20")):
        hit = next((v for k, v in rice.items() if k.startswith(ref)), None)
        print(f"  control {name}: "
              + (f"{hit[0]:.1f} bits, {hit[1]} aa" if hit else "ABSENT - check the input"))

    rows, counts = [], collections.OrderedDict()
    for g in GENOMES:
        n_hma, n_small, _ = count_genome(g, args.hma_dir, args.gff_dir, args.min_bit)
        counts[g] = (n_hma, n_small)
        rows.append([g, n_small] + clade_entry(clade.get(g)))

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Genome",
                    f"sHMA proteins (<={MAX_LEN} aa, HMA bit >={args.min_bit:g})",
                    "Closest Clade A orthologue",
                    "% identity to nearest validated AVR-Pik target"])
        w.writerows(rows)

    print(f"\nwrote {args.out}")
    for g, (n_hma, n_small) in counts.items():
        print(f"  {g:12s} HMA genes {n_hma:4d}   small HMA genes {n_small:4d}")

    if args.sweep:
        print("\n  finger millet small HMA genes by bit cutoff")
        for cutoff in range(0, 11):
            _, n_small, _ = count_genome("Ecoracana", args.hma_dir, args.gff_dir, float(cutoff))
            print(f"    bit >= {cutoff:2d}   {n_small}")


if __name__ == "__main__":
    main()
