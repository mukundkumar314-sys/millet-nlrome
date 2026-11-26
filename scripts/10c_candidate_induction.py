#!/usr/bin/env python3
"""Stage 10c - response of pre-specified candidate genes to blast inoculation.

Rebuilds Supplementary Table S23 from the salmon quantifications of stage 10b.

BioProject PRJNA648385 provides one library per genotype, timepoint and
treatment, so genotype serves as the replicate unit for the treatment contrast,
giving twelve paired inoculated-against-mock comparisons across four genotypes
and three timepoints. The 0 h libraries are mock-only and form a baseline.

Testing is confined to the 24 genes named in target_sets.json - the cloned blast
resistance gene orthologues, the members of finger millet head-to-head pairs and
the two Clade A HMA proteins - so the analysis supports no genome-wide
differential expression claim.
"""
import os, csv, json, math, argparse, collections, statistics

ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
PSEUDO = 1.0          # added to TPM before taking the log ratio


def read_quant(path):
    """Transcript -> (TPM, NumReads) for one salmon run."""
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[0]] = (float(f[3]), float(f[4]))
    return out


def gene_of(transcript, t2g):
    """Transcript identifier to locus tag, via the map built from the CDS FASTA headers."""
    return t2g.get(transcript)


def wilcoxon_signed_rank(values):
    """Exact two-sided signed-rank p value; ties dropped, as scipy does by default."""
    vals = [v for v in values if v != 0]
    n = len(vals)
    if n == 0:
        return 1.0
    order = sorted(range(n), key=lambda i: abs(vals[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(vals[order[j + 1]]) == abs(vals[order[i]]):
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    w_plus = sum(r for r, v in zip(ranks, vals) if v > 0)
    w_minus = sum(r for r, v in zip(ranks, vals) if v < 0)
    w = min(w_plus, w_minus)
    # exact distribution by enumeration over sign assignments
    counts = {0: 1}
    for r in range(1, n + 1):
        nxt = collections.Counter()
        for s, c in counts.items():
            nxt[s] += c
            nxt[s + r] += c
        counts = nxt
    total = 2 ** n
    tail = sum(c for s, c in counts.items() if s <= w)
    return min(1.0, 2.0 * tail / total)


def benjamini_hochberg(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = pvals[i] * n / (n - rank + 1)
        prev = min(prev, val)
        q[i] = min(1.0, prev)
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rna-dir", default=os.path.join(ROOT, "rna"))
    ap.add_argument("--out", default=os.path.join(ROOT, "rna", "TableS23_candidate_induction.csv"))
    args = ap.parse_args()

    design = json.load(open(os.path.join(args.rna_dir, "design.json")))
    targets = json.load(open(os.path.join(args.rna_dir, "target_sets.json")))
    t2g_path = os.path.join(args.rna_dir, "t2g_all.json")
    t2g = json.load(open(t2g_path)) if os.path.exists(t2g_path) else {}

    genes = [(cat, g) for cat, gs in targets.items() for g in gs]
    wanted = {g for _, g in genes}
    print(f"  {len(design)} libraries, {len(wanted)} pre-specified genes")

    # gene-level TPM and read counts per run
    tpm = collections.defaultdict(dict)
    reads = collections.defaultdict(dict)
    for d in design:
        qf = os.path.join(args.rna_dir, "quant", d["run"], "quant.sf")
        if not os.path.exists(qf):
            print(f"  missing {qf}")
            continue
        for tx, (t, r) in read_quant(qf).items():
            g = gene_of(tx, t2g)
            if g in wanted:
                tpm[g][d["run"]] = tpm[g].get(d["run"], 0.0) + t
                reads[g][d["run"]] = reads[g].get(d["run"], 0.0) + r

    # pair each inoculated library with the mock library of the same genotype and timepoint
    by_key = {(d["geno"], d["time"], d["treat"]): d["run"] for d in design}
    pairs = [(by_key[(g, t, "inoc")], by_key[(g, t, "ctrl")], g, t)
             for (g, t, tr) in list(by_key) if tr == "inoc"
             and (g, t, "ctrl") in by_key]
    pairs.sort(key=lambda p: (p[2], p[3]))
    print(f"  {len(pairs)} paired inoculated-against-mock comparisons")

    rows, pvals = [], []
    for cat, g in genes:
        lfc = []
        for inoc, ctrl, _, _ in pairs:
            a = tpm[g].get(inoc, 0.0)
            b = tpm[g].get(ctrl, 0.0)
            lfc.append(math.log2((a + PSEUDO) / (b + PSEUDO)))
        n_up = sum(1 for v in lfc if v > 0)
        p = wilcoxon_signed_rank(lfc)
        mean_ctrl = statistics.mean([tpm[g].get(c, 0.0) for _, c, _, _ in pairs])
        mean_inoc = statistics.mean([tpm[g].get(i, 0.0) for i, _, _, _ in pairs])
        paired_runs = {r for p in pairs for r in p[:2]}
        total_reads = math.floor(sum(v for r, v in reads[g].items() if r in paired_runs))
        rows.append([cat, g, f"{statistics.median(lfc):.2f}", f"{n_up}/{len(pairs)}",
                     f"{mean_ctrl:.2f}", f"{mean_inoc:.2f}", total_reads,
                     f"{p:.4f}"])
        pvals.append(p)

    for row, q in zip(rows, benjamini_hochberg(pvals)):
        row.append(f"{q:.3f}")

    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Category", "Gene", "Median log2FC (inoc vs mock)", "Pairs up / total",
                    "Mean TPM mock", "Mean TPM inoculated", "Total reads across the 24 paired libraries",
                    "Wilcoxon p", "BH q"])
        w.writerows(rows)
    print(f"\nwrote {args.out}")
    for r in rows:
        if float(r[8]) < 0.05:
            print(f"  passes BH 0.05: {r[1]}  median log2FC {r[2]}  q {r[8]}")


if __name__ == "__main__":
    main()
