#!/usr/bin/env python3
"""Stage 09e - effector family occurrence by host lineage, and the annotation gap.

Rebuilds Supplementary Table S22 from the hmmsearch table produced by stage 09c

    hmmsearch --tblout fam.tbl --cut_ga <families>.hmm all_prot.fa

where all_prot.fa holds every protein of the downloaded Pyricularia oryzae
proteomes, each header prefixed with its assembly accession and a pipe.

The canonical AVR-Pik family gene is additionally scored from the assemblies
with tblastn (stage 09d), because gene prediction frequently overlooks small
effector genes: in the Eleusine-infecting lineage the gene is annotated in one
isolate yet present in the assembly of every one.
"""
import os, csv, json, argparse, statistics, collections

ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
LINEAGES = ["Oryza", "Eleusine", "Triticum"]
FAMILIES = ["AVR-Pik_HID", "MAX_effector", "AvrPiz-t", "AvrPi9", "Avr_Pia"]
HMM_NAME = {"AVRPikFam": "AVR-Pik_HID", "MAX": "MAX_effector", "AvrPiz-t": "AvrPiz-t",
            "AvrPi9": "AvrPi9", "Avr-Pia": "Avr_Pia", "Avr_Pia": "Avr_Pia"}
MIN_PIDENT = 90.0      # a genomic hit counts as the canonical gene above this identity
MIN_LEN = 100          # over at least this many aligned residues


def load_lineages(groups_path, acc_path):
    """Accession -> host lineage, with the reference isolate 70-15 placed in Oryza."""
    groups = json.load(open(groups_path))
    lin = {}
    for name, accs in groups.items():
        for a in accs:
            lin[a] = name
    accs = [l.strip() for l in open(acc_path) if l.strip()]
    for a in accs:
        if lin.get(a, "other/unknown") == "other/unknown":
            lin[a] = "Oryza"           # 70-15: its BioSample carries no host field
    return {a: lin[a] for a in accs}


def parse_tblout(path):
    """Assembly accession -> family -> copy number."""
    counts = collections.defaultdict(collections.Counter)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 4:
                continue
            target, query = f[0], f[2]
            fam = HMM_NAME.get(query, query)
            acc = target.split("|", 1)[0]
            counts[acc][fam] += 1
    return counts


def parse_tblastn(path, lineage):
    """Assemblies in which the canonical AVR-Pik family gene is detectable."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            acc, pident, length = f[0], float(f[3]), int(f[4])
            if pident >= MIN_PIDENT and length >= MIN_LEN:
                best[acc] = max(best.get(acc, 0.0), pident)
    return {a: p for a, p in best.items() if a in lineage}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mag-dir", default=os.path.join(ROOT, "mag"))
    ap.add_argument("--counts", default=None,
                    help="max_counts.json from stage 09c; parsed from the hmmsearch "
                         "table instead when absent")
    ap.add_argument("--out", default=os.path.join(ROOT, "mag", "TableS22_effector_families.csv"))
    args = ap.parse_args()

    lineage = load_lineages(os.path.join(args.mag_dir, "groups.json"),
                            os.path.join(args.mag_dir, "acc.txt"))
    counts_path = args.counts or os.path.join(args.mag_dir, "max_counts.json")
    if os.path.exists(counts_path):
        per_family = json.load(open(counts_path))
        counts = collections.defaultdict(collections.Counter)
        for fam, d in per_family.items():
            for a, n in d.items():
                counts[a][fam] = n
        print(f"  counts read from {counts_path}")
    else:
        counts = parse_tblout(os.path.join(args.mag_dir, "fam.tbl"))
        print("  counts parsed from fam.tbl")

    examined = collections.Counter(lineage.values())
    print("  isolates examined: " + ", ".join(f"{k} {examined[k]}" for k in LINEAGES)
          + f", total {sum(examined.values())}")

    rows = []
    for fam in FAMILIES:
        for lin in LINEAGES:
            accs = [a for a, l in lineage.items() if l == lin]
            per = [counts[a][fam] for a in accs]
            carriers = [c for c in per if c > 0]
            rows.append([fam, lin, len(accs), len(carriers), sum(per),
                         int(statistics.median(carriers)) if carriers else 0])

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Pfam family", "Host lineage", "Isolates examined",
                    "Isolates with >=1 member", "Total copies", "Median copies per carrier"])
        w.writerows(rows)
    print(f"\nwrote {args.out}")

    tb = os.path.join(args.mag_dir, "tblastn_all.tsv")
    if os.path.exists(tb):
        hits = parse_tblastn(tb, lineage)
        el = {a: p for a, p in hits.items() if lineage[a] == "Eleusine"}
        n_el = sum(1 for a, l in lineage.items() if l == "Eleusine")
        annotated = sum(1 for a, l in lineage.items()
                        if l == "Eleusine" and counts[a]["AVR-Pik_HID"] > 0)
        print(f"\n  annotation gap, Eleusine-infecting lineage")
        print(f"    canonical gene in the assembly : {len(el)} of {n_el}")
        print(f"    family member in the proteome  : {annotated} of {n_el}")
        if el:
            print(f"    genomic identity range         : {min(el.values()):.1f} to {max(el.values()):.1f}%")


if __name__ == "__main__":
    main()
