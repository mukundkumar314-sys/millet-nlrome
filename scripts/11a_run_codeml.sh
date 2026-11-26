#!/usr/bin/env bash
set -euo pipefail
cd "${PROJECT_ROOT:-/home/claude}"/paml
for og in OG0021 OG0025 OG0083 OG0213; do
  for mod in M1a M2a M7 M8; do
    case $mod in
      M1a) NC=1;; M2a) NC=2;; M7) NC=7;; M8) NC=8;;
    esac
    d=run/$og/$mod; mkdir -p $d
    cat > $d/codeml.ctl <<CTL
      seqfile = ../../../aln/$og.codon.phy
     treefile = ../../../aln/$og.treefile
      outfile = out.txt
        noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = $NC
        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.5
        getSE = 0
 RateAncestor = 0
    cleandata = 0
CTL
    (cd $d && codeml codeml.ctl > log.txt 2>&1)
    echo "done $og $mod"
  done
done
echo CODEML_DONE
