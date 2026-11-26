import json, os, re, glob
from collections import defaultdict, Counter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")

NLR=f"{ROOT}/nlr"
recs=json.load(open(f"{NLR}/nbarc_genes.json"))
byid={f"{r['species']}|{r['gene']}|{r['prot']}":r for r in recs}

# ---------------- domain category curation ----------------
CORE_NB   = {"NB-ARC","WHD_DRP","WHD_plant_disease"}
CORE_CC   = {"Rx_N"}
CORE_TIR  = {"TIR","TIR_2"}
CORE_RPW8 = {"RPW8"}
# NB-ARC-adjacent P-loop NTPase / STAND folds -> structural, never an ID
STAND = {"AAA_22","AAA_16","AAA","AAA_14","AAA_2","AAA_5","AAA_lid_9","AAA_19","AAA_28",
         "NACHT","nSTAND_NTPase5","Sigma54_activat","ATPase_2","ClpB_D2-small","MLKL_N",
         "FNIP","NPHP3_N","Pox_A32","ABC_tran","Mg_chelatase"}
TE_FAM = {"RVT_1","RVT_2","RVT_3","RT_RNaseH","RT_RNaseH_2","rve","Integrase_H2C2",
          "gag_pre-integrs","Retrotran_gag_2","Retrotran_gag_3","Ty3_capsid","SH3_retrovirus",
          "SH3_Tf2-1","RVP_2","zf-RVT","zf-CCHC","Transposase_21","Transposase_24",
          "Transpos_assoc","DDE_Tnp_4","TniB","Transp_inhibit","POL3_N","Pol_BBD",
          "DUF4218","DUF4216","Dimer_Tnp_hAT","Peptidase_A17"}
def is_lrr(n): return n.startswith("LRR") or n in {"LRRNT_2"}

# ---------------- load Pfam hits ----------------
hits=defaultdict(list)
for f in glob.glob(f"{NLR}/pfam_out/*.domtbl"):
    for ln in open(f):
        if ln.startswith('#'): continue
        c=ln.split()
        if len(c)<23: continue
        tname,qname,qacc=c[0],c[3],c[4]
        hits[tname].append(dict(fam=qname,acc=qacc,ev=float(c[12]),score=float(c[13]),
                                afrom=int(c[17]),ato=int(c[18])))

def category(fam):
    if fam in CORE_NB:   return "NB"
    if fam in CORE_CC:   return "CC"
    if fam in CORE_TIR:  return "TIR"
    if fam in CORE_RPW8: return "RPW8"
    if is_lrr(fam):      return "LRR"
    if fam in STAND:     return "STAND"
    if fam in TE_FAM:    return "TE"
    return "ID"

rows=[]; idcount=Counter(); idspec=defaultdict(set); te_flag=0
for gid,r in byid.items():
    hh=sorted(hits.get(gid,[]),key=lambda x:x['afrom'])
    cats={category(h['fam']) for h in hh}
    ids=sorted({h['fam'] for h in hh if category(h['fam'])=="ID"})
    tes=sorted({h['fam'] for h in hh if category(h['fam'])=="TE"})
    # architecture prefix
    if   "TIR"  in cats: pre="T"
    elif "RPW8" in cats: pre="R"
    elif "CC"   in cats: pre="C"
    else: pre=""
    arch = pre + "N" + ("L" if "LRR" in cats else "") + ("X" if ids else "")
    dom_order = "-".join(h['fam'] for h in hh if category(h['fam']) in {"NB","CC","LRR","TIR","RPW8","ID"})
    for i in ids: idcount[i]+=1; idspec[i].add(r['species'])
    if tes: te_flag+=1
    rows.append(dict(id=gid, species=r['species'], gene=r['gene'], prot=r['prot'],
                     chrom=r['chrom'], start=r['start'], end=r['end'], strand=r['strand'],
                     length=len(r['seq']), nb_score=r['nb_score'], nb_eval=r['nb_eval'],
                     arch=arch, has_CC=("CC" in cats), has_LRR=("LRR" in cats),
                     n_id=len(ids), ids=";".join(ids), te_assoc=";".join(tes),
                     domain_order=dom_order))

json.dump(rows,open(f"{NLR}/nlr_classified.json","w"))

SP=["Ecoracana","Sitalica","Sviridis","Phallii","Dexilis","Ecrusgalli","Osativa","Sbicolor","Zmays","Bdistachyon"]
print("=== NLR architecture counts per species ===")
archs=sorted({r['arch'] for r in rows})
print(f"{'species':<13}"+"".join(f"{a:>7}" for a in archs)+f"{'TOTAL':>8}")
for sp in SP:
    sub=[r for r in rows if r['species']==sp]; c=Counter(r['arch'] for r in sub)
    print(f"{sp:<13}"+"".join(f"{c.get(a,0):>7}" for a in archs)+f"{len(sub):>8}")
c=Counter(r['arch'] for r in rows)
print(f"{'ALL':<13}"+"".join(f"{c.get(a,0):>7}" for a in archs)+f"{len(rows):>8}")

print(f"\n=== Integrated domains: {len(idcount)} distinct families ===")
print(f"{'domain':<28}{'genes':>7}{'species':>9}")
for d,n in idcount.most_common():
    print(f"{d:<28}{n:>7}{len(idspec[d]):>9}")
print(f"\nTE-associated NB-ARC genes (flagged): {te_flag}")
