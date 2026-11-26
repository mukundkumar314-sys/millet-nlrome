import json, glob, gzip, re
from collections import defaultdict, Counter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
NLR=f"{ROOT}/nlr"
exec(open(f"{ROOT}/scripts/02_classify.py").read().split("# ---------------- load Pfam hits")[0].split("recs=json.load")[0].replace('import json, os, re, glob','import json, os, re, glob'))
CORE_NB={"NB-ARC","WHD_DRP","WHD_plant_disease"}; CORE_CC={"Rx_N"}
CORE_TIR={"TIR","TIR_2"}; CORE_RPW8={"RPW8"}
STAND={"AAA_22","AAA_16","AAA","AAA_14","AAA_2","AAA_5","AAA_lid_9","AAA_19","AAA_28","NACHT",
       "nSTAND_NTPase5","Sigma54_activat","ATPase_2","ClpB_D2-small","MLKL_N","FNIP","NPHP3_N",
       "Pox_A32","ABC_tran","Mg_chelatase"}
TE_FAM={"RVT_1","RVT_2","RVT_3","RT_RNaseH","RT_RNaseH_2","rve","Integrase_H2C2","gag_pre-integrs",
        "Retrotran_gag_2","Retrotran_gag_3","Ty3_capsid","SH3_retrovirus","SH3_Tf2-1","RVP_2",
        "zf-RVT","zf-CCHC","Transposase_21","Transposase_24","Transpos_assoc","DDE_Tnp_4","TniB",
        "Transp_inhibit","POL3_N","Pol_BBD","DUF4218","DUF4216","Dimer_Tnp_hAT","Peptidase_A17"}
def is_lrr(n): return n.startswith("LRR") or n in {"LRRNT_2"}
def category(f):
    if f in CORE_NB: return "NB"
    if f in CORE_CC: return "CC"
    if f in CORE_TIR: return "TIR"
    if f in CORE_RPW8: return "RPW8"
    if is_lrr(f): return "LRR"
    if f in STAND: return "STAND"
    if f in TE_FAM: return "TE"
    return "ID"

hits=defaultdict(list)
for f in glob.glob(f"{NLR}/pfam_proso/*.domtbl"):
    for ln in open(f):
        if ln.startswith('#'): continue
        c=ln.split()
        if len(c)<23: continue
        hits[c[0]].append(dict(fam=c[3],afrom=int(c[17]),ato=int(c[18]),ev=float(c[13])))
pg=json.load(open(f"{NLR}/proso_genes.json"))
new=[]
for r in pg:
    gid=f"{r['species']}|{r['gene']}|{r['prot']}"
    hh=sorted(hits.get(gid,[]),key=lambda x:x['afrom'])
    cats={category(h['fam']) for h in hh}
    ids=sorted({h['fam'] for h in hh if category(h['fam'])=="ID"})
    tes=sorted({h['fam'] for h in hh if category(h['fam'])=="TE"})
    pre="T" if "TIR" in cats else ("R" if "RPW8" in cats else ("C" if "CC" in cats else ""))
    arch=pre+"N"+("L" if "LRR" in cats else "")+("X" if ids else "")
    new.append(dict(id=gid,species=r["species"],gene=r["gene"],prot=r["prot"],chrom=r["chrom"],
                    start=r["start"],end=r["end"],strand=r["strand"],length=len(r["seq"]),
                    nb_score=r["nb_score"],nb_eval=r["nb_eval"],arch=arch,
                    has_CC=("CC" in cats),has_LRR=("LRR" in cats),n_id=len(ids),
                    ids=";".join(ids),te_assoc=";".join(tes),
                    domain_order="~".join(h['fam'] for h in hh if category(h['fam']) in {"NB","CC","LRR","TIR","RPW8","ID"}),
                    domains=[{"fam":h["fam"],"s":h["afrom"],"e":h["ato"],"ev":h["ev"]} for h in hh]))
# chromosome labels
acc2chr={}
for ln in open(f"{ROOT}/gff/Pmiliaceum.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': acc2chr[f[4]]=f[2]
for r in new: r["chrom_lab"]=acc2chr.get(r["chrom"],r["chrom"])

rows=json.load(open(f"{NLR}/nlr_classified.json"))
rows=[r for r in rows if r["species"]!="Pmiliaceum"]+new
json.dump(rows,open(f"{NLR}/nlr_classified.json","w"))
print("proso NLRs added:",len(new),"  total NLRs:",len(rows))
print("proso architectures:",dict(Counter(r["arch"] for r in new)))
print("proso ID families:",Counter(i for r in new for i in r["ids"].split(";") if i).most_common(8))

# gene coords + genome size
allg=json.load(open(f"{NLR}/all_genes.json"))
g=[]
with gzip.open(f"{ROOT}/gff/Pmiliaceum.gff.gz",'rt') as fh:
    for ln in fh:
        if ln.startswith('#'): continue
        f=ln.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='gene': continue
        m=re.search(r'locus_tag=([^;]+)',f[8]) or re.search(r'ID=gene-([^;]+)',f[8])
        if m: g.append((m.group(1),acc2chr.get(f[0],f[0]),int(f[3]),int(f[4])))
allg["Pmiliaceum"]=g
json.dump(allg,open(f"{NLR}/all_genes.json","w"))
sz=json.load(open(f"{NLR}/genome_sizes.json"))
tot=0
for ln in open(f"{ROOT}/gff/Pmiliaceum.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)>=10 and f[8].isdigit(): tot+=int(f[8])
sz["Pmiliaceum"]=tot
json.dump(sz,open(f"{NLR}/genome_sizes.json","w"))
print(f"proso genes={len(g)}  assembly={tot/1e6:.1f} Mb  NLR density={len(new)/(tot/1e6)*100:.1f} per 100 Mb")
with open(f"{NLR}/nbarc_candidates.fa","a") as fh:
    for r in pg: fh.write(f">{r['species']}|{r['gene']}|{r['prot']}\n{r['seq']}\n")
print("appended proso to nbarc_candidates.fa; total:",sum(1 for l in open(f"{NLR}/nbarc_candidates.fa") if l.startswith('>')))
