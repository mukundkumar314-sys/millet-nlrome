import gzip, re, os, sys, json
from collections import defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")

PROT=f"{ROOT}/prot"; GFF=f"{ROOT}/gff"; NLR=f"{ROOT}/nlr"
SPECIES=["Ecoracana","Sitalica","Sviridis","Phallii","Dexilis","Ecrusgalli",
         "Osativa","Sbicolor","Zmays","Bdistachyon"]

def read_fasta(p):
    name=None; buf=[]
    with open(p) as fh:
        for ln in fh:
            if ln.startswith('>'):
                if name: yield name, ''.join(buf)
                name=ln[1:].rstrip('\n'); buf=[]
            else: buf.append(ln.strip())
    if name: yield name, ''.join(buf)

# ---- finger millet: protein -> gene + coords from NCBI GFF ----
def ecoracana_map():
    p2g={}; gcoord={}
    with gzip.open(f"{GFF}/Ecoracana.gff.gz",'rt') as fh:
        for ln in fh:
            if ln.startswith('#'): continue
            f=ln.rstrip('\n').split('\t')
            if len(f)<9: continue
            if f[2]=='CDS':
                a=f[8]
                pid=re.search(r'protein_id=([^;]+)',a)
                par=re.search(r'(?:locus_tag|gene)=([^;]+)',a)
                if pid and par:
                    p2g[pid.group(1)]=par.group(1)
            if f[2]=='gene':
                a=f[8]
                g=re.search(r'(?:locus_tag|ID=gene-)([^;]+)',a)
                if g:
                    gid=g.group(1)
                    gcoord[gid]=(f[0], int(f[3]), int(f[4]), f[6])
    return p2g, gcoord

# ---- parse Ensembl pep header ----
ENS=re.compile(r'(\S+)\s+pep\s+\S+:(?:[^:]+):([^:]+):(\d+):(\d+):(-?1)')
def parse_ens(h):
    m=ENS.match(h)
    g=re.search(r'gene:(\S+)',h)
    if not (m and g): return None
    return dict(prot=m.group(1), chrom=m.group(2), start=int(m.group(3)),
                end=int(m.group(4)), strand='+' if m.group(5)=='1' else '-', gene=g.group(1))

def nbarc_best(sp):
    """best NB-ARC domain per protein: (score, evalue, hmm_from, hmm_to, env_from, env_to)"""
    best={}
    f=f"{NLR}/{sp}.nbarc.domtbl"
    if not os.path.exists(f): return best
    for ln in open(f):
        if ln.startswith('#'): continue
        c=ln.split()
        if len(c)<23: continue
        pid=c[0]; sc=float(c[13]); ev=float(c[12])
        rec=(sc,ev,int(c[15]),int(c[16]),int(c[19]),int(c[20]))
        if pid not in best or sc>best[pid][0]: best[pid]=rec
    return best

summary=[]; all_recs=[]
for sp in SPECIES:
    pep=f"{PROT}/{sp}.pep.fa"
    seqs=dict(read_fasta(pep))
    nb=nbarc_best(sp)
    if sp=="Ecoracana":
        p2g,gcoord=ecoracana_map()
        recs=[]
        for h,s in seqs.items():
            pid=h.split()[0]
            if pid not in nb: continue
            gid=p2g.get(pid)
            if gid is None: continue
            co=gcoord.get(gid)
            if co is None: continue
            recs.append(dict(species=sp,prot=pid,gene=gid,chrom=co[0],start=co[1],
                             end=co[2],strand=co[3],seq=s,
                             nb_score=nb[pid][0],nb_eval=nb[pid][1],
                             nb_hfrom=nb[pid][2],nb_hto=nb[pid][3],
                             nb_efrom=nb[pid][4],nb_eto=nb[pid][5]))
    else:
        recs=[]
        for h,s in seqs.items():
            pid=h.split()[0]
            if pid not in nb: continue
            d=parse_ens(h)
            if d is None: continue
            recs.append(dict(species=sp,prot=pid,gene=d['gene'],chrom=d['chrom'],
                             start=d['start'],end=d['end'],strand=d['strand'],seq=s,
                             nb_score=nb[pid][0],nb_eval=nb[pid][1],
                             nb_hfrom=nb[pid][2],nb_hto=nb[pid][3],
                             nb_efrom=nb[pid][4],nb_eto=nb[pid][5]))
    # collapse: longest protein per gene
    bygene={}
    for r in recs:
        g=r['gene']
        if g not in bygene or len(r['seq'])>len(bygene[g]['seq']): bygene[g]=r
    genes=list(bygene.values())
    all_recs.extend(genes)
    summary.append((sp,len(recs),len(genes)))

print(f"{'species':<14}{'nbarc_prot':>11}{'nbarc_gene':>11}")
for sp,a,b in summary: print(f"{sp:<14}{a:>11}{b:>11}")
print(f"{'TOTAL':<14}{sum(a for _,a,_ in summary):>11}{sum(b for _,_,b in summary):>11}")

os.makedirs(NLR,exist_ok=True)
with open(f"{NLR}/nbarc_genes.json","w") as fh: json.dump(all_recs,fh)
with open(f"{NLR}/nbarc_candidates.fa","w") as fh:
    for r in all_recs: fh.write(f">{r['species']}|{r['gene']}|{r['prot']}\n{r['seq']}\n")
print("\nwrote nbarc_candidates.fa and nbarc_genes.json")
