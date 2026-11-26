import gzip,re,json
from collections import defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
PROT=f"{ROOT}/prot"; GFF=f"{ROOT}/gff"; NLR=f"{ROOT}/nlr"
def read_fasta(p):
    name=None;buf=[]
    for ln in open(p):
        if ln.startswith('>'):
            if name: yield name,''.join(buf)
            name=ln[1:].rstrip('\n');buf=[]
        else: buf.append(ln.strip())
    if name: yield name,''.join(buf)
p2g={};gcoord={}
with gzip.open(f"{GFF}/Pmiliaceum.gff.gz",'rt') as fh:
    for ln in fh:
        if ln.startswith('#'): continue
        f=ln.rstrip('\n').split('\t')
        if len(f)<9: continue
        if f[2]=='CDS':
            pid=re.search(r'protein_id=([^;]+)',f[8]); par=re.search(r'(?:locus_tag|gene)=([^;]+)',f[8])
            if pid and par: p2g[pid.group(1)]=par.group(1)
        if f[2]=='gene':
            g=re.search(r'locus_tag=([^;]+)',f[8]) or re.search(r'ID=gene-([^;]+)',f[8])
            if g: gcoord[g.group(1)]=(f[0],int(f[3]),int(f[4]),f[6])
best={}
for ln in open(f"{NLR}/Pmiliaceum.nbarc.domtbl"):
    if ln.startswith('#'): continue
    c=ln.split()
    if len(c)<23: continue
    sc=float(c[13])
    if c[0] not in best or sc>best[c[0]][0]:
        best[c[0]]=(sc,float(c[12]),int(c[15]),int(c[16]),int(c[19]),int(c[20]))
recs=[]
for h,s in read_fasta(f"{PROT}/Pmiliaceum.pep.fa"):
    pid=h.split()[0]
    if pid not in best: continue
    g=p2g.get(pid); co=gcoord.get(g) if g else None
    if co is None: continue
    recs.append(dict(species="Pmiliaceum",prot=pid,gene=g,chrom=co[0],start=co[1],end=co[2],
                     strand=co[3],seq=s,nb_score=best[pid][0],nb_eval=best[pid][1],
                     nb_hfrom=best[pid][2],nb_hto=best[pid][3],nb_efrom=best[pid][4],nb_eto=best[pid][5]))
bg={}
for r in recs:
    if r['gene'] not in bg or len(r['seq'])>len(bg[r['gene']]['seq']): bg[r['gene']]=r
genes=list(bg.values())
print(f"proso: {len(recs)} proteins -> {len(genes)} genes")
json.dump(genes,open(f"{NLR}/proso_genes.json","w"))
with open(f"{NLR}/proso_candidates.fa","w") as fh:
    for r in genes: fh.write(f">{r['species']}|{r['gene']}|{r['prot']}\n{r['seq']}\n")
