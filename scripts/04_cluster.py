import json, random, statistics as st
from collections import defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
random.seed(42)
NLR=f"{ROOT}/nlr"
rows=json.load(open(f"{NLR}/nlr_classified.json"))
allg=json.load(open(f"{NLR}/all_genes.json"))

# finger millet accession -> chromosome label (A/B subgenome)
acc2chr={}
for ln in open(f"{ROOT}/gff/Ecoracana.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': acc2chr[f[4]]=f[2]
import json as _j
_pm=_j.load(open(f"{ROOT}/nlr/proso_acc2chr.json"))
for r in rows:
    if r['species']=='Ecoracana': r['chrom_lab']=acc2chr.get(r['chrom'], r['chrom'])
    elif r['species']=='Pmiliaceum': r['chrom_lab']=_pm.get(r['chrom'], r['chrom'])
    else: r['chrom_lab']=r['chrom']
for sp,gl in allg.items():
    if sp=='Ecoracana':
        allg[sp]=[(g,acc2chr.get(c,c),s,e) for g,c,s,e in gl]

SP=["Ecoracana","Pmiliaceum","Sitalica","Sviridis","Phallii","Dexilis","Ecrusgalli","Osativa","Sbicolor","Zmays","Bdistachyon"]
GRID=[10,20,30,50,75,100,150,200,300,500]  # kb

def consec(pos):
    out=[]
    for c,ps in pos.items():
        ps=sorted(ps)
        out += [ps[i+1]-ps[i] for i in range(len(ps)-1)]
    return out

def frac_le(d,g): return sum(1 for x in d if x<=g*1000)/len(d) if d else 0.0

res={}
print(f"{'species':<13}{'NLR':>5}{'g(kb)':>7}{'enr':>7}{'clusters':>9}{'in_clust':>9}{'%clust':>8}{'max':>5}")
for sp in SP:
    sub=[r for r in rows if r['species']==sp]
    pos=defaultdict(list)
    for r in sub: pos[r['chrom_lab']].append((r['start']+r['end'])//2)
    obs=consec(pos)
    # null: random genes, same per-chrom counts
    bychr=defaultdict(list)
    for g,c,s,e in allg[sp]: bychr[c].append((s+e)//2)
    exp_fracs={g:[] for g in GRID}
    for _ in range(100):
        npos=defaultdict(list)
        for c,lst in pos.items():
            pool=bychr.get(c,[])
            if len(pool)>=len(lst): npos[c]=random.sample(pool,len(lst))
        nd=consec(npos)
        for g in GRID: exp_fracs[g].append(frac_le(nd,g))
    # pick largest g with enrichment >= 2
    chosen=GRID[0]; chosen_enr=0
    for g in GRID:
        e=st.mean(exp_fracs[g]) or 1e-9
        enr=frac_le(obs,g)/e
        if enr>=2.0: chosen, chosen_enr = g, enr
    # single-linkage clusters at chosen g
    nclu=0; inclu=0; mx=0; clusters=[]
    for c,lst in pos.items():
        ps=sorted(lst); cur=[ps[0]]
        for i in range(1,len(ps)):
            if ps[i]-ps[i-1] <= chosen*1000: cur.append(ps[i])
            else:
                if len(cur)>=2: clusters.append((c,len(cur))); nclu+=1; inclu+=len(cur); mx=max(mx,len(cur))
                cur=[ps[i]]
        if len(cur)>=2: clusters.append((c,len(cur))); nclu+=1; inclu+=len(cur); mx=max(mx,len(cur))
    res[sp]=dict(n=len(sub), g=chosen, enr=round(chosen_enr,2), nclusters=nclu,
                 in_cluster=inclu, pct=round(100*inclu/len(sub),1), maxsize=mx,
                 clusters=clusters)
    print(f"{sp:<13}{len(sub):>5}{chosen:>7}{chosen_enr:>7.2f}{nclu:>9}{inclu:>9}{100*inclu/len(sub):>7.1f}%{mx:>5}")

json.dump({k:{kk:vv for kk,vv in v.items()} for k,v in res.items()},
          open(f"{NLR}/clusters.json","w"))
json.dump(rows,open(f"{NLR}/nlr_classified.json","w"))

# ---- finger millet subgenome split ----
print("\n=== Finger millet A vs B subgenome ===")
fm=[r for r in rows if r['species']=='Ecoracana']
import re
A=[r for r in fm if re.fullmatch(r'\dA',str(r['chrom_lab']))]
B=[r for r in fm if re.fullmatch(r'\dB',str(r['chrom_lab']))]
U=[r for r in fm if r not in A and r not in B]
print(f"A subgenome: {len(A)}   B subgenome: {len(B)}   unplaced/scaffold: {len(U)}")
from collections import Counter
print("A arch:",dict(Counter(r['arch'] for r in A)))
print("B arch:",dict(Counter(r['arch'] for r in B)))
print("\nper-chromosome NLR counts (finger millet):")
cc=Counter(r['chrom_lab'] for r in fm)
for k in sorted([x for x in cc if re.fullmatch(r'\d[AB]',str(x))]): print(f"  {k}: {cc[k]}")
