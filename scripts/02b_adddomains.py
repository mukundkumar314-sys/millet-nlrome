import json, glob
from collections import defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
NLR=f"{ROOT}/nlr"
rows=json.load(open(f"{NLR}/nlr_classified.json"))
hits=defaultdict(list)
for f in glob.glob(f"{NLR}/pfam_out/*.domtbl"):
    for ln in open(f):
        if ln.startswith('#'): continue
        c=ln.split()
        if len(c)<23: continue
        hits[c[0]].append((c[3],int(c[17]),int(c[18]),float(c[13])))
n=0
for r in rows:
    hh=sorted(hits.get(r["id"],[]),key=lambda x:x[1])
    r["domains"]=[{"fam":f,"s":a,"e":b,"ev":e} for f,a,b,e in hh]
    if hh: n+=1
json.dump(rows,open(f"{NLR}/nlr_classified.json","w"))
print("rows with domain coords:",n,"/",len(rows))
ex=[r for r in rows if r["n_id"]>0][0]
print("example:",ex["id"],[(d["fam"],d["s"],d["e"]) for d in ex["domains"]][:6])
