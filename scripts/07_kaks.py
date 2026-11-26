import json, random, warnings, itertools
warnings.filterwarnings("ignore")
random.seed(3)
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.codonalign.codonseq import cal_dn_ds, CodonSeq
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")

cds=json.load(open(f"{ROOT}/nlr/nlr_cds_byid.json"))
rows=json.load(open(f"{ROOT}/nlr/nlr_classified.json"))
pan=json.load(open(f"{ROOT}/nlr/pan_nlrome.json"))
byid={r["id"]:r for r in rows}
key={f'{r["species"]}|{r["gene"]}':r for r in rows}

al=PairwiseAligner(); al.substitution_matrix=substitution_matrices.load("BLOSUM62")
al.open_gap_score=-11; al.extend_gap_score=-1; al.mode="global"

def prot(s):
    return str(Seq(s[:len(s)//3*3]).translate()).rstrip("*")

def kaks(c1,c2):
    p1,p2=prot(c1),prot(c2)
    if len(p1)<80 or len(p2)<80: return None
    try: aln=al.align(p1,p2)[0]
    except Exception: return None
    A,B=str(aln[0]),str(aln[1])
    oa=ob=0; ca=[]; cb=[]
    for x,y in zip(A,B):
        if x!="-" and y!="-":
            s1=c1[oa*3:oa*3+3]; s2=c2[ob*3:ob*3+3]
            if len(s1)==3 and len(s2)==3 and "N" not in (s1+s2).upper():
                ca.append(s1); cb.append(s2)
        if x!="-": oa+=1
        if y!="-": ob+=1
    if len(ca)<80: return None
    try: dn,ds=cal_dn_ds(CodonSeq("".join(ca)),CodonSeq("".join(cb)),method="YN00")
    except Exception: return None
    if ds is None or dn is None or ds<=0.01 or ds>2 or dn>2: return None
    return dn,ds,dn/ds

MAXP=12
pergene={}; perog=[]
for og in pan:
    mem=[m for m in og["members"] if m in cds]
    if len(mem)<2: continue
    allp=list(itertools.combinations(mem,2))
    if len(allp)>MAXP: allp=random.sample(allp,MAXP)
    vals=[]
    for a,b in allp:
        r=kaks(cds[a],cds[b])
        if r is None: continue
        vals.append(r[2])
        pergene.setdefault(a,[]).append(r[2]); pergene.setdefault(b,[]).append(r[2])
    if vals:
        perog.append(dict(og=og["og"],cat=og["cat"],n=len(mem),
                          w=sorted(vals)[len(vals)//2],npairs=len(vals)))
json.dump(perog,open(f"{ROOT}/nlr/kaks_og.json","w"))
json.dump({k:sum(v)/len(v) for k,v in pergene.items()},open(f"{ROOT}/nlr/kaks_gene.json","w"))
print("orthogroups with Ka/Ks:",len(perog))
print("genes with Ka/Ks     :",len(pergene))
import statistics as st
from collections import Counter
print("by category:",Counter(o["cat"] for o in perog))
