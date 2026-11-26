import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
from collections import Counter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; pan=d["pan"]
pik=json.load(open(f"{ROOT}/refnlr/pik_data.json"))
hma_sa=json.load(open(f"{ROOT}/refnlr/hma_standalone.json"))
hma_loci=json.load(open(f"{ROOT}/genome/hma_loci.json"))
ORD=ORD11; HMA_C="#D55E00"
DC={"CC":"#999933","NB":"#0072B2","WHD":"#56B4E9","LRR":"#009E73","HMA":HMA_C}
def dcat(f):
    if f=="Rx_N": return "CC"
    if f=="NB-ARC": return "NB"
    if f in ("WHD_DRP","WHD_plant_disease"): return "WHD"
    if f.startswith("LRR"): return "LRR"
    return None

fig=plt.figure(figsize=(W2,W2*0.60))
gs=fig.add_gridspec(2,12,height_ratios=[0.80,1.0],hspace=0.60,wspace=4.0,
                    left=0.080,right=0.980,top=0.905,bottom=0.085)

# a: rice Pik locus
ax=fig.add_subplot(gs[0,0:4]); ax.axis("off")
p1=pik["Os11t0688832-01"]; p2=pik["Os11t0689100-01"]
lo=min(p1["start"],p2["start"]); hi=max(p1["end"],p2["end"]); span=hi-lo
X=lambda v:(v-lo)/span
ax.plot([0,1],[0.44,0.44],color="#c8c8c8",lw=1.0,zorder=1)
for p,nm,c in [(p1,"$Pik$-1  sensor",HMA_C),(p2,"$Pik$-2  helper","#0072B2")]:
    a,b=X(p["start"]),X(p["end"])
    ax.add_patch(plt.Rectangle((a,0.32),b-a,0.24,facecolor=c,alpha=0.85,lw=0,zorder=3))
    ar=-1 if p["strand"]=="-" else 1; xs=b if ar==1 else a
    ax.annotate("",xy=(xs+0.05*ar,0.44),xytext=(xs,0.44),zorder=4,
                arrowprops=dict(arrowstyle="-|>",color=c,lw=1.0,mutation_scale=6))
    ax.text((a+b)/2,0.62,nm,ha="center",fontsize=5.2,color=INK)
ax.annotate("",xy=(X(p2["start"]),0.44),xytext=(X(p1["end"]),0.44),
            arrowprops=dict(arrowstyle="<->",color=INK2,lw=0.6,mutation_scale=5),zorder=5)
ax.text((X(p1["end"])+X(p2["start"]))/2,0.50,"1,100 bp",ha="center",fontsize=4.7,color=INK2)
ax.text(0.5,0.08,"$O.$ $sativa$ chr 11 : 27,978,368–27,989,128",ha="center",fontsize=4.8,color=MUTED)
ax.text(0.5,0.86,"The $Pik$ sensor–helper pair",ha="center",fontsize=5.6,color=INK)
ax.set_xlim(-0.07,1.07); ax.set_ylim(0,1)
panel(ax,"a",dx=-0.065,dy=1.04)

# b: sensor vs helper orthogroup breadth
ax=fig.add_subplot(gs[0,4:8]); ygrid(ax,axis="x")
hog=[og for og in pan if any(m.split("|")[1]=="Os11g0689100" for m in og["members"])][0]
sog=[og for og in pan if any(m.split("|")[1]=="Os11g0688832" for m in og["members"])][0]
ax.barh([1],[hog["nspecies"]],height=0.5,color="#0072B2",zorder=3)
ax.barh([0],[sog["nspecies"]],height=0.5,color=HMA_C,zorder=3)
ax.set_yticks([0,1]); ax.set_yticklabels(["$Pik$-1 sensor\northogroup","$Pik$-2 helper\northogroup"],fontsize=5.2)
ax.set_xlabel("Genomes containing an orthologue (of 11)")
ax.set_xlim(0,11.6); ax.set_xticks(range(0,12,2))
ax.text(hog["nspecies"]-0.25,1,f"{hog['nspecies']} species · {hog['n']} genes",va="center",ha="right",fontsize=4.9,color="white")
ax.text(sog["nspecies"]+0.25,0,"1 species · 1 gene — rice only",va="center",fontsize=4.9,color=HMA_C)
ax.text(0.5,1.085,"Conserved helper, lineage-restricted sensor",transform=ax.transAxes,
        ha="center",fontsize=5.2,color=INK)
panel(ax,"b",dx=-0.33,dy=1.20)

# c: HMA detection benchmark
ax=fig.add_subplot(gs[0,8:12]); ygrid(ax,axis="x")
labs=["Pfam HMA (PF00403)\ngathering threshold","Pfam HMA (PF00403)\nE ≤ 10","Pikm1/RGA5 HMA seed\n(this study)"]
vals=[0,1,1]; tp=[0,0,1]
y=np.arange(3)
ax.barh(y,vals,height=0.52,color=["#d9d9d9","#d9d9d9",HMA_C],zorder=3)
for i in range(3):
    t="rice $Pik$-1 recovered" if tp[i] else ("0 hits" if vals[i]==0 else "1 hit, not $Pik$-1")
    ax.text(vals[i]+0.06,i,t,va="center",fontsize=4.9,color=INK if tp[i] else MUTED)
ax.set_yticks(y); ax.set_yticklabels(labs,fontsize=4.8); ax.invert_yaxis()
ax.set_xlabel("HMA-integrated NLRs detected (n)"); ax.set_xlim(0,2.4); ax.set_xticks([0,1,2])
ax.text(0.5,1.055,"Benchmark on 3,341 NLRs from 11 genomes",transform=ax.transAxes,
        ha="center",fontsize=4.9,color=MUTED)
panel(ax,"c",dx=-0.42,dy=1.20)

# d: genome-level HMA loci vs distance to nearest NLR (finger millet)
ax=fig.add_subplot(gs[1,0:5]); ygrid(ax,axis="x")
dist=[]
import json as _j
rowsf=[r for r in rows if r["species"]=="Ecoracana"]
byacc={}
for r in rowsf: byacc.setdefault(r["chrom"],[]).append((r["start"],r["end"]))
for l in hma_loci:
    best=None
    for s,e in byacc.get(l["acc"],[]):
        dd=0 if not(l["e"]<s or l["s"]>e) else (s-l["e"] if l["e"]<s else l["s"]-e)
        if best is None or dd<best: best=dd
    if best is not None: dist.append(best/1000)
dist=sorted(dist)
ax.barh(np.arange(len(dist)),dist,height=0.66,color="#9ab8d4",zorder=3)
ax.axvline(10,color=HMA_C,lw=0.9,ls=(0,(3,2)),zorder=4)
ax.text(10,len(dist)-0.2,"  10 kb",fontsize=4.9,color=HMA_C,va="top")
ax.set_xscale("log"); ax.set_xlabel("Distance from HMA locus to nearest NLR gene (kb, log scale)")
ax.set_ylabel("HMA loci in the $E.$ $coracana$ genome")
ax.set_yticks([]); ax.set_xlim(8,4000)
ax.text(0.985,0.06,f"n = {len(dist)} HMA loci (tblastn, bit ≥ 50)\nnone inside or within 10 kb of an NLR",
        transform=ax.transAxes,ha="right",fontsize=4.9,color=INK)
panel(ax,"d",dx=-0.115)

# e: architectures
ax=fig.add_subplot(gs[1,5:12]); ax.axis("off")
order=[("Os11t0688832-01","$O. sativa$  $Pik$-1",True),
       ("TKV99356","$S. viridis$  ortholog",False),
       ("scaffold389.6.cds","$E. crus-galli$  ortholog",False),
       ("EES09001","$S. bicolor$  ortholog",False)]
maxlen=max(pik[k]["length"] for k,_,_ in order); H=0.40; y=len(order)
for k,lab,ish in order:
    r=pik[k]; y-=1; W=r["length"]/maxlen
    ax.plot([0,W],[y,y],color="#c8c8c8",lw=1.0,zorder=1)
    drawn=[]
    for fam,a0,b0 in r["domains"]:
        c=dcat(fam)
        if c is None: continue
        a,b=a0/maxlen,b0/maxlen
        if any(not(b<=p or a>=q) for p,q in drawn): continue
        drawn.append((a,b))
        ax.add_patch(plt.Rectangle((a,y-H/2),max(b-a,0.004),H,facecolor=DC[c],edgecolor="white",lw=0.4,zorder=3))
    if ish:
        a,b=237/maxlen,287/maxlen
        ax.add_patch(plt.Rectangle((a,y-H/2),b-a,H,facecolor=HMA_C,edgecolor="white",lw=0.4,zorder=4))
        ax.text((a+b)/2,y+H/2+0.14,"HMA",ha="center",fontsize=4.9,color=HMA_C,fontweight="bold")
    else:
        ax.text(262/maxlen,y,"×",ha="center",va="center",fontsize=6.5,color="#9a9a9a",zorder=5)
    ax.text(-0.012,y,lab,ha="right",va="center",fontsize=5.0,color=INK)
    ax.text(W+0.012,y,f"{r['length']} aa",ha="left",va="center",fontsize=4.8,color=INK2)
ax.axvline(262/maxlen,color="#dcdcdc",lw=0.5,ls=(0,(2,2)),zorder=0)
ax.set_xlim(-0.34,1.16); ax.set_ylim(-0.80,len(order)-0.15)
hs=[plt.Rectangle((0,0),1,1,facecolor=DC[k]) for k in ["CC","NB","WHD","LRR","HMA"]]
ax.legend(hs,["CC (Rx_N)","NB-ARC","WHD / ARC2","LRR","integrated HMA"],
          loc="lower center",bbox_to_anchor=(0.40,-0.16),ncol=5,fontsize=5.0)
ax.text(0.40,1.03,"× marks the position of the HMA integration in rice $Pik$-1",
        transform=ax.transAxes,ha="center",fontsize=4.8,color=MUTED)
panel(ax,"e",dx=-0.015,dy=1.10)

fig.savefig(f"{ROOT}/figures/Figure6.pdf"); fig.savefig(f"{ROOT}/figures/Figure6.png")
print(f"Figure 6 saved | helper OG {hog['og']}: {hog['nspecies']} spp, {hog['n']} genes | sensor OG {sog['og']}: {sog['nspecies']} spp, {sog['n']} genes")
print(f"finger millet HMA loci: n={len(dist)}  min distance to NLR = {min(dist):.1f} kb  median = {np.median(dist):.0f} kb")
