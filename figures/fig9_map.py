import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json, re
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; cl=d["clusters"]
fm=[r for r in rows if r["species"]=="Ecoracana"]
CHL={}
for ln in open(f"{ROOT}/gff/Ecoracana.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\r\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': CHL[f[2]]=int(f[8])
cand=json.load(open(f"{ROOT}/nlr/fm_blastR_candidates.json"))
cmap={}
for c in cand: cmap.setdefault(c["fm_gene"],set()).add(c["rgene"])
pairs=[p for p in json.load(open(f"{ROOT}/nlr/h2h_pairs.json")) if p["species"]=="Ecoracana"]
pg=set()
for p in pairs: pg.add(p["g1"]); pg.add(p["g2"])
g=cl["Ecoracana"]["g"]*1000
CA="#D55E00"; CB="#0072B2"; CAND="#117733"; PAIRC="#332288"

fig=plt.figure(figsize=(W2,W2*0.60))
ax=fig.add_axes([0.055,0.085,0.925,0.845])
P=[str(i) for i in range(1,10)]; maxL=max(CHL.values()); BW=0.26
for i,p in enumerate(P):
    for k,(sub,c) in enumerate([("A",CA),("B",CB)]):
        ch=f"{p}{sub}"; L=CHL[ch]/1e6; y=i+(-0.20 if k==0 else 0.20)
        ax.add_patch(plt.Rectangle((0,y-BW/2),L,BW,facecolor="#f2f2f2",
                                   edgecolor="#c4c4c4",lw=0.35,zorder=2))
        sel=sorted([r for r in fm if str(r.get("chrom_lab"))==ch],key=lambda r:r["start"])
        pos=[(r["start"]+r["end"])//2 for r in sel]
        cur=[pos[0]] if pos else []
        for q in range(1,len(pos)):
            if pos[q]-pos[q-1]<=g: cur.append(pos[q])
            else:
                if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),
                    max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.20,lw=0,zorder=3))
                cur=[pos[q]]
        if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),
            max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.20,lw=0,zorder=3))
        for r in sel:
            x0=((r["start"]+r["end"])//2)/1e6
            ax.plot([x0,x0],[y-BW/2,y+BW/2],color=c,lw=0.42,zorder=4,alpha=0.9)
        for r in sel:
            if r["gene"] in pg:
                x0=((r["start"]+r["end"])//2)/1e6
                ax.plot([x0],[y+BW/2+0.055],marker="v",ms=1.7,color=PAIRC,zorder=6,
                        markeredgewidth=0)
        for r in sel:
            if r["gene"] in cmap:
                x0=((r["start"]+r["end"])//2)/1e6
                ax.plot([x0,x0],[y-BW/2-0.012,y+BW/2+0.012],color=CAND,lw=1.5,zorder=7)
                nm=sorted(cmap[r["gene"]])[0].split("/")[0]
                ax.annotate(nm,(x0,y-BW/2-0.055),ha="center",va="top",fontsize=4.3,
                            color=CAND,style="italic",zorder=8)
        ax.text(-1.6,y,ch,ha="right",va="center",fontsize=5.2,color=INK)
        ax.text(L+1.2,y,str(len(sel)),ha="left",va="center",fontsize=4.7,color=c)
ax.set_xlim(-8,maxL/1e6+8); ax.set_ylim(-0.72,len(P)-0.28)
ax.set_xlabel("Chromosome position (Mb)"); ax.set_yticks([]); ax.invert_yaxis()
for sp_ in ["left","right","top"]: ax.spines[sp_].set_visible(False)
from matplotlib.lines import Line2D
hs=[plt.Rectangle((0,0),1,1,facecolor=CA),plt.Rectangle((0,0),1,1,facecolor=CB),
    plt.Rectangle((0,0),1,1,facecolor="#9a9a9a",alpha=0.25),
    Line2D([],[],color=CAND,lw=1.6),
    Line2D([],[],marker="v",ls="",ms=2.6,color=PAIRC)]
ax.legend(hs,["NLR on A subgenome","NLR on B subgenome","NLR cluster (≥3 genes)",
              "orthologue of a cloned blast R gene","member of a head-to-head pair"],
          loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=5,fontsize=5.0)
ax.text(0.998,0.012,"number at the right of each chromosome = NLR genes on that chromosome",
        transform=ax.transAxes,ha="right",fontsize=4.6,color=MUTED)
fig.savefig(f"{ROOT}/figures/Figure9.pdf"); fig.savefig(f"{ROOT}/figures/Figure9.png")
print("Figure 9 saved; candidates mapped:",len(cmap),"; paired NLRs:",len(pg))
