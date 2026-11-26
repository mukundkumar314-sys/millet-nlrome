import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
from collections import defaultdict
from scipy import stats
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
ORD=ORD11
ACC="#D55E00"; NEU="#9ab8d4"
pd_=json.load(open(f"{ROOT}/nlr/robust_pairdist.json"))
kg=json.load(open(f"{ROOT}/nlr/kaks_groups.json"))
pergene=json.load(open(f"{ROOT}/nlr/kaks_gene.json"))
byid={r["id"]:r for r in rows}

fig=plt.figure(figsize=(W2,W2*0.31))
gs=fig.add_gridspec(1,14,wspace=3.6,left=0.070,right=0.985,top=0.86,bottom=0.26)

# a: ID-NLR rate per species
ax=fig.add_subplot(gs[0,0:5]); ygrid(ax)
x=np.arange(len(ORD))
rate=[100*sum(1 for r in rows if r["species"]==s and r["n_id"]>0)/sum(1 for r in rows if r["species"]==s) for s in ORD]
cols=[ACC if s=="Ecoracana" else NEU for s in ORD]
ax.bar(x,rate,width=0.66,color=cols,zorder=3)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=42,ha="right",style="italic")
ax.set_ylabel("NLRs carrying an ID (%)"); ax.set_ylim(0,13.4)
i=ORD.index("Ecoracana")
ax.annotate("lowest of 11\nOR 0.23, $P$ = 8×10$^{-4}$",(i,rate[i]),xytext=(0,22),
            textcoords="offset points",ha="center",fontsize=4.8,color=ACC,
            arrowprops=dict(arrowstyle="-",color=ACC,lw=0.5))
panel(ax,"a",dx=-0.145)

# b: pair-distance robustness
ax=fig.add_subplot(gs[0,5:8]); ygrid(ax)
ks=sorted(int(k) for k in pd_)
xs=np.arange(len(ks))
od=[pd_[str(k)]["od"] for k in ks]
ax.bar(xs,od,width=0.60,color=ACC,zorder=3)
ax.axhline(1,color=MUTED,lw=0.6,ls=(0,(3,2)),zorder=2)
for i,k in enumerate(ks):
    ax.text(i,od[i]+0.08,f"{pd_[str(k)]['p']:.0e}".replace("e-0","e-"),ha="center",fontsize=4.2,color=MUTED)
ax.set_xticks(xs); ax.set_xticklabels([f"{k//1000}" for k in ks])
ax.set_xlabel("Pair distance cut-off (kb)"); ax.set_ylabel("ID enrichment in pairs (odds ratio)")
ax.set_ylim(0,4.1)
ax.text(0.5,-0.42,"$P$ value above each bar",transform=ax.transAxes,ha="center",fontsize=4.6,color=MUTED)
panel(ax,"b",dx=-0.34)

# c: Ka/Ks by orthogroup class
ax=fig.add_subplot(gs[0,8:11]); ygrid(ax)
data=[kg["core"],kg["shell"],kg["spec"]]
labs=["core +\nsoft-core","shell","species-\nspecific"]
bp=ax.boxplot(data,widths=0.52,patch_artist=True,showfliers=False,
              medianprops=dict(color=INK,lw=0.8),whiskerprops=dict(lw=0.5,color=INK2),
              capprops=dict(lw=0.5,color=INK2),boxprops=dict(lw=0.4,edgecolor=INK2))
for b,c in zip(bp["boxes"],["#332288","#56B4E9",ACC]): b.set_facecolor(c); b.set_alpha(0.85)
for i,v in enumerate(data):
    ax.scatter(np.random.normal(i+1,0.055,len(v)),v,s=1.1,color=INK2,alpha=0.30,zorder=4,linewidths=0)
ax.set_xticks([1,2,3]); ax.set_xticklabels(labs,fontsize=4.9)
ax.set_ylabel("Ka/Ks (orthogroup median)"); ax.set_ylim(0,1.25)
ax.text(0.5,1.0,"$P$ = 2×10$^{-5}$      $P$ = 2×10$^{-6}$",transform=ax.transAxes,
        ha="center",fontsize=4.6,color=MUTED)
panel(ax,"c",dx=-0.34)

# d: Ka/Ks ID vs non-ID
ax=fig.add_subplot(gs[0,11:14]); ygrid(ax)
idv=[v for k,v in pergene.items() if byid[k]["n_id"]>0]
nov=[v for k,v in pergene.items() if byid[k]["n_id"]==0]
bp=ax.boxplot([nov,idv],widths=0.50,patch_artist=True,showfliers=False,
              medianprops=dict(color=INK,lw=0.8),whiskerprops=dict(lw=0.5,color=INK2),
              capprops=dict(lw=0.5,color=INK2),boxprops=dict(lw=0.4,edgecolor=INK2))
bp["boxes"][0].set_facecolor(NEU); bp["boxes"][1].set_facecolor(ACC)
for i,v in enumerate([nov,idv]):
    ax.scatter(np.random.normal(i+1,0.05,len(v)),v,s=0.8,color=INK2,alpha=0.18,zorder=4,linewidths=0)
ax.set_xticks([1,2]); ax.set_xticklabels([f"no ID\n(n={len(nov)})",f"ID-NLR\n(n={len(idv)})"],fontsize=4.9)
ax.set_ylabel("Ka/Ks (gene mean)"); ax.set_ylim(0,1.25)
ax.text(0.5,1.0,"$P$ = 2×10$^{-10}$",transform=ax.transAxes,ha="center",fontsize=4.8,color=MUTED)
panel(ax,"d",dx=-0.34)

fig.savefig(f"{ROOT}/figures/Figure10.pdf"); fig.savefig(f"{ROOT}/figures/Figure10.png")
print("Figure 10 saved")
