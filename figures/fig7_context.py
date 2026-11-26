import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
ORD=ORD11
gd=json.load(open(f"{ROOT}/nlr/gene_density.json"))
te=json.load(open(f"{ROOT}/genome/te_enrichment10.json"))
sy=json.load(open(f"{ROOT}/genome/synteny_ec_si.json"))
tr=json.load(open(f"{ROOT}/genome/truncation_test.json"))
ACC="#D55E00"; NEU="#9ab8d4"

fig=plt.figure(figsize=(W2,W2*0.32))
gs=fig.add_gridspec(1,13,wspace=3.4,left=0.075,right=0.985,top=0.86,bottom=0.26)

# a: gene density
ax=fig.add_subplot(gs[0,0:5]); ygrid(ax)
x=np.arange(len(ORD))
ax.bar(x-0.19,[gd[s]["nlr"] for s in ORD],width=0.36,color=ACC,label="NLR neighbourhoods",zorder=3)
ax.bar(x+0.19,[gd[s]["rand"] for s in ORD],width=0.36,color=NEU,label="matched random genes",zorder=3)
for i,s in enumerate(ORD):
    if gd[s]["p"]<0.001: ax.text(i,max(gd[s]["nlr"],gd[s]["rand"])+0.7,"***",ha="center",fontsize=4.6,color=MUTED)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=42,ha="right",style="italic")
ax.set_ylabel("Genes per 200 kb window (median)")
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=5.2)
panel(ax,"a",dx=-0.145,dy=1.17)

# b: TE enrichment
ax=fig.add_subplot(gs[0,5:8]); ygrid(ax)
v=[100*te["nlr_hits"]/te["nlr_flanks"],100*te["ctrl_hits"]/te["ctrl_flanks"]]
ax.bar([0,1],v,width=0.55,color=[ACC,NEU],zorder=3)
for i,val in enumerate(v): ax.text(i,val+1.0,f"{val:.1f}%",ha="center",fontsize=5.2,color=INK)
ax.set_xticks([0,1]); ax.set_xticklabels(["NLR\nflanks","random gene\nflanks"],fontsize=5.2)
ax.set_ylabel("Flanks carrying a TE domain (%)"); ax.set_ylim(0,72)
ax.text(0.5,1.0,f"OR = {te['odds']:.2f},  $P$ = {te['p']:.0e}",transform=ax.transAxes,
        ha="center",fontsize=5.0,color=MUTED)
ax.text(0.5,-0.42,"$E.$ $coracana$, 10 kb flanks",transform=ax.transAxes,ha="center",
        fontsize=4.8,color=MUTED)
panel(ax,"b",dx=-0.33,dy=1.17)

# c: synteny depletion
ax=fig.add_subplot(gs[0,8:10]); ygrid(ax)
v=[100*sy["nlr_in_block"]/sy["nlr_anchored"],100*sy["other_in"]/sy["other"]]
ax.bar([0,1],v,width=0.55,color=[ACC,NEU],zorder=3)
for i,val in enumerate(v): ax.text(i,val+1.5,f"{val:.0f}%",ha="center",fontsize=5.2,color=INK)
ax.set_xticks([0,1]); ax.set_xticklabels(["NLR","other\ngenes"],fontsize=5.2)
ax.set_ylabel("Inside a collinear block (%)"); ax.set_ylim(0,100)
ax.text(0.5,1.0,f"OR = {sy['odds']:.2f}\n$P$ = {sy['p']:.0e}",transform=ax.transAxes,
        ha="center",fontsize=5.0,color=MUTED)
ax.text(0.5,-0.42,"$E.$ $coracana$ vs $S.$ $italica$",transform=ax.transAxes,ha="center",
        fontsize=4.8,color=MUTED)
panel(ax,"c",dx=-0.48,dy=1.17)

# d: NLR-Annotator re-annotation of the same loci
ax=fig.add_subplot(gs[0,10:13]); ygrid(ax)
import json as _j
from collections import Counter as _C, defaultdict as _dd
cmp=_j.load(open(f"{ROOT}/genome/nlra_compare.json"))
tab=_dd(_C)
for g,ann,na,nm,n in cmp:
    if na: tab[ann.rstrip("X")][na]+=1
order=["N","CN","NL","CNL"]
COLS={"N":"#d9d9d9","CN":"#b0b0b0","NL":"#7fb3d5","CNL":"#0072B2"}
xx=np.arange(len(order)); bot=np.zeros(len(order))
for call in ["N","CN","NL","CNL"]:
    v=np.array([100*tab[a].get(call,0)/max(sum(tab[a].values()),1) for a in order])
    ax.bar(xx,v,bottom=bot,width=0.62,color=COLS[call],label=call,zorder=3,
           edgecolor="white",linewidth=0.5)
    bot+=v
ax.set_xticks(xx); ax.set_xticklabels([f"{a}\n(n={sum(tab[a].values())})" for a in order],fontsize=5.0)
ax.set_ylabel("NLR-Annotator call (%)"); ax.set_ylim(0,100)
ax.set_xlabel("Architecture in the annotation",labelpad=1,fontsize=5.5)
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=4,fontsize=5.0,title=None)
gain=sum(1 for g,a,na,nm,n in cmp if na and "L" not in a.rstrip("X") and "L" in na)
nol=sum(1 for g,a,na,nm,n in cmp if na and "L" not in a.rstrip("X"))
ax.text(0.5,-0.42,f"{gain}/{nol} ({100*gain/nol:.0f}%) of LRR-less models gain an LRR",
        transform=ax.transAxes,ha="center",fontsize=4.8,color=MUTED)
panel(ax,"d",dx=-0.33,dy=1.22)

fig.savefig(f"{ROOT}/figures/Figure7.pdf"); fig.savefig(f"{ROOT}/figures/Figure7.png")
print("Figure 7 saved")
print("NLR-Annotator comparison panel written")
