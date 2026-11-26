import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt
from collections import Counter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; cl=d["clusters"]; sz=d["sizes"]

ORD=ORD11
cnt={s:sum(1 for r in rows if r["species"]==s) for s in ORD}
dens={s:cnt[s]/(sz[s]/1e6)*100 for s in ORD}

fig=plt.figure(figsize=(W2,W2*0.52))
gs=fig.add_gridspec(2,12,hspace=0.72,wspace=2.6,left=0.075,right=0.985,top=0.91,bottom=0.11)

# ---- a: NLR count ----
ax=fig.add_subplot(gs[0,0:5]); ygrid(ax)
x=np.arange(len(ORD))
ax.bar(x,[cnt[s] for s in ORD],width=0.68,color=[COL[s] for s in ORD],zorder=3)
for i,s in enumerate(ORD):
    ax.text(i,cnt[s]+8,str(cnt[s]),ha="center",va="bottom",fontsize=5,color=INK)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=38,ha="right",style="italic")
ax.set_ylabel("NLR genes (n)"); ax.set_ylim(0,510)
panel(ax,"a",dx=-0.155)

# ---- b: architecture composition ----
ax=fig.add_subplot(gs[0,5:9]); ygrid(ax)
bot=np.zeros(len(ORD))
for a in ARCH_ORDER:
    v=np.array([sum(1 for r in rows if r["species"]==s and archgroup(r["arch"])==a)/cnt[s]*100 for s in ORD])
    ax.bar(x,v,bottom=bot,width=0.68,color=ARCH_COL[a],label=a,zorder=3,
           edgecolor="white",linewidth=0.5)
    bot+=v
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=38,ha="right",style="italic")
ax.set_ylabel("Architecture (% of NLRs)"); ax.set_ylim(0,100)
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.03),ncol=5,fontsize=5.2)
panel(ax,"b",dx=-0.175,dy=1.115)

# ---- c: genome size vs NLR density ----
ax=fig.add_subplot(gs[0,9:12]); ygrid(ax); ax.grid(axis="x",color=GRID,lw=0.4,zorder=0)
for s in ORD:
    ax.scatter(sz[s]/1e6,dens[s],s=26,color=COL[s],zorder=3,
               edgecolor="white",linewidth=0.6)
ax.set_xscale("log"); ax.set_xlim(210,3000); ax.set_ylim(-6,128)
ax.set_xticks([250,500,1000,2000]); ax.set_xticklabels(["250","500","1,000","2,000"])
ax.minorticks_off()
ax.set_xlabel("Assembly size (Mb, log scale)"); ax.set_ylabel("NLRs per 100 Mb")
for s_,dxy,ha in [("Bdistachyon",(6,-2),"left"),("Ecoracana",(0,-11),"center"),("Zmays",(0,7),"center")]:
    ax.annotate(LAB[s_],(sz[s_]/1e6,dens[s_]),xytext=dxy,textcoords="offset points",
                fontsize=4.8,style="italic",color=INK,ha=ha)
import numpy as _np
r=_np.corrcoef(_np.log10([sz[s_]/1e6 for s_ in ORD]),[dens[s_] for s_ in ORD])[0,1]
ax.text(0.97,0.95,f"$r$ = {r:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=5,color=MUTED)
panel(ax,"c",dx=-0.33)

# ---- d: clustering ----
ax=fig.add_subplot(gs[1,0:5]); ygrid(ax)
pct=[cl[s]["pct"] for s in ORD]
ax.bar(x,pct,width=0.68,color=[COL[s] for s in ORD],zorder=3)
for i,s in enumerate(ORD):
    ax.text(i,pct[i]+1.4,f"{cl[s]['g']}",ha="center",va="bottom",fontsize=4.6,color=MUTED)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=38,ha="right",style="italic")
ax.set_ylabel("NLRs in physical clusters (%)"); ax.set_ylim(0,97)
ax.text(0.015,1.00,"value above bar = cluster distance threshold $g$ (kb)",transform=ax.transAxes,
        fontsize=4.6,color=MUTED,va="top")
panel(ax,"d",dx=-0.155)

# ---- e: cluster size distribution ----
ax=fig.add_subplot(gs[1,5:9]); ygrid(ax)
allsz=[]
for s in ORD: allsz += [n for _,n in cl[s]["clusters"]]
bins=np.arange(1.5,16.5,1)
h,_=np.histogram(allsz,bins=bins)
ax.bar(bins[:-1]+0.5,h,width=0.7,color="#4d4d4d",zorder=3)
big=sum(1 for v in allsz if v>15)
if big: ax.bar([16.5],[big],width=0.7,color="#D55E00",zorder=3)
ax.set_xlabel("NLRs per cluster"); ax.set_ylabel("Clusters (n)")
ax.set_xticks(list(range(2,16,2))+[16.5]); ax.set_xticklabels([str(v) for v in range(2,16,2)]+[">15"])
panel(ax,"e",dx=-0.175)

# ---- f: singleton vs clustered ----
ax=fig.add_subplot(gs[1,9:12]); ygrid(ax)
clus=np.array([cl[s]["in_cluster"] for s in ORD]); sing=np.array([cnt[s] for s in ORD])-clus
ax.barh(x,clus,height=0.66,color="#4d4d4d",label="clustered",zorder=3)
ax.barh(x,sing,left=clus,height=0.66,color="#d9d9d9",label="singleton",zorder=3,
        edgecolor="white",linewidth=0.5)
ax.set_yticks(x); ax.set_yticklabels([LAB[s] for s in ORD],style="italic",fontsize=5)
ax.invert_yaxis(); ax.set_xlabel("NLR genes (n)")
ax.legend(loc="lower right",fontsize=5.2)
panel(ax,"f",dx=-0.46)

fig.savefig(f"{ROOT}/figures/Figure1.pdf"); fig.savefig(f"{ROOT}/figures/Figure1.png")
print("Figure 1 saved")
print("counts:",{s:cnt[s] for s in ORD})
print("density per 100Mb:",{s:round(dens[s],1) for s in ORD})
