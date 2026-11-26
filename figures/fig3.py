import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, random, json
from collections import Counter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; pan=d["pan"]
random.seed(1)
ORD=ORD11
CATS=["core","softcore","shell","species-specific"]
CCOL={"core":"#332288","softcore":"#0072B2","shell":"#56B4E9","species-specific":"#D55E00"}

fig=plt.figure(figsize=(W2,W2*0.50))
gs=fig.add_gridspec(2,12,hspace=0.70,wspace=3.0,left=0.075,right=0.985,top=0.92,bottom=0.12)

# a: composition (OGs vs genes)
ax=fig.add_subplot(gs[0,0:3]); ygrid(ax)
og_n=[sum(1 for r in pan if r["cat"]==c) for c in CATS]
gn=[sum(r["n"] for r in pan if r["cat"]==c) for c in CATS]
for k,(vals,lab) in enumerate([(og_n,"orthogroups"),(gn,"NLR genes")]):
    tot=sum(vals); bot=0
    for c,v in zip(CATS,vals):
        ax.bar([k],[v/tot*100],bottom=[bot],width=0.6,color=CCOL[c],zorder=3,
               edgecolor="white",linewidth=0.6)
        if v/tot*100>5:
            ax.text(k,bot+v/tot*50,f"{v:,}",ha="center",va="center",fontsize=4.8,
                    color="white" if c in("core","softcore","species-specific") else INK)
        bot+=v/tot*100
ax.set_xticks([0,1]); ax.set_xticklabels(["orthogroups\n(n = %d)"%sum(og_n),"NLR genes\n(n = %d)"%sum(gn)],fontsize=5.2)
ax.set_ylabel("Orthogroup composition (%)"); ax.set_ylim(0,100); ax.set_xlim(-0.6,1.6)
hs=[plt.Rectangle((0,0),1,1,facecolor=CCOL[c]) for c in CATS]
ax.legend(hs,["core (11/11)","soft-core (10/11)","shell (2–9)","species-specific (1)"],
          loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2,fontsize=5.0)
panel(ax,"a",dx=-0.30,dy=1.17)

# b: accumulation curve
ax=fig.add_subplot(gs[0,3:7]); ygrid(ax)
sets=[{r["og"] for r in pan if s in r["species"]} for s in ORD]
PAN=[];CORE=[]
for _ in range(200):
    idx=list(range(len(ORD))); random.shuffle(idx)
    u=set(); i_=None; pv=[]; cv=[]
    for k,j in enumerate(idx):
        u|=sets[j]; i_=sets[j] if i_ is None else (i_&sets[j])
        pv.append(len(u)); cv.append(len(i_))
    PAN.append(pv); CORE.append(cv)
PAN=np.array(PAN); CORE=np.array(CORE); xs=np.arange(1,len(ORD)+1)
for arr,c,lab in [(PAN,"#D55E00","cumulative orthogroups"),(CORE,"#332288","core orthogroups")]:
    m=arr.mean(0); lo=np.percentile(arr,5,axis=0); hi=np.percentile(arr,95,axis=0)
    ax.fill_between(xs,lo,hi,color=c,alpha=0.18,lw=0,zorder=2)
    ax.plot(xs,m,color=c,lw=1.3,zorder=3,label=lab)
    ax.scatter(xs,m,s=9,color=c,zorder=4,edgecolor="white",linewidth=0.5)
ax.set_xlabel("Genomes sampled (n)"); ax.set_ylabel("Orthogroups (n)")
ax.set_xticks(xs); ax.legend(loc="center left",fontsize=5.2)
ax.text(len(ORD),PAN.mean(0)[-1]+28,f"{int(PAN.mean(0)[-1])}",ha="right",fontsize=5,color="#D55E00")
ax.text(len(ORD),CORE.mean(0)[-1]+28,f"{int(CORE.mean(0)[-1])}",ha="right",fontsize=5,color="#332288")
panel(ax,"b",dx=-0.22)

# c: OG presence heatmap (top 35 by size)
ax=fig.add_subplot(gs[0,7:12])
top=[r for r in pan if r["nspecies"]>=2][:35]
M=np.array([[sum(1 for m in r["members"] if m.split("|")[0]==s) for s in ORD] for r in top],dtype=float)
im=ax.imshow(np.log10(M+1),cmap="Blues",aspect="auto")
ax.set_xticks(np.arange(len(ORD))); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_yticks([0,9,19,29,34]); ax.set_yticklabels([top[i]["og"] for i in [0,9,19,29,34]],fontsize=4.6)
ax.set_ylabel("Orthogroup (ranked by size)",fontsize=5.5,labelpad=1)
for s_ in ax.spines.values(): s_.set_visible(False)
cax=ax.inset_axes([0.65,1.035,0.34,0.032])
cb=fig.colorbar(im,cax=cax,orientation="horizontal")
cb.set_label("NLR copies per genome, log$_{10}$(n+1)",fontsize=4.6,labelpad=1.5)
cax.xaxis.set_label_position("top"); cax.xaxis.set_ticks_position("top")
cb.ax.tick_params(labelsize=4.3,width=0.4,length=1.5,pad=1); cb.outline.set_linewidth(0)
panel(ax,"c",dx=-0.12,dy=1.15)

# d: per-species OG stats
ax=fig.add_subplot(gs[1,0:6]); ygrid(ax)
x=np.arange(len(ORD)); bot=np.zeros(len(ORD))
for c in CATS:
    v=np.array([sum(sum(1 for m in r["members"] if m.split("|")[0]==s)
                    for r in pan if r["cat"]==c and s in r["species"])
                for s in ORD],dtype=float)
    ax.bar(x,v,bottom=bot,width=0.68,color=CCOL[c],zorder=3,edgecolor="white",linewidth=0.5)
    bot+=v
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_ylabel("NLR genes (n)")
ax.set_ylim(0,470)
panel(ax,"d",dx=-0.13)

# e: species-specific fraction
ax=fig.add_subplot(gs[1,6:12]); ygrid(ax)
ss=[sum(r["n"] for r in pan if r["cat"]=="species-specific" and r["species"]==[s]) for s in ORD]
tot=[sum(1 for r in rows if r["species"]==s) for s in ORD]
pc=[100*a/b for a,b in zip(ss,tot)]
ax.bar(x,pc,width=0.68,color=[COL[s] for s in ORD],zorder=3)
for i in range(len(ORD)): ax.text(i,pc[i]+0.9,f"{ss[i]}",ha="center",va="bottom",fontsize=4.8,color=INK)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_ylabel("Species-specific NLRs (%)"); ax.set_ylim(0,max(pc)*1.22)
ax.text(0.015,1.0,"value above bar = number of species-specific NLRs",transform=ax.transAxes,
        fontsize=4.6,color=MUTED,va="top")
panel(ax,"e",dx=-0.13)

fig.savefig(f"{ROOT}/figures/Figure3.pdf"); fig.savefig(f"{ROOT}/figures/Figure3.png")
print("Figure 3 saved")
print("pan mean final:",int(PAN.mean(0)[-1]),"core mean final:",int(CORE.mean(0)[-1]))
