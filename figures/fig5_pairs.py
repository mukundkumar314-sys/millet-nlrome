import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
from collections import Counter, defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
ORD=ORD11
pairs=json.load(open(f"{ROOT}/nlr/h2h_pairs.json"))
by={r["gene"]+"|"+r["species"]:r for r in rows}
DEC="#D55E00"

fig=plt.figure(figsize=(W2,W2*0.62))
gs=fig.add_gridspec(2,12,height_ratios=[0.86,1.0],hspace=0.62,wspace=3.2,
                    left=0.080,right=0.980,top=0.925,bottom=0.080)

# a: orientation of neighbouring NLR pairs
ax=fig.add_subplot(gs[0,0:5]); ygrid(ax)
x=np.arange(len(ORD))
orient={}
for sp in ORD:
    sub=[r for r in rows if r["species"]==sp]
    byc=defaultdict(list)
    for r in sub: byc[str(r.get("chrom_lab",r["chrom"]))].append(r)
    h=t=tan=0
    for c,lst in byc.items():
        lst.sort(key=lambda r:r["start"])
        for i in range(len(lst)-1):
            a,b=lst[i],lst[i+1]; gp=b["start"]-a["end"]
            if gp<0 or gp>10000: continue
            if a["strand"]=="-" and b["strand"]=="+": h+=1
            elif a["strand"]=="+" and b["strand"]=="-": t+=1
            else: tan+=1
    orient[sp]=(h,t,tan)
bot=np.zeros(len(ORD))
for k,(lab,c) in enumerate([("head-to-head",DEC),("tail-to-tail","#999999"),("tandem","#d9d9d9")]):
    v=np.array([orient[s][k] for s in ORD],dtype=float)
    ax.bar(x,v,bottom=bot,width=0.68,color=c,label=lab,zorder=3,edgecolor="white",linewidth=0.5)
    bot+=v
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=42,ha="right",style="italic")
ax.set_ylabel("Neighbouring NLR pairs ≤ 10 kb (n)")
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=3,fontsize=5.2)
panel(ax,"a",dx=-0.155,dy=1.14)

# b: ID-carrying head-to-head pairs
ax=fig.add_subplot(gs[0,5:9]); ygrid(ax)
h2h=[orient[s][0] for s in ORD]
idp=[]
for sp in ORD:
    n=0
    for p in pairs:
        if p["species"]!=sp: continue
        a=by.get(p["g1"]+"|"+sp); b=by.get(p["g2"]+"|"+sp)
        if a and b and (a["n_id"]>0 or b["n_id"]>0): n+=1
    idp.append(n)
ax.bar(x,h2h,width=0.68,color="#d9d9d9",zorder=3,label="no integrated domain")
ax.bar(x,idp,width=0.68,color=DEC,zorder=4,label="sensor carries an ID")
for i,s in enumerate(ORD):
    if idp[i]==0: ax.text(i,h2h[i]+0.35,"0",ha="center",va="bottom",fontsize=5.4,color=DEC,fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=42,ha="right",style="italic")
ax.set_ylabel("Head-to-head NLR pairs (n)")
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2,fontsize=5.2)
panel(ax,"b",dx=-0.185,dy=1.14)

# c: which decoy domain, by species
ax=fig.add_subplot(gs[0,9:12])
dec=defaultdict(Counter)
for p in pairs:
    for g in (p["g1"],p["g2"]):
        r=by.get(g+"|"+p["species"])
        if r and r["n_id"]>0:
            for i in r["ids"].split(";"):
                if i: dec[i][p["species"]]+=1
MERGE={"PK_Tyr_Ser-Thr":"protein kinase","Pkinase":"protein kinase",
       "Myb_DNA-binding":"Myb","Myb_DNA-bind_6":"Myb"}
dec2=defaultdict(Counter)
for k,v in dec.items():
    for s,n in v.items(): dec2[MERGE.get(k,k)][s]=max(dec2[MERGE.get(k,k)][s],n)
fams=sorted(dec2,key=lambda k:-sum(dec2[k].values()))[:9]
M=np.array([[dec2[f][s] for s in ORD] for f in fams],dtype=float)
im=ax.imshow(M,cmap="Oranges",aspect="auto",vmin=0,vmax=max(M.max(),1))
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=42,ha="right",style="italic",fontsize=4.8)
ax.set_yticks(np.arange(len(fams))); ax.set_yticklabels(fams,fontsize=4.9)
for i in range(len(fams)):
    for j in range(len(ORD)):
        if M[i,j]>0: ax.text(j,i,int(M[i,j]),ha="center",va="center",fontsize=4.4,
                             color="white" if M[i,j]>M.max()*0.6 else INK)
ax.set_xticks(np.arange(-.5,len(ORD),1),minor=True); ax.set_yticks(np.arange(-.5,len(fams),1),minor=True)
ax.grid(which="minor",color="white",lw=0.7); ax.tick_params(which="minor",length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.text(0.5,1.055,"Decoy domain carried by the paired sensor",transform=ax.transAxes,
        ha="center",fontsize=5.0,color=INK)
panel(ax,"c",dx=-0.33,dy=1.16)

# d: schematics of three pair types
ax=fig.add_subplot(gs[1,0:12]); ax.axis("off")
EX=[("O. sativa  chr 11",  "$Pik$-1   CNL + HMA","$Pik$-2   CNL",1100,"HMA"),
    ("S. bicolor  chr 8",  "SORBI_3008G174100   CNL + WRKY","SORBI_3008G174200   CNL",2223,"WRKY"),
    ("P. miliaceum  chr 17","C2845_PM17G03140   NL + AvrRpt-cleav.","C2845_PM17G03150   CN",5390,"AvrRpt-cleav."),
    ("E. coracana  chr 6B","QOZ80_6BG0483700   CNL","QOZ80_6BG0483710   CNL",8399,None)]
NLRC="#c3d4e6"; L=0.26; G=0.075; X0=0.10
y=len(EX)
ax.text(X0+L/2,len(EX)-0.30,"sensor",ha="center",fontsize=5.6,color=INK)
ax.text(X0+L+G+L/2,len(EX)-0.30,"helper",ha="center",fontsize=5.6,color=INK)
for loc,s1,s2,gap,dom in EX:
    y-=1
    ax.add_patch(plt.Rectangle((X0,y-0.15),L,0.30,facecolor=NLRC,edgecolor="white",lw=0.5,zorder=3))
    ax.add_patch(plt.Rectangle((X0+L+G,y-0.15),L,0.30,facecolor=NLRC,edgecolor="white",lw=0.5,zorder=3))
    if dom:
        ax.add_patch(plt.Rectangle((X0+0.020,y-0.15),0.058,0.30,facecolor=DEC,edgecolor="white",lw=0.4,zorder=4))
        ax.text(X0+0.049,y+0.19,dom,ha="center",fontsize=4.5,color=DEC,fontweight="bold")
    else:
        ax.text(X0+0.049,y+0.19,"no decoy",ha="center",fontsize=4.5,color=MUTED,style="italic")
    ax.annotate("",xy=(X0-0.030,y),xytext=(X0-0.002,y),zorder=5,
                arrowprops=dict(arrowstyle="-|>",color=INK2,lw=0.7,mutation_scale=6))
    ax.annotate("",xy=(X0+2*L+G+0.030,y),xytext=(X0+2*L+G+0.002,y),zorder=5,
                arrowprops=dict(arrowstyle="-|>",color=INK2,lw=0.7,mutation_scale=6))
    ax.plot([X0+L,X0+L+G],[y,y],color="#c0c0c0",lw=0.8,zorder=1)
    ax.text(X0+L+G/2,y+0.19,f"{gap:,} bp",ha="center",fontsize=4.4,color=MUTED)
    ax.text(X0-0.042,y,loc,ha="right",va="center",fontsize=5.0,style="italic",color=INK)
    ax.text(X0+L/2,y-0.26,s1,ha="center",va="top",fontsize=4.6,color=INK2)
    ax.text(X0+L+G+L/2,y-0.26,s2,ha="center",va="top",fontsize=4.6,color=INK2)
ax.set_xlim(-0.16,0.78); ax.set_ylim(-0.70,len(EX)-0.10)
hs=[plt.Rectangle((0,0),1,1,facecolor=c) for c in [NLRC,DEC]]
ax.legend(hs,["NLR (CC-NB-ARC-LRR)","integrated decoy domain"],loc="lower center",
          bbox_to_anchor=(0.34,-0.12),ncol=2,fontsize=5.2)
panel(ax,"d",dx=0.0,dy=1.04)

fig.savefig(f"{ROOT}/figures/Figure5.pdf"); fig.savefig(f"{ROOT}/figures/Figure5.png")
tot=sum(orient[s][0] for s in ORD); tid=sum(idp)
print(f"Figure 5 saved | head-to-head pairs {tot} | ID-carrying {tid} | E. coracana ID-pairs {idp[ORD.index('Ecoracana')]}")
print("decoy families:",fams)
