import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt
from collections import Counter, defaultdict
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
ORD=ORD11
idrows=[r for r in rows if r["n_id"]>0]
idc=Counter(); idsp=defaultdict(set)
for r in idrows:
    for i in r["ids"].split(";"):
        if i: idc[i]+=1; idsp[i].add(r["species"])
SHORT={"PK_Tyr_Ser-Thr":"PK_Tyr_Ser-Thr","C2CH-3rd_BIRD-IDD":"C2CH-3rd_BIRD","C2CH-4th_BIRD-IDD":"C2CH-4th_BIRD",
       "C2H2-2nd_BIRD-IDD":"C2H2-2nd_BIRD"}

fig=plt.figure(figsize=(W2,W2*0.72))
gs=fig.add_gridspec(2,13,height_ratios=[0.80,1.00],hspace=0.30,wspace=4.2,
                    left=0.085,right=0.975,top=0.945,bottom=0.055)

# ---- a ----
ax=fig.add_subplot(gs[0,0:4]); ygrid(ax)
x=np.arange(len(ORD))
n_id=[sum(1 for r in idrows if r["species"]==s) for s in ORD]
n_all=[sum(1 for r in rows if r["species"]==s) for s in ORD]
ax.bar(x,n_id,width=0.68,color=[COL[s] for s in ORD],zorder=3)
for i in range(len(ORD)):
    ax.text(i,n_id[i]+0.9,f"{100*n_id[i]/n_all[i]:.1f}",ha="center",va="bottom",fontsize=4.6,color=MUTED)
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_ylabel("ID-NLR genes (n)"); ax.set_ylim(0,max(n_id)*1.38)
ax.text(0.02,1.005,"value above bar = % of that species' NLRs",transform=ax.transAxes,
        fontsize=4.6,color=MUTED,va="top")
panel(ax,"a",dx=-0.21)

# ---- b ----
top=[k for k,_ in idc.most_common(15)]
M=np.array([[sum(1 for r in idrows if r["species"]==s and fam in r["ids"].split(";"))
             for s in ORD] for fam in top],dtype=float)
ax=fig.add_subplot(gs[0,4:10])
im=ax.imshow(M,cmap="Blues",aspect="auto",vmin=0,vmax=M.max())
ax.set_xticks(np.arange(len(ORD))); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_yticks(np.arange(len(top))); ax.set_yticklabels([SHORT.get(t,t) for t in top],fontsize=4.9)
for i in range(len(top)):
    for j in range(len(ORD)):
        if M[i,j]>0:
            ax.text(j,i,int(M[i,j]),ha="center",va="center",fontsize=4.4,
                    color="white" if M[i,j]>M.max()*0.55 else INK)
ax.set_xticks(np.arange(-.5,len(ORD),1),minor=True); ax.set_yticks(np.arange(-.5,len(top),1),minor=True)
ax.grid(which="minor",color="white",lw=0.7); ax.tick_params(which="minor",length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
cax=ax.inset_axes([0.62,1.035,0.36,0.030])
cb=fig.colorbar(im,cax=cax,orientation="horizontal")
cb.set_label("ID-NLR genes (n)",fontsize=4.8,labelpad=1.5)
cax.xaxis.set_label_position("top"); cax.xaxis.set_ticks_position("top")
cb.ax.tick_params(labelsize=4.4,width=0.4,length=1.5,pad=1); cb.outline.set_linewidth(0)
panel(ax,"b",dx=-0.30,dy=1.14)

# ---- c ----
ax=fig.add_subplot(gs[0,10:13]); ygrid(ax,axis="x")
breadth=Counter(len(v) for v in idsp.values()); ks=sorted(breadth)
ax.barh([str(k) for k in ks],[breadth[k] for k in ks],height=0.66,color="#4d4d4d",zorder=3)
for i,k in enumerate(ks): ax.text(breadth[k]+1.0,i,str(breadth[k]),va="center",fontsize=4.8,color=INK)
ax.set_xlabel("ID families (n)"); ax.set_ylabel("Species sharing the ID family",labelpad=0)
ax.set_xlim(0,max(breadth.values())*1.20)
panel(ax,"c",dx=-0.40)

# ---- d: true-coordinate domain architectures ----
ax=fig.add_subplot(gs[1,0:13]); ax.axis("off")
def pick(fam):
    c=[r for r in idrows if fam in r["ids"].split(";") and r["domains"] and 700<=r["length"]<=1700]
    if not c: c=[r for r in idrows if fam in r["ids"].split(";") and r["domains"]]
    return sorted(c,key=lambda r:r["length"])[len(c)//2] if c else None
WANT=[("Pkinase","kinase"),("WRKY","WRKY"),("Jacalin","Jacalin"),("zf-BED","zf-BED"),
      ("Thioredoxin","thioredoxin"),("PP2","PP2"),("NAM","NAC/NAM"),("WD40","WD40")]
ex=[(lab,pick(f)) for f,lab in WANT]
ex=[(l,r) for l,r in ex if r]
# add the finger-millet dual-ID example
fm=[r for r in idrows if r["species"]=="Ecoracana" and r["n_id"]>=2]
if fm:
    _f=sorted(fm,key=lambda r:r["length"])[0]
    ex.append((" + ".join(_f["ids"].split(";")),_f))
DC={"CC":"#999933","NB":"#0072B2","WHD":"#56B4E9","LRR":"#009E73","ID":"#D55E00"}
def dcat(f):
    if f=="Rx_N": return "CC"
    if f=="NB-ARC": return "NB"
    if f in ("WHD_DRP","WHD_plant_disease"): return "WHD"
    if f.startswith("LRR"): return "LRR"
    if f in ("AAA_22","AAA_16","AAA","AAA_14","NACHT","nSTAND_NTPase5","Sigma54_activat"): return None
    return "ID"
maxlen=max(r["length"] for _,r in ex); H=0.46; y=len(ex)
for lab,r in ex:
    y-=1
    W=r["length"]/maxlen
    ax.plot([0,W],[y,y],color="#c8c8c8",lw=1.0,zorder=1,solid_capstyle="butt")
    drawn=[]
    for dm in r["domains"]:
        c=dcat(dm["fam"])
        if c is None: continue
        a=dm["s"]/maxlen; b=dm["e"]/maxlen
        if any(not(b<=p or a>=q) for p,q in drawn): continue
        drawn.append((a,b))
        ax.add_patch(plt.Rectangle((a,y-H/2),max(b-a,0.004),H,facecolor=DC[c],
                                   edgecolor="white",lw=0.4,zorder=3))
        if c=="ID":
            ax.text((a+b)/2,y+H/2+0.10,dm["fam"],ha="center",va="bottom",fontsize=4.0,color="#D55E00")
    ax.text(-0.012,y,LAB[r["species"]],ha="right",va="center",fontsize=5.0,style="italic",color=INK)
    ax.text(W+0.012,y,f"{lab} · {r['length']} aa",ha="left",va="center",fontsize=4.8,color=INK2)
ax.set_xlim(-0.17,1.33); ax.set_ylim(-0.62,len(ex)-0.28)
hs=[plt.Rectangle((0,0),1,1,facecolor=DC[k]) for k in ["CC","NB","WHD","LRR","ID"]]
ax.legend(hs,["CC (Rx_N)","NB-ARC","WHD / ARC2","LRR","integrated domain"],
          loc="lower center",bbox_to_anchor=(0.45,-0.085),ncol=5,fontsize=5.2)
panel(ax,"d",dx=-0.005,dy=1.035)
fig.savefig(f"{ROOT}/figures/Figure2.pdf"); fig.savefig(f"{ROOT}/figures/Figure2.png")
print("Figure 2 saved; ID-NLRs:",len(idrows),"ID families:",len(idc))
print("finger millet multi-ID NLRs:",len(fm))
