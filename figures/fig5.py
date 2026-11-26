import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
pik=json.load(open(f"{ROOT}/refnlr/pik_data.json"))
hma_sa=json.load(open(f"{ROOT}/refnlr/hma_standalone.json"))
ORD=["Bdistachyon","Osativa","Ecoracana","Dexilis","Ecrusgalli","Phallii","Sviridis","Sitalica","Sbicolor","Zmays"]
HMA_C="#D55E00"
DC={"CC":"#999933","NB":"#0072B2","WHD":"#56B4E9","LRR":"#009E73","HMA":HMA_C}
def dcat(f):
    if f=="Rx_N": return "CC"
    if f=="NB-ARC": return "NB"
    if f in ("WHD_DRP","WHD_plant_disease"): return "WHD"
    if f.startswith("LRR"): return "LRR"
    return None

fig=plt.figure(figsize=(W2,W2*0.50))
gs=fig.add_gridspec(2,12,height_ratios=[0.78,1.0],hspace=0.58,wspace=3.0,
                    left=0.075,right=0.985,top=0.90,bottom=0.105)

# ---- a: rice Pik locus, head-to-head ----
ax=fig.add_subplot(gs[0,0:6]); ax.axis("off")
p1=pik["Os11t0688832-01"]; p2=pik["Os11t0689100-01"]
lo=min(p1["start"],p2["start"]); hi=max(p1["end"],p2["end"]); span=hi-lo
def X(v): return (v-lo)/span
ax.plot([0,1],[0.42,0.42],color="#c8c8c8",lw=1.0,zorder=1)
for p,nm,c in [(p1,"$Pik$-1  (sensor)",HMA_C),(p2,"$Pik$-2  (helper)","#0072B2")]:
    a,b=X(p["start"]),X(p["end"])
    ax.add_patch(plt.Rectangle((a,0.30),b-a,0.24,facecolor=c,alpha=0.85,lw=0,zorder=3))
    ar = -1 if p["strand"]=="-" else 1
    xs = b if ar==1 else a
    ax.annotate("",xy=(xs+0.045*ar,0.42),xytext=(xs,0.42),zorder=4,
                arrowprops=dict(arrowstyle="-|>",color=c,lw=1.0,mutation_scale=6))
    ax.text((a+b)/2,0.60,nm,ha="center",fontsize=5.4,color=INK)
    ax.text((a+b)/2,0.22,f"{p['strand']} strand",ha="center",fontsize=4.6,color=MUTED)
gap=p2["start"]-p1["end"]
ax.annotate("",xy=(X(p2["start"]),0.42),xytext=(X(p1["end"]),0.42),
            arrowprops=dict(arrowstyle="<->",color=INK2,lw=0.6,mutation_scale=5),zorder=5)
ax.text((X(p1["end"])+X(p2["start"]))/2,0.49,f"{gap:,} bp",ha="center",fontsize=4.8,color=INK2)
ax.text(0.5,0.03,f"O. sativa chr 11 : {lo:,}–{hi:,}",ha="center",fontsize=4.9,color=MUTED,style="italic")
ax.text(0.5,0.84,"Head-to-head $Pik$ sensor–helper pair",ha="center",fontsize=5.6,color=INK)
ax.set_xlim(-0.06,1.06); ax.set_ylim(0,1)
panel(ax,"a",dx=-0.045,dy=1.06)

# ---- b: HMA detection benchmark ----
ax=fig.add_subplot(gs[0,6:12]); ygrid(ax,axis="x")
labs=["Pfam HMA (PF00403)\ngathering threshold","Pfam HMA (PF00403)\nE ≤ 10",
      "Pikm1/RGA5 HMA seed\n(this study)"]
vals=[0,1,1]; tp=[0,0,1]
y=np.arange(3)
ax.barh(y,vals,height=0.55,color=["#d9d9d9","#d9d9d9",HMA_C],zorder=3)
for i in range(3):
    txt = "rice $Pik$-1 detected" if tp[i] else ("0 hits" if vals[i]==0 else "1 hit, not $Pik$-1")
    ax.text(vals[i]+0.05,i,txt,va="center",fontsize=4.9,color=INK if tp[i] else MUTED)
ax.set_yticks(y); ax.set_yticklabels(labs,fontsize=4.9)
ax.invert_yaxis(); ax.set_xlabel("HMA-integrated NLRs detected (n)")
ax.set_xlim(0,2.3); ax.set_xticks([0,1,2])
ax.text(0.5,1.09,"Benchmark on 3,012 NLRs; positive control = rice $Pik$-1 (bit score 13.4)",
        transform=ax.transAxes,ha="center",fontsize=4.7,color=MUTED)
panel(ax,"b",dx=-0.37,dy=1.14)

# ---- c: architectures, rice Pik-1 vs millet orthologs ----
ax=fig.add_subplot(gs[1,0:7]); ax.axis("off")
order=[("Os11t0688832-01","O. sativa $Pik$-1",True),
       ("TKV99356","S. viridis ortholog",False),
       ("TKW02560","S. viridis ortholog",False),
       ("scaffold389.6.cds","E. crus-galli ortholog",False),
       ("EES09001","S. bicolor ortholog",False)]
maxlen=max(pik[k]["length"] for k,_,_ in order); H=0.42; y=len(order)
for k,lab,ishma in order:
    r=pik[k]; y-=1; W=r["length"]/maxlen
    ax.plot([0,W],[y,y],color="#c8c8c8",lw=1.0,zorder=1)
    drawn=[]
    for fam,a0,b0 in r["domains"]:
        c=dcat(fam)
        if c is None: continue
        a,b=a0/maxlen,b0/maxlen
        if any(not(b<=p or a>=q) for p,q in drawn): continue
        drawn.append((a,b))
        ax.add_patch(plt.Rectangle((a,y-H/2),max(b-a,0.004),H,facecolor=DC[c],
                                   edgecolor="white",lw=0.4,zorder=3))
    if ishma:
        a,b=237/maxlen,287/maxlen
        ax.add_patch(plt.Rectangle((a,y-H/2),b-a,H,facecolor=HMA_C,edgecolor="white",lw=0.4,zorder=4))
        ax.text((a+b)/2,y+H/2+0.13,"HMA",ha="center",fontsize=4.8,color=HMA_C,fontweight="bold")
    else:
        ax.text(262/maxlen,y,"×",ha="center",va="center",fontsize=6,color="#b0b0b0",zorder=5)
    ax.text(-0.012,y,lab,ha="right",va="center",fontsize=5.0,color=INK)
    ax.text(W+0.012,y,f"{r['length']} aa",ha="left",va="center",fontsize=4.8,color=INK2)
ax.axvline(262/maxlen,color="#d9d9d9",lw=0.5,ls=(0,(2,2)),zorder=0)
ax.set_xlim(-0.30,1.18); ax.set_ylim(-0.75,len(order)-0.20)
hs=[plt.Rectangle((0,0),1,1,facecolor=DC[k]) for k in ["CC","NB","WHD","LRR","HMA"]]
ax.legend(hs,["CC (Rx_N)","NB-ARC","WHD / ARC2","LRR","integrated HMA"],
          loc="lower center",bbox_to_anchor=(0.42,-0.14),ncol=5,fontsize=5.0)
ax.text(0.42,1.02,"× marks the position where the HMA is integrated in rice $Pik$-1",
        transform=ax.transAxes,ha="center",fontsize=4.7,color=MUTED)
panel(ax,"c",dx=-0.03,dy=1.10)

# ---- d: standalone HMA proteins vs HMA-integrated NLRs ----
ax=fig.add_subplot(gs[1,7:12]); ygrid(ax)
x=np.arange(len(ORD))
sa=[hma_sa.get(s,0) for s in ORD]
ax.bar(x,sa,width=0.68,color="#b8c7d9",zorder=3,label="standalone HMA proteins")
inte=[1 if s=="Osativa" else 0 for s in ORD]
ax.bar(x,inte,width=0.68,bottom=sa,color=HMA_C,zorder=4,label="HMA-integrated NLRs")
for i,s in enumerate(ORD):
    if inte[i]: ax.annotate("1",(i,sa[i]+1),xytext=(0,9),textcoords="offset points",
                            ha="center",fontsize=5,color=HMA_C,fontweight="bold",
                            arrowprops=dict(arrowstyle="-",color=HMA_C,lw=0.5))
ax.set_xticks(x); ax.set_xticklabels([LAB[s] for s in ORD],rotation=40,ha="right",style="italic")
ax.set_ylabel("Proteins (n)"); ax.set_ylim(0,max(sa)*1.30)
ax.legend(loc="upper center",bbox_to_anchor=(0.5,1.14),ncol=1,fontsize=5.0)
panel(ax,"d",dx=-0.22,dy=1.17)

fig.savefig(f"{ROOT}/figures/Figure5.pdf"); fig.savefig(f"{ROOT}/figures/Figure5.png")
print("Figure 5 saved")
print("Pik-1/Pik-2 intergenic gap:",gap,"bp; strands",p1["strand"],p2["strand"])
