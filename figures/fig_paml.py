import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import json, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup()
ORANGE="#D55E00"; BLUE="#0072B2"; GREY="#9a9a9a"; INK="#222222"; LIGHT="#d9d9d9"
res=json.load(open(f"{ROOT}/paml/codeml_results.json"))
sm =json.load(open(f"{ROOT}/paml/site_mapping.json"))
dom=json.load(open(f"{ROOT}/paml/og_domains.json"))
OGS=["OG0021","OG0025","OG0083","OG0213"]
LBL={"OG0021":"OG0021\nPi50, Piz-t","OG0025":"OG0025\nPik, Pikm, Pi1",
     "OG0083":"OG0083\nPi2, Pi9","OG0213":"OG0213\nPi5, Pii"}
NSEQ={"OG0021":18,"OG0025":17,"OG0083":8,"OG0213":3}

fig=plt.figure(figsize=(W2,W2*0.42))
gs=fig.add_gridspec(2,2,width_ratios=[0.85,1.35],height_ratios=[1,1],
                    hspace=0.75,wspace=0.30,left=0.085,right=0.985,top=0.88,bottom=0.17)

# a : LRTs
ax=fig.add_subplot(gs[0,0])
x=np.arange(len(OGS)); w=0.36
v1=[res[o]["lrt"]["M1a vs M2a"]["2dlnL"] for o in OGS]
v2=[res[o]["lrt"]["M7 vs M8"]["2dlnL"] for o in OGS]
b1=ax.bar(x-w/2,v1,w,color=BLUE,label="M1a vs M2a")
b2=ax.bar(x+w/2,v2,w,color=ORANGE,label="M7 vs M8")
for xi,o in zip(x,OGS):
    for off,key in [(-w/2,"M1a vs M2a"),(w/2,"M7 vs M8")]:
        p=res[o]["lrt"][key]["p"]; d=res[o]["lrt"][key]["2dlnL"]
        ax.text(xi+off,d+3,("n.s." if p>=0.05 else f"{p:.0e}".replace("e-0","e-")),
                ha="center",fontsize=5.3,rotation=90,va="bottom")
ax.axhline(5.99,color=GREY,ls="--",lw=0.7)
ax.text(-0.45,9,"5% threshold",fontsize=5.4,color=GREY,ha="left")
ax.set_xticks(x); ax.set_xticklabels([LBL[o] for o in OGS],fontsize=5.8)
ax.set_ylabel("2$\\Delta$lnL",fontsize=6.5); ax.set_ylim(0,165)
ax.legend(frameon=False,fontsize=5.8,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.30),handlelength=1.0)
ax.spines[['top','right']].set_visible(False)
ax.set_title("a",loc="left",fontweight="bold",fontsize=9,x=-0.26,y=1.28)

# b : LRR enrichment for OG0021
ax=fig.add_subplot(gs[1,0])
t=sm["OG0021"]["lrr_test"]
obs=[t["in_lrr"]/ (t["in_lrr"]+t["out_lrr"])*100, (t["out_lrr"])/(t["in_lrr"]+t["out_lrr"])*100]
exp=[t["lrr_aa"]/t["prot_aa"]*100, 100-t["lrr_aa"]/t["prot_aa"]*100]
xx=np.arange(2); w2=0.36
ax.bar(xx-w2/2,[exp[0],exp[1]],w2,color=LIGHT,label="expected from length")
ax.bar(xx+w2/2,[obs[0],obs[1]],w2,color=ORANGE,label="observed")
for i,(e,o) in enumerate(zip(exp,obs)):
    ax.text(i-w2/2,e+1.5,f"{e:.0f}%",ha="center",fontsize=5.8)
    ax.text(i+w2/2,o+1.5,f"{o:.0f}%",ha="center",fontsize=5.8,fontweight="bold")
ax.set_xticks(xx); ax.set_xticklabels(["inside LRR","outside LRR"],fontsize=6.2)
ax.set_ylabel("share (%)",fontsize=6.5); ax.set_ylim(0,100)
ax.legend(frameon=False,fontsize=5.8,loc="upper center",ncol=2,bbox_to_anchor=(0.5,1.22),handlelength=1.0)
ax.text(0,1.04,f"OG0021 selected sites: OR = {t['OR']:.2f}, P = {t['p']:.4f}",
        transform=ax.transAxes,fontsize=6,color=ORANGE,fontweight="bold",va="bottom")
ax.spines[['top','right']].set_visible(False)
ax.set_title("b",loc="left",fontweight="bold",fontsize=9,x=-0.26,y=1.20)

# c : site map
ax=fig.add_subplot(gs[:,1])
DC={"Rx_N":"#999933","NB-ARC":BLUE,"WHD_DRP":"#56B4E9","LRR":"#009E73","AAA":"#b0b0b0"}
rows=[o for o in OGS]
for yi,og in enumerate(rows):
    y=len(rows)-1-yi
    rec=[m for m in dom[og]["members"] if m["tag"]=="s00_Ecorac"][0]
    Lp=rec["len"]
    ax.add_patch(Rectangle((0,y-0.13),Lp,0.26,facecolor="#eeeeee",edgecolor="none"))
    for d in rec["domains"]:
        fam="LRR" if d["fam"].startswith("LRR") else d["fam"]
        ax.add_patch(Rectangle((d["s"],y-0.19),d["e"]-d["s"],0.38,
                     facecolor=DC.get(fam,GREY),edgecolor="none",alpha=0.92))
    sig=[s for s in sm[og]["sites"] if s["star"] and s["ref_pos"]]
    for s in sig:
        ax.plot(s["ref_pos"],y+0.33,marker="v",ms=3.6 if s["star"]=="*" else 4.6,
                color=ORANGE if s["star"]=="**" else "#E69F00",
                markeredgecolor=INK,markeredgewidth=0.25,zorder=5)
    note=f"{rec['gene']}  ({NSEQ[og]} seqs)"
    if og=="OG0213": note+="  — too few sequences, unreliable"
    ax.text(-30,y,LBL[og].replace("\n"," · "),fontsize=5.9,ha="right",va="center")
    ax.text(Lp+25,y,note,fontsize=5.6,va="center",color="#666666" if og!="OG0213" else "#b03030")
ax.set_xlim(-330,1250); ax.set_ylim(-0.75,len(rows)-0.25)
ax.set_yticks([]); ax.set_xlabel("position in the finger millet member (aa)",fontsize=6.5)
for s in ["top","right","left"]: ax.spines[s].set_visible(False)
h=[Patch(color=DC[k],label=k) for k in ["Rx_N","NB-ARC","WHD_DRP","LRR"]]
h+=[plt.Line2D([],[],marker="v",ls="",color="#E69F00",mec=INK,mew=0.25,ms=3.6,label="site P > 0.95"),
    plt.Line2D([],[],marker="v",ls="",color=ORANGE,mec=INK,mew=0.25,ms=4.6,label="site P > 0.99")]
ax.legend(handles=h,frameon=False,fontsize=5.7,ncol=6,loc="upper center",
          bbox_to_anchor=(0.5,1.13),handlelength=1.0,columnspacing=1.0)
ax.set_title("c",loc="left",fontweight="bold",fontsize=9,x=-0.02,y=1.11)
plt.savefig(f"{ROOT}/figures/Figure_selection.png",dpi=600,bbox_inches="tight")
plt.savefig(f"{ROOT}/figures/Figure_selection.pdf",bbox_inches="tight")
print("saved")
