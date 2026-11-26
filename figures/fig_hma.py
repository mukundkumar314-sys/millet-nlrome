import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import json, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from Bio import Phylo
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup()
ORANGE="#D55E00"; BLUE="#0072B2"; GREY="#9a9a9a"; INK="#222222"

rbh=json.load(open(f"{ROOT}/hma/cladeA_rbh.json"))
A1="Os_OsHIPP19_Os04g0469000"; A2="Os_OsHIPP20_Os04g0464100"
C4A="Ec_KAK3144512.1"; C4B="Ec_KAK3142336.1"
NAME={A1:"OsHIPP19",A2:"OsHIPP20",C4A:"QOZ80_4AG0314090",C4B:"QOZ80_4BG0345250"}
clade=set(json.load(open(f"{ROOT}/hma/cladeA.json"))["clade_tips"])

fig=plt.figure(figsize=(W2,W2*0.56))
gs=fig.add_gridspec(2,3,width_ratios=[1.05,0.95,0.95],height_ratios=[1,1],
                    hspace=0.55,wspace=0.52,left=0.045,right=0.985,top=0.93,bottom=0.16)

# ---------------- a : whole tree ----------------
ax=fig.add_subplot(gs[:,0])
t=Phylo.read(f"{ROOT}/hma/shma.treefile","newick"); t.root_at_midpoint()
tips=t.get_terminals(); ypos={x.name:i for i,x in enumerate(tips)}
D={}
def walk(n,acc=0.0):
    acc=acc+(n.branch_length or 0.0); D[n]=acc
    for c in n.clades: walk(c,acc)
walk(t.root)
def ymid(n):
    return ypos[n.name] if n.is_terminal() else (lambda v:(min(v)+max(v))/2)([ymid(c) for c in n.clades])
for n in t.find_clades(order="postorder"):
    y=ymid(n); x=D[n]
    inA=n.is_terminal() and n.name in clade
    if n.branch_length:
        ax.plot([x-(n.branch_length or 0),x],[y,y],lw=1.6 if inA else 0.4,
                color=ORANGE if inA else GREY,solid_capstyle="butt",zorder=3 if inA else 1)
    if not n.is_terminal():
        ys=[ymid(c) for c in n.clades]
        ax.plot([x,x],[min(ys),max(ys)],lw=0.4,color=GREY,zorder=1)
ymin=min(ypos[n] for n in clade); ymax=max(ypos[n] for n in clade)
xmax=max(D.values())
ax.add_patch(Rectangle((-0.03,ymin-1.6),xmax*1.06,(ymax-ymin)+3.2,
            facecolor=ORANGE,alpha=0.10,edgecolor="none",zorder=0))
ax.annotate("Clade A",xy=(xmax*0.80,(ymin+ymax)/2),xytext=(xmax*0.52,(ymin+ymax)/2-16),
            fontsize=7,color=ORANGE,fontweight="bold",
            arrowprops=dict(arrowstyle="-|>",color=ORANGE,lw=1.1,shrinkA=0,shrinkB=2))
ax.text(0.0,1.005,"152 small HMA domains · 57 rice, 95 finger millet\nLG+G4, 1000 ultrafast bootstrap replicates",
        transform=ax.transAxes,fontsize=5.9,color="#666666",va="bottom")
ax.set_ylim(-3,len(tips)+1); ax.set_xlim(-0.05,xmax*1.12); ax.axis("off")
ax.set_title("a",loc="left",fontweight="bold",fontsize=9,x=-0.03,y=1.03)

# ---------------- b : clade A zoom ----------------
ax=fig.add_subplot(gs[0,1])
sub=t.common_ancestor([A1,A2])
st=sub.get_terminals(); sy={x.name:i for i,x in enumerate(st)}
base=D[sub]-(sub.branch_length or 0)
def ymid2(n): return sy[n.name] if n.is_terminal() else (lambda v:(min(v)+max(v))/2)([ymid2(c) for c in n.clades])
for n in sub.find_clades(order="postorder"):
    y=ymid2(n); x=D[n]-base
    if n.branch_length: ax.plot([x-(n.branch_length or 0),x],[y,y],lw=1.5,color=ORANGE,solid_capstyle="butt")
    if not n.is_terminal():
        ys=[ymid2(c) for c in n.clades]; ax.plot([x,x],[min(ys),max(ys)],lw=1.5,color=ORANGE)
        if n.confidence: ax.text(x-0.006,(min(ys)+max(ys))/2-0.16,str(int(n.confidence)),
                                 fontsize=6,color=ORANGE,ha="right",va="top",fontweight="bold")
for n in st:
    x=D[n]-base; y=sy[n.name]
    ax.plot(x,y,"o",ms=4,color=ORANGE,markeredgecolor="white",markeredgewidth=0.5,zorder=4)
    ax.text(x+0.012,y,NAME[n.name],fontsize=6.3,va="center",
            fontweight="bold" if n.name.startswith("Ec_") else "normal",color=INK)
ax.set_xlim(-0.015,(max(D[n] for n in st)-base)*2.6); ax.set_ylim(-0.7,len(st)-0.3); ax.axis("off")
ax.set_title("b",loc="left",fontweight="bold",fontsize=9,x=-0.06,y=1.02)
ax.text(0,1.0,"Clade A, enlarged",transform=ax.transAxes,fontsize=6.4,color=ORANGE,fontweight="bold",va="bottom")
ax.text(0,-0.04,"bold = finger millet; plain = rice",transform=ax.transAxes,fontsize=5.8,color="#666666",va="top")

# ---------------- c : identity ----------------
ax=fig.add_subplot(gs[1,1])
labs=["QOZ80_4AG0314090","QOZ80_4BG0345250"]
v19=[80.6,75.8]; v20=[85.0,80.0]
y=np.arange(2); h=0.33
ax.barh(y+h/2,v19,h,color=BLUE); ax.barh(y-h/2,v20,h,color=ORANGE)
for yy,v in zip(y+h/2,v19): ax.text(v+1.5,yy,f"{v:.1f}",va="center",fontsize=6)
for yy,v in zip(y-h/2,v20): ax.text(v+1.5,yy,f"{v:.1f}",va="center",fontsize=6)
ax.set_yticks(y); ax.set_yticklabels(labs,fontsize=6.2)
ax.set_xlabel("HMA-domain identity (%)",fontsize=6.6); ax.set_xlim(0,108)
ax.set_xticks([0,20,40,60,80,100])
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE,label="vs OsHIPP19"),Patch(color=ORANGE,label="vs OsHIPP20")],
          frameon=False,fontsize=6,loc="upper center",bbox_to_anchor=(0.5,-0.40),ncol=2,handlelength=1.1)
ax.spines[['top','right']].set_visible(False)
ax.set_title("c",loc="left",fontweight="bold",fontsize=9,x=-0.06,y=1.0)

# ---------------- d : across 11 genomes ----------------
ax=fig.add_subplot(gs[:,2])
order=ORD11
pid=[max(rbh[s][0][1],rbh[s][1][1]) for s in order]
col=[ORANGE if s=="Ecoracana" else (BLUE if s=="Osativa" else GREY) for s in order]
y=np.arange(len(order))[::-1]
ax.barh(y,pid,color=col,height=0.68)
for yy,p in zip(y,pid): ax.text(p+1.5,yy,f"{p:.0f}",va="center",fontsize=5.9)
ax.set_yticks(y); ax.set_yticklabels([LAB.get(s,s) for s in order],fontsize=6.2,style="italic")
ax.set_xlabel("identity to nearest validated\nAVR-Pik target (%)",fontsize=6.6)
ax.set_xlim(0,100); ax.set_xticks([0,25,50,75,100])
ax.text(0.0,1.012,"a Clade A orthologue occurs in every genome examined;\nHMA integration into an NLR occurs in rice alone",
        transform=ax.transAxes,fontsize=6,color="#666666",va="bottom")
ax.spines[['top','right']].set_visible(False)
ax.set_title("d",loc="left",fontweight="bold",fontsize=9,x=-0.22,y=1.03)

plt.savefig(f"{ROOT}/figures/Figure_HMA.png",dpi=600,bbox_inches="tight")
plt.savefig(f"{ROOT}/figures/Figure_HMA.pdf",bbox_inches="tight")
print("saved")
