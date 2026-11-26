import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import json, numpy as np, matplotlib.pyplot as plt, collections
from matplotlib.patches import Patch
from Bio import Phylo
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup()
ORANGE="#D55E00"; BLUE="#0072B2"; GREEN="#009E73"; GREY="#9a9a9a"; INK="#222222"
COLL={"Oryza":BLUE,"Eleusine":ORANGE,"Triticum":GREEN}
lin=json.load(open(f"{ROOT}/mag/lineage.json"))
cnt=json.load(open(f"{ROOT}/mag/max_counts.json"))
tot=collections.Counter(lin.values())

fig=plt.figure(figsize=(W2,W2*0.52))
gs=fig.add_gridspec(2,3,width_ratios=[1.0,1.05,1.0],height_ratios=[1,1],
                    hspace=0.62,wspace=0.45,left=0.05,right=0.985,top=0.90,bottom=0.15)

# ---- a : tree ----
ax=fig.add_subplot(gs[:,0])
t=Phylo.read(f"{ROOT}/mag/hid.treefile","newick"); t.root_at_midpoint()
t.ladderize()
tips=t.get_terminals(); ypos={x.name:i for i,x in enumerate(tips)}
D={}
def _dep(n):
    return 0 if n.is_terminal() else 1+max(_dep(c) for c in n.clades)
MAXD=_dep(t.root)
def walk(n,lev=0):
    D[n]= MAXD if n.is_terminal() else lev
    for c in n.clades: walk(c,lev+1)
walk(t.root)
def ymid(n): return ypos[n.name] if n.is_terminal() else (lambda v:(min(v)+max(v))/2)([ymid(c) for c in n.clades])
for n in t.find_clades(order="postorder"):
    y=ymid(n); x=D[n]
    par=[p for p in t.find_clades() if n in p.clades]
    if par:
        xp=D[par[0]]
        c=COLL[{"Ory":"Oryza","Ele":"Eleusine","Tri":"Triticum"}[n.name[:3]]] if n.is_terminal() else GREY
        ax.plot([xp,x],[y,y],lw=1.0 if n.is_terminal() else 0.45,color=c,solid_capstyle="butt")
    if not n.is_terminal():
        ys=[ymid(c_) for c_ in n.clades]; ax.plot([x,x],[min(ys),max(ys)],lw=0.45,color=GREY)
for n in tips:
    ax.plot(D[n],ypos[n.name],"o",ms=3.0,color=COLL[{"Ory":"Oryza","Ele":"Eleusine","Tri":"Triticum"}[n.name[:3]]],
            markeredgecolor="white",markeredgewidth=0.25,zorder=4)
ax.set_xlim(-0.4,MAXD+1.6); ax.set_ylim(-3.5,len(tips)+1); ax.axis("off")
ax.legend(handles=[Patch(color=COLL[k],label=f"{k}-infecting") for k in ["Oryza","Eleusine","Triticum"]],
          frameon=False,fontsize=6,loc="lower left",bbox_to_anchor=(0.0,-0.055),ncol=1)
ax.text(0,1.005,"65 AVR-Pik-like HMA-interaction-domain\nproteins; WAG+F, cladogram",
        transform=ax.transAxes,fontsize=5.8,color="#666666",va="bottom")
ax.set_title("a",loc="left",fontweight="bold",fontsize=9,x=-0.02,y=1.04)

# ---- b : family presence ----
ax=fig.add_subplot(gs[0,1])
fams=["AVR-Pik_HID","AvrPi9","AvrPiz-t","MAX_effector","Avr_Pia"]
lins=["Oryza","Eleusine","Triticum"]
M=np.zeros((len(fams),3))
for i,f in enumerate(fams):
    for j,L in enumerate(lins):
        accs=[a for a in cnt.get(f,{}) if lin[a]==L]
        M[i,j]=100*len(accs)/tot[L]
im=ax.imshow(M,cmap="Oranges",vmin=0,vmax=100,aspect="auto")
for i in range(len(fams)):
    for j in range(3):
        accs=[a for a in cnt.get(fams[i],{}) if lin[a]==lins[j]]
        ax.text(j,i,f"{len(accs)}/{tot[lins[j]]}",ha="center",va="center",fontsize=6,
                color="white" if M[i,j]>55 else INK)
ax.set_xticks(range(3)); ax.set_xticklabels([f"{l}\n(n={tot[l]})" for l in lins],fontsize=6)
ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams,fontsize=6.2)
ax.set_title("b",loc="left",fontweight="bold",fontsize=9,x=-0.42,y=1.02)
ax.text(0,1.03,"isolates carrying at least one member",transform=ax.transAxes,fontsize=6,color="#666666",va="bottom")
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)

# ---- c : copy number ----
ax=fig.add_subplot(gs[1,1])
dat=[[cnt["AVR-Pik_HID"].get(a,0) for a in lin if lin[a]==L and a in cnt["AVR-Pik_HID"]] for L in lins]
bp=ax.boxplot(dat,patch_artist=True,widths=0.55,medianprops=dict(color=INK,lw=1.2),
              flierprops=dict(marker='o',ms=2.5,mfc=GREY,mec='none'))
for p,L in zip(bp['boxes'],lins): p.set_facecolor(COLL[L]); p.set_alpha(0.75); p.set_edgecolor('none')
for i,d in enumerate(dat,1):
    ax.plot(np.random.normal(i,0.055,len(d)),d,'o',ms=2.2,color=INK,alpha=0.45,zorder=3)
ax.set_xticklabels(lins,fontsize=6.3)
ax.set_ylabel("AVR-Pik_HID copies\nper isolate",fontsize=6.4)
ax.spines[['top','right']].set_visible(False)
ax.set_title("c",loc="left",fontweight="bold",fontsize=9,x=-0.42,y=1.0)

# ---- d : annotation gap ----
ax=fig.add_subplot(gs[0,2])
cats=["present in the\ngenome assembly","present in the\nannotated proteome"]
vals=[33,1]
b=ax.barh([1,0],vals,color=[ORANGE,"#d9d9d9"],height=0.5)
for yy,v in zip([1,0],vals): ax.text(v+0.7,yy,f"{v}/33",va="center",fontsize=7,fontweight="bold")
ax.set_yticks([1,0]); ax.set_yticklabels(cats,fontsize=6.4)
ax.set_xlim(0,42); ax.set_ylim(-0.55,1.55); ax.set_xlabel("Eleusine-infecting isolates",fontsize=6.6)
ax.spines[['top','right']].set_visible(False)
ax.text(0,1.10,"the canonical AVR-Pik-family gene is present in every\n"
               "Eleusine-infecting genome, yet annotated in one",
        transform=ax.transAxes,fontsize=6.2,color="#666666",va="bottom")
ax.text(0.0,-0.34,"genomic detection by tblastn, 95.6–100% amino-acid identity\nover the full 113-residue protein",
        transform=ax.transAxes,fontsize=5.9,color="#666666",va="top")
ax.set_title("d",loc="left",fontweight="bold",fontsize=9,x=-0.30,y=1.12)

ax=fig.add_subplot(gs[1,2])
import collections as _c
best=_c.defaultdict(float)
for L in open(f"{ROOT}/mag/tblastn_all.tsv"):
    f=L.rstrip("\n").split("\t")
    if len(f)<9: continue
    if float(f[6])>best.get(f[0]+"_b",0): best[f[0]+"_b"]=float(f[6]); best[f[0]]=float(f[3])
el=[best[a] for a in best if not a.endswith("_b") and lin.get(a)=="Eleusine"]
ax.hist(el,bins=np.arange(94,101,0.5),color=ORANGE,edgecolor="white",linewidth=0.5)
ax.set_xlabel("amino-acid identity to canonical\nAVR-Pik family protein (%)",fontsize=6.4)
ax.set_ylabel("Eleusine-infecting\ngenomes",fontsize=6.4)
ax.spines[['top','right']].set_visible(False)
ax.set_title("e",loc="left",fontweight="bold",fontsize=9,x=-0.30,y=1.02)
ax.text(0,1.03,f"n = {len(el)} genomes, all above 95% identity",transform=ax.transAxes,fontsize=6,color="#666666",va="bottom")

plt.savefig(f"{ROOT}/figures/Figure_effector.png",dpi=600,bbox_inches="tight")
plt.savefig(f"{ROOT}/figures/Figure_effector.pdf",bbox_inches="tight")
print("saved")
