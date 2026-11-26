import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, json
from matplotlib.lines import Line2D
from Bio import Phylo
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]
ORD=ORD11
key={f'{r["species"]}|{r["gene"]}':r for r in rows}
pairs=json.load(open(f"{ROOT}/nlr/h2h_pairs.json"))
inpair=set()
for p in pairs:
    inpair.add(f'{p["species"]}|{p["g1"]}'); inpair.add(f'{p["species"]}|{p["g2"]}')
cand={c["fm_gene"] for c in json.load(open(f"{ROOT}/nlr/fm_blastR_candidates.json"))}

T=Phylo.read(f"{ROOT}/nlr/nbarc_tree11.nwk","newick"); T.ladderize()
tips=T.get_terminals(); N=len(tips)
ang={}; GAP=np.deg2rad(8)
for i,t in enumerate(tips): ang[t]=GAP/2+(2*np.pi-GAP)*i/N
depth=T.depths(unit_branch_lengths=False)
_tipd=sorted(depth[t] for t in tips)
SC=_tipd[int(0.80*len(_tipd))]          # 80th percentile of tip depth -> outer edge
def rad(c): return min(depth[c]/SC,1.0)
for c in T.get_nonterminals(order="postorder"):
    ang[c]=np.mean([ang[x] for x in c.clades])

fig=plt.figure(figsize=(W2,W2*0.56))
gs=fig.add_gridspec(1,13,wspace=1.2,left=0.01,right=0.985,top=0.97,bottom=0.03)

ax=fig.add_subplot(gs[0,0:8],projection="polar")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
ax.axis("off")
segs=[]
for c in T.get_nonterminals():
    r0=rad(c)
    a=[ang[x] for x in c.clades]
    th=np.linspace(min(a),max(a),max(int(abs(max(a)-min(a))*120),2))
    ax.plot(th,[r0]*len(th),color="#8a8a8a",lw=0.22,zorder=2,solid_capstyle="butt")
    for x in c.clades:
        ax.plot([ang[x],ang[x]],[r0,rad(x)],color="#8a8a8a",lw=0.22,zorder=2,solid_capstyle="butt")
R0=1.030
for t in tips:
    sp=t.name.split("|")[0]
    ax.plot([ang[t]],[rad(t)],marker="o",ms=0.9,color=COL.get(sp,"#999999"),zorder=4,
            markeredgewidth=0)
R1=1.085
for t in tips:
    r=key.get(t.name)
    if r and r["n_id"]>0:
        ax.plot([ang[t],ang[t]],[R1,R1+0.045],color="#D55E00",lw=0.45,zorder=5)
R2=1.145
for t in tips:
    if t.name in inpair:
        ax.plot([ang[t],ang[t]],[R2,R2+0.045],color="#332288",lw=0.45,zorder=5)
ax.set_ylim(0,1.27)
hs=[Line2D([],[],marker="o",ls="",ms=2.6,color=COL[s],label=LAB[s]) for s in ORD]
hs+=[Line2D([],[],color="#D55E00",lw=1.1,label="integrated domain"),
     Line2D([],[],color="#332288",lw=1.1,label="head-to-head pair")]
ax.legend(handles=hs,loc="center",bbox_to_anchor=(0.5,0.5),ncol=2,fontsize=4.6,
          frameon=False,handletextpad=0.4,labelspacing=0.28,columnspacing=0.8)
ax.text(np.deg2rad(2),1.31,"a",fontsize=8,fontweight="bold",color=INK,ha="left",va="top")

# ---- b: Pik clade ----
ax=fig.add_subplot(gs[0,8:13])
P=Phylo.read(f"{ROOT}/nlr/pik_clade.treefile","newick"); P.ladderize()
pt=P.get_terminals(); n=len(pt)
ypos={t:i for i,t in enumerate(pt)}
pd=P.depths(unit_branch_lengths=False); pmx=max(pd.values())
for c in P.get_nonterminals(order="postorder"):
    ypos[c]=np.mean([ypos[x] for x in c.clades])
def px(c): return pd[c]/pmx
for c in P.get_nonterminals():
    x0=px(c); ys=[ypos[x] for x in c.clades]
    ax.plot([x0,x0],[min(ys),max(ys)],color="#6a6a6a",lw=0.6,zorder=2)
    for x in c.clades:
        ax.plot([x0,px(x)],[ypos[x],ypos[x]],color="#6a6a6a",lw=0.6,zorder=2)
    if c.confidence is not None and c is not P.root:
        try: v=float(c.confidence)
        except: v=None
        if v is not None and v>=70:
            ax.text(x0-0.012,ypos[c]+0.17,f"{int(v)}",fontsize=4.0,color=MUTED,ha="right")
for t in pt:
    role,sp,gene=t.name.split("__")
    c="#D55E00" if role=="SENSOR" else "#0072B2"
    ax.plot([px(t)],[ypos[t]],marker="o",ms=2.4,color=c,zorder=4,markeredgewidth=0)
    lab=f"{LAB[sp]}  {gene}"
    ax.text(px(t)+0.022,ypos[t],lab,fontsize=4.4,va="center",color=INK)
    if gene in cand:
        ax.plot([1.52],[ypos[t]],marker="*",ms=4.0,color="#332288",zorder=5,markeredgewidth=0)
ax.set_xlim(-0.03,1.95); ax.set_ylim(-0.8,n-0.2); ax.invert_yaxis(); ax.axis("off")
hs=[Line2D([],[],marker="o",ls="",ms=2.6,color="#D55E00",label="$Pik$-1 sensor clade"),
    Line2D([],[],marker="o",ls="",ms=2.6,color="#0072B2",label="$Pik$-2 helper orthogroup"),
    Line2D([],[],marker="*",ls="",ms=3.8,color="#332288",label="blast R-gene candidate")]
ax.legend(handles=hs,loc="lower center",bbox_to_anchor=(0.45,-0.045),ncol=1,fontsize=4.6,frameon=False)
ax.text(0.45,1.0,"$Pik$ sensor and helper clades",transform=ax.transAxes,ha="center",
        fontsize=5.6,color=INK)
ax.text(-0.02,1.03,"b",transform=ax.transAxes,fontsize=8,fontweight="bold",color=INK)

fig.savefig(f"{ROOT}/figures/Figure8.pdf"); fig.savefig(f"{ROOT}/figures/Figure8.png")
print("Figure 8 saved; tips:",N,"pik clade tips:",n)
