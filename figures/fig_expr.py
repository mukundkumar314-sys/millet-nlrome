import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import json, numpy as np, matplotlib.pyplot as plt, collections, math, statistics, glob
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup()
ORANGE="#D55E00"; BLUE="#0072B2"; GREY="#9a9a9a"; INK="#222222"; LIGHT="#d9d9d9"
R=json.load(open(f"{ROOT}/rna/induction.json"))
G=json.load(open(f"{ROOT}/rna/global_sanity.json"))
design={x['run']:x for x in json.load(open(f"{ROOT}/rna/design.json"))}
t2g=json.load(open(f"{ROOT}/rna/t2g.json")); gene_of={}
for g,ts in t2g.items():
    for t in ts: gene_of[t]=g
TPM=collections.defaultdict(lambda: collections.defaultdict(float))
for f in glob.glob(f"{ROOT}/rna/quant/*/quant.sf"):
    r=f.split("/")[-2]
    for i,L in enumerate(open(f)):
        if i==0: continue
        x=L.split("\t"); g=gene_of.get(x[0])
        if g: TPM[g][r]+=float(x[3])
pairs=[]
for geno in ['ACC','BK_FM','GBK','KAT_FM']:
    for t in [18,36,72]:
        c=[r for r in design if design[r]['geno']==geno and design[r]['time']==t and design[r]['treat']=='ctrl'][0]
        i=[r for r in design if design[r]['geno']==geno and design[r]['time']==t and design[r]['treat']=='inoc'][0]
        pairs.append((geno,t,c,i))

fig=plt.figure(figsize=(W2,W2*0.47))
gs=fig.add_gridspec(2,3,width_ratios=[1.25,0.85,0.95],height_ratios=[1,1],
                    hspace=0.70,wspace=0.48,left=0.09,right=0.985,top=0.89,bottom=0.17)

# a : candidate dot plot
ax=fig.add_subplot(gs[:,0])
CATC={"blastR candidate":BLUE,"head-to-head pair":ORANGE,"Clade A HMA target":"#009E73"}
R2=[r for r in R if r['total_reads_inoc']+r['total_reads_ctrl']>0]
R2.sort(key=lambda r:(r['cat'],r['median_log2FC']))
y=np.arange(len(R2))
for i,r in enumerate(R2):
    sig = r['q']<0.05
    ax.plot([0,r['median_log2FC']],[i,i],color=LIGHT,lw=0.8,zorder=1)
    ax.plot(r['median_log2FC'],i,'o',ms=6 if sig else 4,
            color=CATC[r['cat']],markeredgecolor=INK if sig else "white",
            markeredgewidth=0.9 if sig else 0.3,zorder=3)
    if sig: ax.text(r['median_log2FC']+0.08,i,f"q = {r['q']:.3f}",fontsize=6,va="center",fontweight="bold")
ax.axvline(0,color=INK,lw=0.6)
ax.set_yticks(y); ax.set_yticklabels([r['gene'] for r in R2],fontsize=5.6)
ax.set_xlabel("median log$_2$ fold change, inoculated vs mock\n(12 paired comparisons)",fontsize=6.5)
ax.set_xlim(-0.6,1.25)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=CATC[k],label=k) for k in CATC],frameon=False,fontsize=5.9,
          loc="lower right",bbox_to_anchor=(1.02,-0.015))
ax.spines[['top','right']].set_visible(False)
ax.set_title("a",loc="left",fontweight="bold",fontsize=9,x=-0.42,y=1.01)
ax.text(0,1.012,"24 pre-specified candidate genes; filled outline marks q < 0.05",
        transform=ax.transAxes,fontsize=6,color="#666666",va="bottom")

# b : the significant gene, per pair
ax=fig.add_subplot(gs[0,1])
g="QOZ80_6BG0485990"
lf=[math.log2((TPM[g][i]+1)/(TPM[g][c]+1)) for _,_,c,i in pairs]
cols=[{'ACC':BLUE,'BK_FM':ORANGE,'GBK':'#009E73','KAT_FM':'#CC79A7'}[p[0]] for p in pairs]
ax.bar(range(12),lf,color=cols,width=0.74)
ax.axhline(0,color=INK,lw=0.6)
ax.set_xticks(range(12)); ax.set_xticklabels([f"{p[1]}" for p in pairs],fontsize=5.5)
ax.set_xlabel("hours post inoculation\n(4 genotypes × 3 timepoints)",fontsize=6.3)
ax.set_ylabel("log$_2$ fold change",fontsize=6.4)
ax.set_title("b",loc="left",fontweight="bold",fontsize=9,x=-0.30,y=1.02)
ax.text(0,1.20,f"{g}: induced in 11 of 12 comparisons",transform=ax.transAxes,fontsize=6.2,color=ORANGE,va="bottom",fontweight="bold")
ax.legend(handles=[Patch(color=v,label=k) for k,v in {'ACC':BLUE,'BK_FM':ORANGE,'GBK':'#009E73','KAT_FM':'#CC79A7'}.items()],
          frameon=False,fontsize=5.6,ncol=4,loc="upper center",bbox_to_anchor=(0.5,1.20),handlelength=0.9,columnspacing=0.8)
ax.spines[['top','right']].set_visible(False)

# c : the head-to-head pair
ax=fig.add_subplot(gs[1,1])
g1,g2="QOZ80_6BG0485970","QOZ80_6BG0485990"
cm=[statistics.mean([TPM[g][c] for _,_,c,_ in pairs]) for g in (g1,g2)]
im=[statistics.mean([TPM[g][i] for _,_,_,i in pairs]) for g in (g1,g2)]
x=np.arange(2); w=0.34
ax.bar(x-w/2,cm,w,color=LIGHT,label="mock")
ax.bar(x+w/2,im,w,color=ORANGE,label="inoculated")
ax.set_xticks(x); ax.set_xticklabels(["...0485970\nPi5/Pii orthologue","...0485990\npartner"],fontsize=6.0)
ax.set_ylabel("mean TPM",fontsize=6.4)
ax.legend(frameon=False,fontsize=6,loc="upper left")
ax.set_title("c",loc="left",fontweight="bold",fontsize=9,x=-0.30,y=1.02)
ax.text(0,1.03,"the two members of one head-to-head pair\n(QOZ80_6BG0485970 and QOZ80_6BG0485990)",transform=ax.transAxes,fontsize=5.9,color="#666666",va="bottom")
ax.spines[['top','right']].set_visible(False)

# d : global sanity
ax=fig.add_subplot(gs[:,2])
ax.bar([0,1],[G['cons_up'],G['cons_dn']],color=[ORANGE,BLUE],width=0.6)
for xx,v in zip([0,1],[G['cons_up'],G['cons_dn']]):
    ax.text(xx,v+14,f"{v}",ha="center",fontsize=7,fontweight="bold")
ax.set_xticks([0,1]); ax.set_xticklabels(["induced in\n≥10 of 12","repressed in\n≥10 of 12"],fontsize=6.2)
ax.set_ylabel(f"genes (of {G['expressed']:,} expressed)",fontsize=6.4)
ax.set_ylim(0,G['cons_up']*1.25)
ax.set_title("d",loc="left",fontweight="bold",fontsize=9,x=-0.30,y=1.02)
ax.text(0,1.03,"inoculation produced a clear transcriptional\nresponse, confirming the comparison is informative",
        transform=ax.transAxes,fontsize=6,color="#666666",va="bottom")
ax.spines[['top','right']].set_visible(False)
ax.text(0.5,-0.22,"median log$_2$ fold change across all expressed genes = 0.00",
        transform=ax.transAxes,fontsize=5.9,color="#666666",ha="center",va="top")

plt.savefig(f"{ROOT}/figures/Figure_expression.png",dpi=600,bbox_inches="tight")
plt.savefig(f"{ROOT}/figures/Figure_expression.pdf",bbox_inches="tight")
print("saved")
