import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, re, json
from collections import Counter, defaultdict
from scipy import stats
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; cl=d["clusters"]
fm=[r for r in rows if r["species"]=="Ecoracana"]
CHL={}
for ln in open(f"{ROOT}/gff/Ecoracana.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': CHL[f[2]]=int(f[8])
PAIRS=[f"{i}" for i in range(1,10)]
CA="#D55E00"; CB="#0072B2"
g=cl["Ecoracana"]["g"]*1000

fig=plt.figure(figsize=(W2,W2*0.60))
gs=fig.add_gridspec(2,12,height_ratios=[1.18,0.82],hspace=0.46,wspace=3.0,
                    left=0.075,right=0.985,top=0.93,bottom=0.095)

# ---- a: ideogram ----
ax=fig.add_subplot(gs[0,0:12])
maxL=max(CHL.values()); BW=0.30
for i,p in enumerate(PAIRS):
    for k,(sub,c) in enumerate([("A",CA),("B",CB)]):
        ch=f"{p}{sub}"; L=CHL[ch]/1e6
        y=i+(-0.19 if k==0 else 0.19)
        ax.add_patch(plt.Rectangle((0,y-BW/2),L,BW,facecolor="#f0f0f0",
                                   edgecolor="#c0c0c0",lw=0.35,zorder=2))
        pos=sorted((r["start"]+r["end"])//2 for r in fm if r.get("chrom_lab")==ch)
        # clusters
        cur=[pos[0]] if pos else []
        for q in range(1,len(pos)):
            if pos[q]-pos[q-1]<=g: cur.append(pos[q])
            else:
                if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),
                    max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.22,lw=0,zorder=3))
                cur=[pos[q]]
        if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),
            max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.22,lw=0,zorder=3))
        for x0 in pos:
            ax.plot([x0/1e6,x0/1e6],[y-BW/2,y+BW/2],color=c,lw=0.45,zorder=4,alpha=0.92)
        ax.text(-1.5,y,ch,ha="right",va="center",fontsize=5.0,color=INK)
        ax.text(L+1.2,y,f"{len(pos)}",ha="left",va="center",fontsize=4.6,color=c)
ax.set_xlim(-7,maxL/1e6+7); ax.set_ylim(-0.65,len(PAIRS)-0.35)
ax.set_xlabel("Chromosome position (Mb)"); ax.set_yticks([]); ax.invert_yaxis()
for sp_ in ["left","right","top"]: ax.spines[sp_].set_visible(False)
hs=[plt.Rectangle((0,0),1,1,facecolor=CA),plt.Rectangle((0,0),1,1,facecolor=CB)]
ax.legend(hs,["A subgenome","B subgenome"],loc="lower center",
          bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=5.4)
ax.text(0.995,0.015,"shaded block = NLR cluster (≥3 genes)",transform=ax.transAxes,
        ha="right",fontsize=4.6,color=MUTED)
panel(ax,"a",dx=-0.052,dy=1.085)

# ---- b: homoeologous counts ----
ax=fig.add_subplot(gs[1,0:4]); ygrid(ax)
A=[sum(1 for r in fm if r.get("chrom_lab")==f"{p}A") for p in PAIRS]
B=[sum(1 for r in fm if r.get("chrom_lab")==f"{p}B") for p in PAIRS]
x=np.arange(len(PAIRS))
ax.bar(x-0.19,A,width=0.36,color=CA,label="A",zorder=3)
ax.bar(x+0.19,B,width=0.36,color=CB,label="B",zorder=3)
ax.set_xticks(x); ax.set_xticklabels(PAIRS); ax.set_xlabel("Homoeologous chromosome group")
ax.set_ylabel("NLR genes (n)"); ax.legend(fontsize=5.2,ncol=2,loc="upper left")
r,pv=stats.pearsonr(A,B)
ax.set_ylim(0,34)
ax.text(0.98,0.995,f"$r$ = {r:.2f}, $P$ = {pv:.1e}",transform=ax.transAxes,ha="right",va="top",
        fontsize=4.9,color=MUTED)
panel(ax,"b",dx=-0.18)

# ---- c: subgenome size vs NLR ----
ax=fig.add_subplot(gs[1,4:8]); ygrid(ax)
LA=sum(CHL[f"{p}A"] for p in PAIRS)/1e6; LB=sum(CHL[f"{p}B"] for p in PAIRS)/1e6
nA,nB=sum(A),sum(B)
vals=[[LA,LB],[nA,nB],[nA/LA*100,nB/LB*100]]
labs=["Assembly\nsize (Mb)","NLR genes\n(n)","NLR density\n(per 100 Mb)"]
xx=np.arange(3)
for k in range(3):
    mx=max(vals[k])
    ax.bar(xx[k]-0.19,vals[k][0]/mx*100,width=0.36,color=CA,zorder=3)
    ax.bar(xx[k]+0.19,vals[k][1]/mx*100,width=0.36,color=CB,zorder=3)
    ax.text(xx[k]-0.19,vals[k][0]/mx*100+1.6,f"{vals[k][0]:.0f}" if k<2 else f"{vals[k][0]:.1f}",
            ha="center",fontsize=4.7,color=INK)
    ax.text(xx[k]+0.19,vals[k][1]/mx*100+1.6,f"{vals[k][1]:.0f}" if k<2 else f"{vals[k][1]:.1f}",
            ha="center",fontsize=4.7,color=INK)
ax.set_xticks(xx); ax.set_xticklabels(labs,fontsize=5.0)
ax.set_ylabel("Relative to larger subgenome (%)"); ax.set_ylim(0,118)
chi=stats.chisquare([nA,nB],[ (nA+nB)*LA/(LA+LB), (nA+nB)*LB/(LA+LB) ])
ax.text(0.5,1.015,f"NLR number vs assembly size: $\\chi^2$ = {chi.statistic:.1f}, $P$ = {chi.pvalue:.3f}",
        transform=ax.transAxes,ha="center",fontsize=4.7,color=MUTED)
panel(ax,"c",dx=-0.27)

# ---- d: architecture A vs B ----
ax=fig.add_subplot(gs[1,8:12]); ygrid(ax)
ARCH=["N","CN","NL","CNL","ID-NLR"]
cA=Counter(archgroup(r["arch"]) for r in fm if re.fullmatch(r'\dA',str(r.get("chrom_lab"))))
cB=Counter(archgroup(r["arch"]) for r in fm if re.fullmatch(r'\dB',str(r.get("chrom_lab"))))
xx=np.arange(len(ARCH))
ax.bar(xx-0.19,[cA[a] for a in ARCH],width=0.36,color=CA,label="A",zorder=3)
ax.bar(xx+0.19,[cB[a] for a in ARCH],width=0.36,color=CB,label="B",zorder=3)
ax.set_xticks(xx); ax.set_xticklabels(ARCH); ax.set_ylabel("NLR genes (n)")
ax.set_xlabel("Domain architecture"); ax.set_ylim(0,70); ax.legend(fontsize=5.2,ncol=2,loc="upper left")
tab=np.array([[cA[a] for a in ARCH[:4]],[cB[a] for a in ARCH[:4]]])
c2,p2,_,_=stats.chi2_contingency(tab)
ax.text(0.98,0.995,f"$\\chi^2$ = {c2:.1f}, $P$ = {p2:.3f}",transform=ax.transAxes,
        ha="right",va="top",fontsize=4.9,color=MUTED)
panel(ax,"d",dx=-0.20)

fig.savefig(f"{ROOT}/figures/Figure4.pdf"); fig.savefig(f"{ROOT}/figures/Figure4.png")
print("Figure 4 saved")
print(f"A: {nA} NLRs / {LA:.0f} Mb = {nA/LA*100:.1f} per 100 Mb")
print(f"B: {nB} NLRs / {LB:.0f} Mb = {nB/LB*100:.1f} per 100 Mb")
print(f"homoeolog corr r={r:.3f} P={pv:.2e}")
print(f"density chi2={chi.statistic:.2f} P={chi.pvalue:.4f}")
print(f"architecture A vs B chi2={c2:.2f} P={p2:.4f}")
print("A arch:",dict(cA)); print("B arch:",dict(cB))
