import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt, re, json
from collections import Counter
from scipy import stats
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup(); d=load(); rows=d["rows"]; cl=d["clusters"]
CA="#D55E00"; CB="#0072B2"

# finger millet
CHL={}
for ln in open(f"{ROOT}/gff/Ecoracana.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\r\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': CHL[f[2]]=int(f[8])
fm=[r for r in rows if r["species"]=="Ecoracana"]
# proso
PS=json.load(open(f"{ROOT}/nlr/proso_subgenomes.json")); sub=PS["subgenome"]; homo=PS["pairs"]
PL={}
for ln in open(f"{ROOT}/gff/Pmiliaceum.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\r\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': PL[f[0]]=int(f[8])
pm=[r for r in rows if r["species"]=="Pmiliaceum"]

fig=plt.figure(figsize=(W2,W2*0.60))
gs=fig.add_gridspec(2,12,height_ratios=[1.16,0.84],hspace=0.48,wspace=3.0,
                    left=0.075,right=0.985,top=0.925,bottom=0.095)

# a: finger millet ideogram
ax=fig.add_subplot(gs[0,0:12])
P=[str(i) for i in range(1,10)]; maxL=max(CHL.values()); BW=0.30
g=cl["Ecoracana"]["g"]*1000
for i,p in enumerate(P):
    for k,(s_,c) in enumerate([("A",CA),("B",CB)]):
        ch=f"{p}{s_}"; L=CHL[ch]/1e6; y=i+(-0.19 if k==0 else 0.19)
        ax.add_patch(plt.Rectangle((0,y-BW/2),L,BW,facecolor="#f0f0f0",edgecolor="#c0c0c0",lw=0.35,zorder=2))
        pos=sorted((r["start"]+r["end"])//2 for r in fm if str(r.get("chrom_lab"))==ch)
        cur=[pos[0]] if pos else []
        for q in range(1,len(pos)):
            if pos[q]-pos[q-1]<=g: cur.append(pos[q])
            else:
                if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.22,lw=0,zorder=3))
                cur=[pos[q]]
        if len(cur)>=3: ax.add_patch(plt.Rectangle((cur[0]/1e6,y-BW/2),max((cur[-1]-cur[0])/1e6,0.6),BW,facecolor=c,alpha=0.22,lw=0,zorder=3))
        for x0 in pos: ax.plot([x0/1e6]*2,[y-BW/2,y+BW/2],color=c,lw=0.45,zorder=4,alpha=0.92)
        ax.text(-1.5,y,ch,ha="right",va="center",fontsize=5.0,color=INK)
        ax.text(L+1.2,y,str(len(pos)),ha="left",va="center",fontsize=4.6,color=c)
ax.set_xlim(-7,maxL/1e6+7); ax.set_ylim(-0.65,len(P)-0.35)
ax.set_xlabel("Chromosome position (Mb)"); ax.set_yticks([]); ax.invert_yaxis()
for sp_ in ["left","right","top"]: ax.spines[sp_].set_visible(False)
hs=[plt.Rectangle((0,0),1,1,facecolor=CA),plt.Rectangle((0,0),1,1,facecolor=CB)]
ax.legend(hs,["A subgenome","B subgenome"],loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=5.4)
ax.text(0.995,0.015,"$Eleusine$ $coracana$  (Chloridoideae) · shaded block = NLR cluster (≥3 genes)",
        transform=ax.transAxes,ha="right",fontsize=4.7,color=MUTED)
panel(ax,"a",dx=-0.052,dy=1.085)

# b: finger millet homoeolog counts
ax=fig.add_subplot(gs[1,0:3]); ygrid(ax)
A=[sum(1 for r in fm if str(r.get("chrom_lab"))==f"{p}A") for p in P]
B=[sum(1 for r in fm if str(r.get("chrom_lab"))==f"{p}B") for p in P]
x=np.arange(len(P))
ax.bar(x-0.19,A,width=0.36,color=CA,label="A",zorder=3); ax.bar(x+0.19,B,width=0.36,color=CB,label="B",zorder=3)
ax.set_xticks(x); ax.set_xticklabels(P,fontsize=5); ax.set_xlabel("Homoeologous group")
ax.set_ylabel("NLR genes (n)"); ax.set_ylim(0,34); ax.legend(fontsize=5.2,ncol=2,loc="upper left")
r1,p1=stats.pearsonr(A,B)
ax.text(0.98,0.995,f"$r$ = {r1:.2f}, $P$ = {p1:.1e}",transform=ax.transAxes,ha="right",va="top",fontsize=4.8,color=MUTED)
ax.set_title("$E.$ $coracana$",fontsize=5.6,pad=2)
panel(ax,"b",dx=-0.28)

# c: proso homoeolog counts
ax=fig.add_subplot(gs[1,3:6]); ygrid(ax)
hp=sorted(homo)
pA=[sum(1 for r in pm if str(r.get("chrom_lab"))==(a if sub[a]=="A" else b)) for a,b,n in hp]
pB=[sum(1 for r in pm if str(r.get("chrom_lab"))==(b if sub[a]=="A" else a)) for a,b,n in hp]
x=np.arange(len(hp))
ax.bar(x-0.19,pA,width=0.36,color=CA,zorder=3); ax.bar(x+0.19,pB,width=0.36,color=CB,zorder=3)
ax.set_xticks(x); ax.set_xticklabels([str(i+1) for i in x],fontsize=5); ax.set_xlabel("Homoeologous group")
ax.set_ylabel("NLR genes (n)"); ax.set_ylim(0,34)
r2,p2=stats.pearsonr(pA,pB)
ax.text(0.98,0.995,f"$r$ = {r2:.2f}, $P$ = {p2:.1e}",transform=ax.transAxes,ha="right",va="top",fontsize=4.8,color=MUTED)
ax.set_title("$P.$ $miliaceum$",fontsize=5.6,pad=2)
panel(ax,"c",dx=-0.28)

# d: subgenome size vs NLR, both species
ax=fig.add_subplot(gs[1,6:9]); ygrid(ax)
LA=sum(CHL[f"{p}A"] for p in P)/1e6; LB=sum(CHL[f"{p}B"] for p in P)/1e6
pLA=sum(PL[c] for c,v in sub.items() if v=="A")/1e6; pLB=sum(PL[c] for c,v in sub.items() if v=="B")/1e6
dens=[[sum(A)/LA*100,sum(B)/LB*100],[sum(pA)/pLA*100,sum(pB)/pLB*100]]
xx=np.arange(2)
for k in range(2):
    ax.bar(xx[k]-0.19,dens[k][0],width=0.36,color=CA,zorder=3)
    ax.bar(xx[k]+0.19,dens[k][1],width=0.36,color=CB,zorder=3)
    ax.text(xx[k]-0.19,dens[k][0]+0.7,f"{dens[k][0]:.1f}",ha="center",fontsize=4.8,color=INK)
    ax.text(xx[k]+0.19,dens[k][1]+0.7,f"{dens[k][1]:.1f}",ha="center",fontsize=4.8,color=INK)
ax.set_xticks(xx); ax.set_xticklabels(["$E.$ $coracana$","$P.$ $miliaceum$"],fontsize=5.2)
ax.set_ylabel("NLR density (per 100 Mb)"); ax.set_ylim(0,52)
c1=stats.chisquare([sum(A),sum(B)],[(sum(A)+sum(B))*LA/(LA+LB),(sum(A)+sum(B))*LB/(LA+LB)])
c2=stats.chisquare([sum(pA),sum(pB)],[(sum(pA)+sum(pB))*pLA/(pLA+pLB),(sum(pA)+sum(pB))*pLB/(pLA+pLB)])
ax.text(0.5,1.01,f"$P$ = {c1.pvalue:.2f}          $P$ = {c2.pvalue:.2f}",transform=ax.transAxes,
        ha="center",fontsize=4.7,color=MUTED)
panel(ax,"d",dx=-0.30)

# e: architecture A vs B, both species
ax=fig.add_subplot(gs[1,9:12]); ygrid(ax)
ARCH=["N","CN","NL","CNL"]
cA=Counter(archgroup(r["arch"]) for r in fm if re.fullmatch(r'\dA',str(r.get("chrom_lab"))))
cB=Counter(archgroup(r["arch"]) for r in fm if re.fullmatch(r'\dB',str(r.get("chrom_lab"))))
qA=Counter(archgroup(r["arch"]) for r in pm if sub.get(str(r.get("chrom_lab")))=="A")
qB=Counter(archgroup(r["arch"]) for r in pm if sub.get(str(r.get("chrom_lab")))=="B")
xx=np.arange(len(ARCH))
ax.bar(xx-0.30,[cA[a] for a in ARCH],width=0.19,color=CA,label="$E. c.$ A",zorder=3)
ax.bar(xx-0.10,[cB[a] for a in ARCH],width=0.19,color=CB,label="$E. c.$ B",zorder=3)
ax.bar(xx+0.10,[qA[a] for a in ARCH],width=0.19,color=CA,alpha=0.52,label="$P. m.$ A",zorder=3)
ax.bar(xx+0.30,[qB[a] for a in ARCH],width=0.19,color=CB,alpha=0.52,label="$P. m.$ B",zorder=3)
ax.set_xticks(xx); ax.set_xticklabels(ARCH); ax.set_xlabel("Domain architecture")
ax.set_ylabel("NLR genes (n)"); ax.set_ylim(0,78); ax.legend(fontsize=4.7,ncol=2,loc="upper left")
x1,pp1,_,_=stats.chi2_contingency(np.array([[cA[a] for a in ARCH],[cB[a] for a in ARCH]]))
x2,pp2,_,_=stats.chi2_contingency(np.array([[qA[a] for a in ARCH],[qB[a] for a in ARCH]]))
ax.text(0.98,0.995,f"$P$ = {pp1:.2f} | {pp2:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=4.7,color=MUTED)
panel(ax,"e",dx=-0.30)

fig.savefig(f"{ROOT}/figures/Figure4.pdf"); fig.savefig(f"{ROOT}/figures/Figure4.png")
print(f"Figure 4 saved")
print(f"Ecoracana A {sum(A)}/{LA:.0f}Mb={sum(A)/LA*100:.1f}  B {sum(B)}/{LB:.0f}Mb={sum(B)/LB*100:.1f}  r={r1:.2f} P={p1:.1e} chi2P={c1.pvalue:.3f}")
print(f"Pmiliaceum A {sum(pA)}/{pLA:.0f}Mb={sum(pA)/pLA*100:.1f}  B {sum(pB)}/{pLB:.0f}Mb={sum(pB)/pLB*100:.1f}  r={r2:.2f} P={p2:.1e} chi2P={c2.pvalue:.3f}")
print(f"arch chi2 P: Ecor {pp1:.3f}  Pmil {pp2:.3f}")
