import sys; sys.path.insert(0,f'{ROOT}/scripts')
from figstyle import *
import numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Ellipse
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
setup()
BLUE="#2B6CB0"; GREY="#9AA0A6"; PINK="#D5547F"; GREEN="#1B8A5A"; INK="#1a1a1a"
PALEG="#EAF4EE"; MUT="#6b7280"; DASH="#B9BEC6"; FADE="#8A9099"

fig=plt.figure(figsize=(W2, W2*0.56))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,56); ax.axis("off")

def arw(p1,p2,c=INK,lw=0.9,rad=0.0,ls="-"):
    ax.add_patch(FancyArrowPatch(p1,p2,connectionstyle=f"arc3,rad={rad}",color=c,lw=lw,
        linestyle=ls,arrowstyle="-|>",mutation_scale=7,zorder=5,shrinkA=2,shrinkB=2))

def chain(x,y,segs,h=2.5,gap=0.45,fs=5.2):
    cx=x
    for lab,w,fc,ec,dsh in segs:
        kw=dict(boxstyle="round,pad=0,rounding_size=0.45",zorder=3)
        if dsh: ax.add_patch(FancyBboxPatch((cx,y-h/2),w,h,fc="white",ec=DASH,lw=0.9,ls=(0,(2,1.5)),**kw))
        else:   ax.add_patch(FancyBboxPatch((cx,y-h/2),w,h,fc=fc,ec=ec,lw=0.6,**kw))
        if lab: ax.text(cx+w/2,y,lab,ha="center",va="center",fontsize=fs,zorder=4,
                        color=DASH if dsh else "white",fontweight="bold" if fc==PINK else "normal")
        cx+=w+gap
    return cx-gap

def pool(x,y,w,h,title,sub,hl):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0,rounding_size=0.7",
        fc=PALEG,ec=GREEN,lw=0.9,zorder=2))
    for i,xx in enumerate(np.linspace(x+2.0,x+w-2.0,7)):
        on = i==hl
        ax.add_patch(FancyBboxPatch((xx-0.8,y+h-2.2),1.6,1.4,
            boxstyle="round,pad=0,rounding_size=0.3",
            fc=PINK if on else "#F2C9D7", ec=PINK if on else "#E0AABF",
            lw=0.9 if on else 0.4, zorder=3))
    ax.text(x+w/2,y+h-3.1,title,ha="center",va="top",fontsize=5.4,color=INK,fontweight="bold")
    ax.text(x+w/2,y+h-4.7,sub,ha="center",va="top",fontsize=4.7,color=MUT,linespacing=1.35)

def fungus(x,y,lab):
    ax.plot([x-3.0,x-1.1,x+0.15],[y+2.4,y+1.9,y+1.2],color="#B59A7A",lw=1.4,solid_capstyle="round",zorder=3)
    ax.add_patch(Ellipse((x-3.4,y+2.6),2.1,1.4,fc="#D8C3A5",ec="#B59A7A",lw=0.5,zorder=3))
    ax.add_patch(FancyBboxPatch((x-0.5,y+0.05),1.0,1.2,boxstyle="round,pad=0,rounding_size=0.25",
        fc="#B59A7A",ec="#8E7355",lw=0.5,zorder=3))
    ax.text(x-4.9,y+2.6,lab,fontsize=5.3,color=INK,va="center",ha="right",fontweight="bold",style="italic")

def cross(x,y,s=0.9,c="#C0463F",lw=1.4):
    ax.plot([x-s,x+s],[y-s,y+s],color=c,lw=lw,zorder=6,solid_capstyle="round")
    ax.plot([x-s,x+s],[y+s,y-s],color=c,lw=lw,zorder=6,solid_capstyle="round")

for x0,ttl,sub,col in [(2.0,"Oryza sativa","complete recognition circuit",BLUE),
                       (51.0,"Eleusine coracana","sensor slot never filled",FADE)]:
    ax.add_patch(FancyBboxPatch((x0,8.6),47,45.6,boxstyle="round,pad=0,rounding_size=1.0",
        fc="white",ec="#D7DBE0",lw=0.9,zorder=1))
    ax.text(x0+1.6,52.6,ttl,fontsize=7.4,style="italic",fontweight="bold",color=INK,va="top")
    ax.text(x0+45.4,52.2,sub,fontsize=5.8,color=col,va="top",fontweight="bold",ha="right")
    ax.add_patch(Rectangle((x0+0.9,45.6),45.2,1.2,fc="#D3E8D8",ec="none",zorder=2))

def panel(X, host):
    rice = host=="rice"
    fungus(X+38.0,46.1,"Magnaporthe oryzae" if rice else "Eleusine-infecting lineage")
    arw((X+38.2,45.3),(X+35.5,43.3),lw=0.9)
    ax.add_patch(Circle((X+34.6,42.6),0.8,fc=PINK,ec="#B23E68",lw=0.5,zorder=4))
    if rice:
        ax.text(X+33.4,42.6,"AVR-Pik",fontsize=5.4,ha="right",va="center",
                fontweight="bold",color=INK)
    else:
        ax.text(X+33.4,42.6,"AVR-Pik-like effector (33 of 33 genomes)",fontsize=5.4,
                ha="right",va="center",fontweight="bold",color=INK)

    pool(X+2.5,29.5,17,8.2,"standalone HMA family",
         "OsHIPP19 and OsHIPP20\nare members" if rice else
         "95 members; Clade A pair at\n80–85% identity to OsHIPP19/20",3)

    if rice:
        arw((X+34.4,41.6),(X+29.5,24.4),lw=0.9,rad=0.16)
        ax.text(X+36.0,33.0,"effector binds\nthe decoy",fontsize=5.1,color=INK,
                ha="left",va="center",linespacing=1.35)
    else:
        arw((X+34.0,41.5),(X+11.2,38.4),lw=0.9,rad=0.06)
        ax.text(X+21.5,38.6,"binds the host HMA target",fontsize=5.1,color=INK,
                ha="center",va="center")

    ax.text(X+20.6,28.5,("one member captured\ninto an NLR" if rice else
                         "no member captured\ninto an NLR"),
            fontsize=5.1,color=PINK if rice else MUT,ha="right",fontweight="bold",
            va="top",linespacing=1.35)
    if rice:
        arw((X+19.8,30.2),(X+24.4,23.7),c=PINK,lw=1.2,rad=-0.22)
    else:
        ax.add_patch(FancyArrowPatch((X+19.8,30.2),(X+24.4,23.7),connectionstyle="arc3,rad=-0.22",
            color=DASH,lw=1.1,linestyle=(0,(2,1.5)),arrowstyle="-|>",mutation_scale=7,
            zorder=5,shrinkA=2,shrinkB=2))
        cross(X+23.3,27.3)

    if rice:
        chain(X+20.5,22.9,[("CC",4.0,BLUE,"#1E4E8C",False),("HMA",4.4,PINK,"#B23E68",False),
                           ("NB-ARC",6.4,BLUE,"#1E4E8C",False),("LRR",4.8,BLUE,"#1E4E8C",False)])
        ax.text(X+20.5,20.4,"Pik-1 sensor",fontsize=5.4,color=BLUE,
                fontweight="bold",va="top")
        chain(X+20.5,16.4,[("CC",4.0,GREY,"#7B8186",False),("NB-ARC",6.4,GREY,"#7B8186",False),
                           ("LRR",4.8,GREY,"#7B8186",False)])
        ax.text(X+20.5,13.9,"Pik-2 helper",fontsize=5.4,color=MUT,fontweight="bold",va="top")
        arw((X+36.0,21.5),(X+36.0,17.8),c=GREEN,lw=1.8)
        ax.text(X+37.1,19.65,"immunity",fontsize=5.8,color=GREEN,fontweight="bold",
                va="center",ha="left")
        cols=[BLUE,GREY]; cap="head-to-head, 1,100 bp"
    else:
        chain(X+20.5,22.9,[("CC",4.0,"white",DASH,True),("",4.4,"white",DASH,True),
                           ("NB-ARC",6.4,"white",DASH,True),("LRR",4.8,"white",DASH,True)])
        ax.text(X+20.5,20.4,"sensor NLR — absent",fontsize=5.4,color=FADE,
                fontweight="bold",va="top")
        chain(X+20.5,16.4,[("CC",4.0,GREEN,"#146B45",False),("NB-ARC",6.4,GREEN,"#146B45",False),
                           ("LRR",4.8,GREEN,"#146B45",False)])
        ax.text(X+20.5,13.9,"Pik-2 helper orthologue — present",fontsize=5.4,color=GREEN,
                fontweight="bold",va="top")
        ax.add_patch(FancyArrowPatch((X+36.0,21.5),(X+36.0,17.8),color=DASH,lw=1.5,
            linestyle=(0,(2,1.5)),arrowstyle="-|>",mutation_scale=7,zorder=5,shrinkA=2,shrinkB=2))
        ax.text(X+37.1,19.65,"no activation",fontsize=5.6,color=FADE,fontweight="bold",
                va="center",ha="left")
        cols=[DASH,GREEN]; cap="7 head-to-head pairs, none with a decoy"

    ax.plot([X+20.5,X+31.5],[11.4,11.4],color=INK,lw=0.7)
    for xa,d,c in [(X+22.9,-1,cols[0]),(X+29.1,1,cols[1])]:
        ax.add_patch(mpatches.FancyArrow(xa,11.4,d*1.9,0,width=0.5,head_width=1.2,head_length=0.85,
            fc=c,ec="none",zorder=4,length_includes_head=True))
    ax.text(X+26.0,10.1,cap,fontsize=4.9,color=MUT,ha="center",va="top")

panel(2.0,"rice"); panel(51.0,"millet")

ax.add_patch(FancyBboxPatch((2.0,1.2),96,5.8,boxstyle="round,pad=0,rounding_size=0.9",
    fc="#F7F9FB",ec="#D7DBE0",lw=0.9,zorder=1))
items=[("Pathogen effector","present throughout the Eleusine lineage",True),
       ("Host HMA target","retained at 80–85% identity",True),
       ("Integrated decoy","never captured into an NLR",False)]
for i,(t,s,ok) in enumerate(items):
    x=6.0+i*31.8
    if ok:
        ax.plot([x-2.4,x-1.6],[4.9,4.2],color=GREEN,lw=1.3,solid_capstyle="round")
        ax.plot([x-1.6,x-0.3],[4.2,5.7],color=GREEN,lw=1.3,solid_capstyle="round")
    else: cross(x-1.4,4.9,0.75,"#9AA0A6",1.3)
    ax.text(x+0.6,5.0,t,fontsize=6.2,fontweight="bold",color=GREEN if ok else FADE,va="center")
    ax.text(x+0.6,2.8,s,fontsize=5.2,color=MUT,va="center")
    if i<2: ax.plot([x+29.0,x+29.0],[2.1,6.2],color="#D7DBE0",lw=0.8)

plt.savefig(f"{ROOT}/figures/Figure_model.png",dpi=600,bbox_inches="tight",facecolor="white")
plt.savefig(f"{ROOT}/figures/Figure_model.pdf",bbox_inches="tight",facecolor="white")
print("saved")
