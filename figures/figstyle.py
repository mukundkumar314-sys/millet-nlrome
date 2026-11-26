import matplotlib as mpl, matplotlib.pyplot as plt, json
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
mpl.use("Agg")
MM=1/25.4
W2=180*MM; W1=88*MM; W15=120*MM          # Nature Comms column widths
PAL=["#D55E00","#0072B2","#E69F00","#009E73","#56B4E9","#CC79A7","#999933","#332288","#8C564B","#44AA99","#882255"]
SP=["Ecoracana","Pmiliaceum","Sitalica","Sviridis","Phallii","Dexilis","Ecrusgalli","Sbicolor","Zmays","Osativa","Bdistachyon"]
COL=dict(zip(SP,PAL))
LAB={"Ecoracana":"E. coracana","Pmiliaceum":"P. miliaceum","Sitalica":"S. italica","Sviridis":"S. viridis","Phallii":"P. hallii",
     "Dexilis":"D. exilis","Ecrusgalli":"E. crus-galli","Sbicolor":"S. bicolor","Zmays":"Z. mays",
     "Osativa":"O. sativa","Bdistachyon":"B. distachyon"}
COMMON={"Ecoracana":"finger millet","Pmiliaceum":"proso millet","Sitalica":"foxtail millet","Sviridis":"green millet",
        "Phallii":"Hall's panicgrass","Dexilis":"white fonio","Ecrusgalli":"barnyard millet",
        "Sbicolor":"sorghum","Zmays":"maize","Osativa":"rice","Bdistachyon":"purple false brome"}
CLADE={"Ecoracana":"Chloridoideae","Pmiliaceum":"Panicoideae","Sitalica":"Panicoideae","Sviridis":"Panicoideae",
       "Phallii":"Panicoideae","Dexilis":"Panicoideae","Ecrusgalli":"Panicoideae",
       "Sbicolor":"Panicoideae","Zmays":"Panicoideae","Osativa":"Oryzoideae","Bdistachyon":"Pooideae"}
# ordinal grey ramp for NLR architecture + one accent for ID-NLRs
ARCH_ORDER=["N","CN","NL","CNL","ID-NLR"]
ARCH_COL={"N":"#d9d9d9","CN":"#b0b0b0","NL":"#808080","CNL":"#4d4d4d","ID-NLR":"#D55E00"}
INK="#1a1a1a"; INK2="#4d4d4d"; MUTED="#808080"; GRID="#e6e6e6"

def setup():
    mpl.rcParams.update({
        "font.family":"sans-serif","font.sans-serif":["Liberation Sans","DejaVu Sans"],
        "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
        "font.size":6,"axes.labelsize":6.5,"axes.titlesize":7,
        "xtick.labelsize":5.5,"ytick.labelsize":5.5,"legend.fontsize":5.5,
        "axes.linewidth":0.5,"xtick.major.width":0.5,"ytick.major.width":0.5,
        "xtick.major.size":2,"ytick.major.size":2,"xtick.minor.size":1,"ytick.minor.size":1,
        "axes.edgecolor":INK2,"xtick.color":INK2,"ytick.color":INK2,
        "axes.labelcolor":INK,"text.color":INK,
        "axes.spines.top":False,"axes.spines.right":False,
        "axes.grid":False,"legend.frameon":False,"legend.handlelength":1.2,
        "legend.handletextpad":0.4,"legend.columnspacing":0.9,"legend.labelspacing":0.3,
        "figure.dpi":300,"savefig.dpi":600,"savefig.bbox":"tight","savefig.pad_inches":0.01,
        "lines.linewidth":1.0,"patch.linewidth":0,
    })

def panel(ax,letter,dx=-0.085,dy=1.045):
    ax.text(dx,dy,letter,transform=ax.transAxes,fontsize=8,fontweight="bold",
            va="top",ha="left",color=INK)

def ygrid(ax,axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis,color=GRID,lw=0.4,zorder=0)

ORD11=["Bdistachyon","Osativa","Ecoracana","Dexilis","Pmiliaceum","Phallii","Sviridis","Sitalica","Ecrusgalli","Sbicolor","Zmays"]

def load():
    d={}
    d["rows"]=json.load(open(f"{ROOT}/nlr/nlr_classified.json"))
    d["clusters"]=json.load(open(f"{ROOT}/nlr/clusters.json"))
    d["pan"]=json.load(open(f"{ROOT}/nlr/pan_nlrome.json"))
    d["sizes"]=json.load(open(f"{ROOT}/nlr/genome_sizes.json"))
    d["genes"]=json.load(open(f"{ROOT}/nlr/all_genes.json"))
    return d

def archgroup(a):
    return "ID-NLR" if a.endswith("X") else a
