import re,json,os
from scipy.stats import chi2
def lnL(p):
    if not os.path.exists(p): return None,None
    t=open(p).read()
    m=re.search(r'lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?[\d.]+)',t)
    if not m: return None,None
    return float(m.group(2)), int(m.group(1))
def beb(p):
    if not os.path.exists(p): return []
    t=open(p).read()
    m=re.search(r'Bayes Empirical Bayes \(BEB\).*?\n\n(.*?)\n\n',t,re.S)
    if not m: return []
    out=[]
    for L in m.group(1).split("\n"):
        g=re.match(r'\s*(\d+)\s+([A-Z*])\s+([\d.]+)(\*{0,2})',L)
        if g: out.append({"site":int(g.group(1)),"aa":g.group(2),"pp":float(g.group(3)),"star":g.group(4)})
    return out
res={}
for og in ["OG0021","OG0025","OG0083","OG0213"]:
    d={}
    for mod in ["M1a","M2a","M7","M8"]:
        l,np_=lnL(f"run/{og}/{mod}/out.txt"); d[mod]={"lnL":l,"np":np_}
    r={}
    for a,b,df in [("M1a","M2a",2),("M7","M8",2)]:
        la,lb=d[a]["lnL"],d[b]["lnL"]
        if la is None or lb is None: r[f"{a}v{b}"]=None; continue
        D=2*(lb-la); p=chi2.sf(D,df) if D>0 else 1.0
        r[f"{a}v{b}"]={"2dlnL":D,"df":df,"p":p}
    res[og]={"lnL":d,"lrt":r,"beb_M2a":beb(f"run/{og}/M2a/out.txt"),"beb_M8":beb(f"run/{og}/M8/out.txt")}
json.dump(res,open("codeml_results.json","w"),indent=1)
for og,v in res.items():
    print(f"\n=== {og} ===")
    for k,x in v["lrt"].items():
        print(f"  {k}: 2dlnL={x['2dlnL']:.2f} p={x['p']:.3g}" if x else f"  {k}: incomplete")
    for mod in ["beb_M2a","beb_M8"]:
        s=[b for b in v[mod] if b["star"]]
        print(f"  {mod}: {len(v[mod])} sites listed, {len(s)} with posterior >0.95")
