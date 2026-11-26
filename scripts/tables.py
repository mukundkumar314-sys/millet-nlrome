import json, re
from collections import Counter, defaultdict
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")

NLR=f"{ROOT}/nlr"
rows=json.load(open(f"{NLR}/nlr_classified.json"))
pan=json.load(open(f"{NLR}/pan_nlrome.json"))
cl=json.load(open(f"{NLR}/clusters.json"))
sz=json.load(open(f"{NLR}/genome_sizes.json"))
genes=json.load(open(f"{NLR}/all_genes.json"))
hma_sa=json.load(open(f"{ROOT}/refnlr/hma_standalone.json"))
pik=json.load(open(f"{ROOT}/refnlr/pik_data.json"))

ORD=["Ecoracana","Pmiliaceum","Sitalica","Sviridis","Phallii","Dexilis","Ecrusgalli","Sbicolor","Zmays","Osativa","Bdistachyon"]
FULL={"Ecoracana":"Eleusine coracana","Pmiliaceum":"Panicum miliaceum","Sitalica":"Setaria italica","Sviridis":"Setaria viridis",
      "Phallii":"Panicum hallii","Dexilis":"Digitaria exilis","Ecrusgalli":"Echinochloa crus-galli",
      "Sbicolor":"Sorghum bicolor","Zmays":"Zea mays","Osativa":"Oryza sativa",
      "Bdistachyon":"Brachypodium distachyon"}
COMMON={"Ecoracana":"finger millet","Pmiliaceum":"proso millet","Sitalica":"foxtail millet","Sviridis":"green foxtail",
        "Phallii":"Hall's panicgrass","Dexilis":"white fonio","Ecrusgalli":"barnyard millet",
        "Sbicolor":"sorghum","Zmays":"maize","Osativa":"rice","Bdistachyon":"purple false brome"}
CLADE={"Ecoracana":"Chloridoideae","Pmiliaceum":"Panicoideae","Sitalica":"Panicoideae","Sviridis":"Panicoideae",
       "Phallii":"Panicoideae","Dexilis":"Panicoideae","Ecrusgalli":"Panicoideae",
       "Sbicolor":"Panicoideae","Zmays":"Panicoideae","Osativa":"Oryzoideae","Bdistachyon":"Pooideae"}
ASM={"Ecoracana":("Eleusine_coracana_v1.0","GCA_032690845.1","NCBI GenBank"),
     "Pmiliaceum":("Pm_0390_v2","GCA_003046395.2","NCBI GenBank"),
     "Sitalica":("Setaria_italica_v2.0","GCA_000263155.2","Ensembl Plants 60"),
     "Sviridis":("Setaria_viridis_v2.0","—","Ensembl Plants 60"),
     "Phallii":("PhalliiHAL_v2.1","—","Ensembl Plants 60"),
     "Dexilis":("Fonio_CM05836","—","Ensembl Plants 60"),
     "Ecrusgalli":("ec_v3","GCA_020466025.1","Ensembl Plants 60"),
     "Sbicolor":("Sorghum_bicolor_NCBIv3","—","Ensembl Plants 60"),
     "Zmays":("Zm-B73-REFERENCE-NAM-5.0","—","Ensembl Plants 60"),
     "Osativa":("IRGSP-1.0","—","Ensembl Plants 60"),
     "Bdistachyon":("Brachypodium_distachyon_v3.0","—","Ensembl Plants 60")}

wb=Workbook(); wb.remove(wb.active)
HFILL=PatternFill("solid",fgColor="D9E2F3")
THIN=Side(style="thin",color="BFBFBF")
def sheet(name,headers,data,widths=None,note=None):
    ws=wb.create_sheet(name)
    r0=1
    if note:
        ws.cell(1,1,note).font=Font(name="Arial",size=9,italic=True,color="555555")
        r0=3
    for j,h in enumerate(headers,1):
        c=ws.cell(r0,j,h); c.font=Font(name="Arial",size=10,bold=True)
        c.fill=HFILL; c.alignment=Alignment(horizontal="center",vertical="center",wrap_text=True)
        c.border=Border(bottom=THIN)
    for i,row in enumerate(data,r0+1):
        for j,v in enumerate(row,1):
            c=ws.cell(i,j,v); c.font=Font(name="Arial",size=10)
    for j,h in enumerate(headers,1):
        w=(widths[j-1] if widths else max(11,min(34,len(str(h))+4)))
        ws.column_dimensions[get_column_letter(j)].width=w
    ws.freeze_panes=ws.cell(r0+1,1)
    return ws

# ---- S1 ----
d=[]
for s in ORD:
    n=sum(1 for r in rows if r["species"]==s)
    d.append([FULL[s],COMMON[s],CLADE[s],ASM[s][0],ASM[s][1],ASM[s][2],
              round(sz[s]/1e6,1),len(genes[s]),n,round(n/(sz[s]/1e6)*100,1)])
sheet("S1_Genomes",
 ["Species","Common name","Subfamily","Assembly","Accession","Source","Assembly size (Mb)",
  "Annotated protein-coding genes (n)","NLR genes (n)","NLR density (per 100 Mb)"],d,
 [24,20,15,30,18,20,16,18,13,15],
 "Supplementary Table S1 | Genome assemblies and annotations analysed in this study.")

# ---- S2 ----
d=[]
for r in sorted(rows,key=lambda r:(ORD.index(r["species"]),str(r.get("chrom_lab",r["chrom"])),r["start"])):
    d.append([FULL[r["species"]],r["gene"],r["prot"],str(r.get("chrom_lab",r["chrom"])),
              r["start"],r["end"],r["strand"],r["length"],r["arch"],
              "yes" if r["has_CC"] else "no","yes" if r["has_LRR"] else "no",
              r["n_id"],r["ids"],r["te_assoc"],round(r["nb_score"],1),r["nb_eval"]])
sheet("S2_NLR_catalogue",
 ["Species","Gene ID","Protein ID","Chromosome","Start","End","Strand","Protein length (aa)",
  "Architecture","CC (Rx_N)","LRR","Integrated domains (n)","Integrated domain families",
  "TE-associated domains","NB-ARC bit score","NB-ARC E-value"],d,
 [22,24,22,13,12,12,8,13,12,11,8,13,30,22,13,13],
 "Supplementary Table S2 | Catalogue of 3,012 NLR genes identified across ten Poaceae genomes.")

# ---- S3 (formulas over S2) ----
ARCHS=["N","CN","NL","CNL","NX","CNX","NLX","CNLX"]
ws=wb.create_sheet("S3_Architecture")
ws.cell(1,1,"Supplementary Table S3 | NLR domain architecture composition per species. Values are COUNTIFS formulas over sheet S2.").font=Font(name="Arial",size=9,italic=True,color="555555")
hdr=["Species"]+ARCHS+["Total NLRs","ID-NLRs (n)","ID-NLRs (%)"]
for j,h in enumerate(hdr,1):
    c=ws.cell(3,j,h); c.font=Font(name="Arial",size=10,bold=True); c.fill=HFILL
    c.alignment=Alignment(horizontal="center",wrap_text=True); c.border=Border(bottom=THIN)
N=len(rows)+3
for i,s in enumerate(ORD,4):
    ws.cell(i,1,FULL[s]).font=Font(name="Arial",size=10)
    for j,a in enumerate(ARCHS,2):
        ws.cell(i,j,f'=COUNTIFS(S2_NLR_catalogue!$A$4:$A${N},$A{i},S2_NLR_catalogue!$I$4:$I${N},"{a}")').font=Font(name="Arial",size=10)
    tc=get_column_letter(len(ARCHS)+2)
    ws.cell(i,len(ARCHS)+2,f'=SUM(B{i}:{get_column_letter(len(ARCHS)+1)}{i})').font=Font(name="Arial",size=10)
    ws.cell(i,len(ARCHS)+3,f'=SUM(F{i}:I{i})').font=Font(name="Arial",size=10)
    c=ws.cell(i,len(ARCHS)+4,f'=IFERROR({get_column_letter(len(ARCHS)+3)}{i}/{tc}{i},0)')
    c.font=Font(name="Arial",size=10); c.number_format="0.0%"
for j in range(1,len(hdr)+1):
    ws.column_dimensions[get_column_letter(j)].width=(24 if j==1 else 11)
ws.freeze_panes="A4"

# ---- S4 ----
d=[[FULL[s],cl[s]["n"],cl[s]["g"],cl[s]["enr"],cl[s]["nclusters"],cl[s]["in_cluster"],
    cl[s]["pct"]/100,cl[s]["maxsize"],cl[s]["n"]-cl[s]["in_cluster"]] for s in ORD]
ws=sheet("S4_Clustering",
 ["Species","NLR genes (n)","Cluster distance threshold g (kb)","Observed/expected enrichment at g",
  "Clusters (n)","NLRs in clusters (n)","NLRs in clusters (%)","Largest cluster (NLRs)","Singleton NLRs (n)"],d,
 [24,13,18,18,12,14,14,14,14],
 "Supplementary Table S4 | Physical clustering of NLR genes. g was estimated per species as the largest distance at which observed NLR-NLR spacing was >=2-fold enriched over 100 randomisations of annotated genes.")
for i in range(4,4+len(ORD)): ws.cell(i,7).number_format="0.0%"

# ---- S5 ----
idrows=[r for r in rows if r["n_id"]>0]
idc=Counter(); idsp=defaultdict(set)
for r in idrows:
    for i in r["ids"].split(";"):
        if i: idc[i]+=1; idsp[i].add(r["species"])
d=[]
for fam,n in idc.most_common():
    per=[sum(1 for r in idrows if r["species"]==s and fam in r["ids"].split(";")) for s in ORD]
    d.append([fam,n,len(idsp[fam])]+per)
sheet("S5_Integrated_domains",["Integrated domain family (Pfam)","NLR genes (n)","Species (n)"]+[FULL[s] for s in ORD],d,
 [30,12,10]+[17]*len(ORD),
 "Supplementary Table S5 | Integrated domain (ID) families detected in NLR proteins by full Pfam-A scan (gathering thresholds).")

# ---- S6 ----
d=[[r["og"],r["n"],r["nspecies"],r["cat"],", ".join(FULL[x] for x in r["species"]),
    ", ".join(m.split("|")[1] for m in r["members"][:25])+(" ..." if len(r["members"])>25 else "")]
   for r in pan]
sheet("S6_PanNLRome",["Orthogroup","NLR genes (n)","Species (n)","Category","Species","Member gene IDs (first 25)"],d,
 [12,12,11,17,46,60],
 "Supplementary Table S6 | Pan-NLRome orthogroups. core = present in 10/10 genomes; soft-core = 9/10; shell = 2-8; species-specific = 1.")

# ---- S7 ----
CHL={}
for ln in open(f"{ROOT}/gff/Ecoracana.assembly_report.txt"):
    if ln.startswith('#'): continue
    f=ln.rstrip('\n').split('\t')
    if len(f)>=10 and f[1]=='assembled-molecule': CHL[f[2]]=int(f[8])
fm=[r for r in rows if r["species"]=="Ecoracana"]
d=[]
for p in range(1,10):
    for sub in "AB":
        ch=f"{p}{sub}"; L=CHL[ch]
        sel=[r for r in fm if r.get("chrom_lab")==ch]
        c=Counter(r["arch"] for r in sel)
        d.append([ch,sub,round(L/1e6,1),len(sel),round(len(sel)/(L/1e6)*100,1),
                  c.get("CNL",0)+c.get("CNLX",0),c.get("CN",0)+c.get("CNX",0),
                  c.get("NL",0)+c.get("NLX",0),c.get("N",0),
                  sum(1 for r in sel if r["n_id"]>0)])
sheet("S7_Fingermillet_chrom",
 ["Chromosome","Subgenome","Length (Mb)","NLR genes (n)","NLR density (per 100 Mb)",
  "CNL","CN","NL","N","ID-NLRs"],d,[13,12,12,13,15,9,9,9,9,10],
 "Supplementary Table S7 | NLR distribution across the A and B subgenomes of Eleusine coracana (allotetraploid, 2n = 4x = 36).")

# ---- S8 ----
d=[["Rice Pik-1 (sensor)","Oryza sativa","Os11g0688832","11",pik["Os11t0688832-01"]["start"],
    pik["Os11t0688832-01"]["end"],"-","CNL + integrated HMA",1125,"HMA at aa 237-287 (bit score 13.4)"],
   ["Rice Pik-2 (helper)","Oryza sativa","Os11g0689100","11",pik["Os11t0689100-01"]["start"],
    pik["Os11t0689100-01"]["end"],"+","CNL",1044,"head-to-head with Pik-1, 1,100 bp intergenic"],
   ["Pik-1 ortholog","Setaria viridis","SEVIR (TKV99356)","8",pik["TKV99356"]["start"],
    pik["TKV99356"]["end"],"+","CNL",1126,"no integrated HMA"],
   ["Pik-1 ortholog","Setaria viridis","SEVIR (TKW02560)","8",pik["TKW02560"]["start"],
    pik["TKW02560"]["end"],"-","CNL",1027,"no integrated HMA"],
   ["Pik-1 ortholog","Echinochloa crus-galli","scaffold389.6","CH04",pik["scaffold389.6.cds"]["start"],
    pik["scaffold389.6.cds"]["end"],"-","CNL",1021,"no integrated HMA"],
   ["Pik-1 ortholog","Sorghum bicolor","EES09001","5",pik["EES09001"]["start"],
    pik["EES09001"]["end"],"+","CNL",1003,"no integrated HMA"]]
sheet("S8_Pik_HMA",["Locus","Species","Gene","Chromosome","Start","End","Strand",
                    "Architecture","Length (aa)","Note"],d,
 [20,24,22,13,12,12,8,22,12,42],
 "Supplementary Table S8 | The Pik sensor-helper locus and its orthologues. HMA detection used custom profiles built from the integrated HMA domains of rice Pikm1-TS and RGA5.")
d2=[[FULL[s],hma_sa.get(s,0),1 if s=="Osativa" else 0] for s in ORD]
sheet("S9_HMA_proteins",["Species","Standalone HMA proteins (n)","HMA-integrated NLRs (n)"],d2,
 [24,22,22],
 "Supplementary Table S9 | Standalone heavy-metal-associated (HMA) proteins versus HMA-integrated NLRs. Standalone HMA proteins (candidate effector targets) are present in every genome; the integration into an NLR is detected only in rice.")




# ---- S10 head-to-head pairs ----
pairs=json.load(open(f"{ROOT}/nlr/h2h_pairs.json"))
by={r["gene"]+"|"+r["species"]:r for r in rows}
d=[]
for p in pairs:
    a=by.get(p["g1"]+"|"+p["species"]); b=by.get(p["g2"]+"|"+p["species"])
    if not a or not b: continue
    ids=";".join(x for x in [a["ids"],b["ids"]] if x)
    d.append([FULL[p["species"]],str(p["chrom"]),p["g1"],p["a1"],p["g2"],p["a2"],p["gap"],
              "yes" if ids else "no", ids])
d.sort(key=lambda r:(r[0],r[1]))
sheet("S10_SensorHelper_pairs",
 ["Species","Chromosome","Gene 1 (minus strand)","Architecture 1","Gene 2 (plus strand)",
  "Architecture 2","Intergenic distance (bp)","Sensor carries an ID","Integrated domains"],d,
 [24,14,26,13,26,13,16,14,28],
 "Supplementary Table S10 | Head-to-head (divergently transcribed) NLR gene pairs within 10 kb, the genomic signature of sensor-helper NLR pairs. The rice Pik-1/Pik-2 pair is recovered as a positive control.")

# ---- S11 Ka/Ks ----
kk=json.load(open(f"{ROOT}/genome/fm_kaks.json"))
d=[[r["og"],r["a"],r["ca"],r["b"],r["cb"],r["codons"],round(r["dn"],4),round(r["ds"],4),round(r["w"],3),r["arch"]]
   for r in sorted(kk,key=lambda r:-r["w"])]
sheet("S11_KaKs_homoeologs",
 ["Orthogroup","Gene (A subgenome)","Chr A","Gene (B subgenome)","Chr B","Aligned codons",
  "Ka (dN)","Ks (dS)","Ka/Ks","Architecture"],d,
 [12,26,9,26,9,13,11,11,10,12],
 "Supplementary Table S11 | Pairwise Ka/Ks (YN00) between homoeologous NLR gene pairs of the Eleusine coracana A and B subgenomes.")

# ---- S12 HMA genomic loci ----
hl=json.load(open(f"{ROOT}/genome/hma_loci.json"))
fmr=[r for r in rows if r["species"]=="Ecoracana"]
byacc={}
for r in fmr: byacc.setdefault(r["chrom"],[]).append((r["start"],r["end"],r["gene"],r["arch"]))
d=[]
for l in sorted(hl,key=lambda x:-x["bit"]):
    best=None
    for s_,e_,g_,a_ in byacc.get(l["acc"],[]):
        dd=0 if not(l["e"]<s_ or l["s"]>e_) else (s_-l["e"] if l["e"]<s_ else l["s"]-e_)
        if best is None or dd<best[0]: best=(dd,g_,a_)
    d.append([str(l["chrom"]),l["s"],l["e"],l["q"],round(l["pid"],1),round(l["bit"],0),
              f"{l['ev']:.1e}",best[1] if best else "-",best[2] if best else "-",
              best[0] if best else "-"])
sheet("S12_HMA_genomic_loci",
 ["Chromosome","Start","End","Query seed","Identity (%)","Bit score","E-value",
  "Nearest NLR gene","Its architecture","Distance (bp)"],d,
 [13,13,13,13,12,11,12,24,14,14],
 "Supplementary Table S12 | HMA-domain loci detected by tblastn in the Eleusine coracana genome assembly, with distance to the nearest annotated NLR gene. No locus falls inside or within 10 kb of an NLR.")



# ---- S13 blast R-gene orthologues ----
cand=json.load(open(f"{ROOT}/nlr/fm_blastR_candidates.json"))
d=[[c["rgene"],c["rice"],c["og"],c["fm_gene"],c["prot"],c["chrom"],c["start"],c["end"],
    c["arch"],c["length"]] for c in sorted(cand,key=lambda c:(c["chrom"],c["start"]))]
sheet("S13_BlastR_candidates",
 ["Cloned rice blast R gene","Rice NLR locus","Orthogroup","E. coracana gene","Protein ID",
  "Chromosome","Start","End","Architecture","Length (aa)"],d,
 [22,20,12,24,18,13,12,12,13,12],
 "Supplementary Table S13 | Finger millet NLRs that fall in the same orthogroup as a cloned rice blast resistance gene. Candidate genes for targeted work on blast resistance.")

# ---- S14 genomic context ----
gd=json.load(open(f"{ROOT}/nlr/gene_density.json"))
d=[[FULL[s],round(gd[s]["nlr"],0),round(gd[s]["rand"],0),round(gd[s]["ratio"],2),f'{gd[s]["p"]:.2e}']
   for s in ORD if s in gd]
sheet("S14_Gene_density",
 ["Species","NLR neighbourhoods (median genes per 200 kb)","Matched random genes (median)",
  "Ratio","Mann-Whitney P"],d,[24,22,22,10,15],
 "Supplementary Table S14 | Gene density in 200 kb windows centred on NLR loci versus matched random genes.")

# ---- S15 ID classes ----
idc=json.load(open(f"{ROOT}/nlr/id_classes.json"))
d=[[k,v[0],len(v[1]),", ".join(FULL.get(x,x) for x in v[1])]
   for k,v in sorted(idc.items(),key=lambda kv:-kv[1][0])]
sheet("S15_ID_classes",["Integrated domain class","NLR genes (n)","Species (n)","Species"],d,
 [28,12,11,70],
 "Supplementary Table S15 | Integrated domains collapsed from 86 raw Pfam families into 62 non-redundant classes (e.g. Pkinase and PK_Tyr_Ser-Thr merged).")

# ---- S16 annotation quality ----
tr=json.load(open(f"{ROOT}/genome/truncation_test.json"))
te=json.load(open(f"{ROOT}/genome/te_enrichment10.json"))
sy=json.load(open(f"{ROOT}/genome/synteny_ec_si.json"))
d=[["Annotated CN genes with LRR within 10 kb downstream",f'{tr["hits"].get("CN",0)}/{tr["tot"]["CN"]}',
    f'{100*tr["hits"].get("CN",0)/tr["tot"]["CN"]:.1f}%',f'OR {tr["cn_vs_cnl"][0]:.1f} vs CNL',f'{tr["cn_vs_cnl"][1]:.1e}'],
   ["Annotated N genes with LRR within 10 kb downstream",f'{tr["hits"].get("N",0)}/{tr["tot"]["N"]}',
    f'{100*tr["hits"].get("N",0)/tr["tot"]["N"]:.1f}%',f'OR {tr["n_vs_cnl"][0]:.1f} vs CNL',f'{tr["n_vs_cnl"][1]:.1e}'],
   ["Annotated CNL genes with LRR within 10 kb downstream",f'{tr["hits"].get("CNL",0)}/{tr["tot"]["CNL"]}',
    f'{100*tr["hits"].get("CNL",0)/tr["tot"]["CNL"]:.1f}%',"reference","—"],
   ["NLR flanks carrying a TE protein domain (10 kb)",f'{te["nlr_hits"]}/{te["nlr_flanks"]}',
    f'{100*te["nlr_hits"]/te["nlr_flanks"]:.1f}%',f'OR {te["odds"]:.2f}',f'{te["p"]:.1e}'],
   ["Control flanks carrying a TE protein domain (10 kb)",f'{te["ctrl_hits"]}/{te["ctrl_flanks"]}',
    f'{100*te["ctrl_hits"]/te["ctrl_flanks"]:.1f}%',"reference","—"],
   ["NLRs inside a collinear block (vs S. italica)",f'{sy["nlr_in_block"]}/{sy["nlr_anchored"]}',
    f'{100*sy["nlr_in_block"]/sy["nlr_anchored"]:.1f}%',f'OR {sy["odds"]:.2f}',f'{sy["p"]:.1e}'],
   ["Other genes inside a collinear block",f'{sy["other_in"]}/{sy["other"]}',
    f'{100*sy["other_in"]/sy["other"]:.1f}%',"reference","—"]]
sheet("S16_Genomic_context",["Test","Count","Percentage","Effect size","P value"],d,
 [52,14,13,22,12],
 "Supplementary Table S16 | Genomic context of NLR loci in Eleusine coracana and a test of annotation completeness. Two-thirds of genes annotated as CN have LRR-coding sequence immediately downstream, so architecture classes derived from automated proteome annotation under-call the LRR.")



# ---- S17 NLR-Annotator re-annotation ----
cmp=json.load(open(f"{ROOT}/genome/nlra_compare.json"))
fmrows={r["gene"]:r for r in rows if r["species"]=="Ecoracana"}
d=[]
for g,ann,na,nm,n in sorted(cmp,key=lambda x:(x[1],x[0])):
    r=fmrows.get(g)
    d.append([g,str(r["chrom_lab"]) if r else "-",ann,na if na else "no overlapping locus",
              nm if na else "-",n,
              "LRR gained" if (na and "L" not in ann.rstrip("X") and "L" in na) else
              ("confirmed" if na and ("L" in ann.rstrip("X"))==("L" in na) else
               ("LRR lost" if na else "—"))])
sheet("S17_NLRAnnotator",
 ["E. coracana gene","Chromosome","Architecture (proteome annotation)",
  "Architecture (NLR-Annotator)","NLR motifs","NLR-Annotator loci in the 40 kb region","Outcome"],d,
 [24,13,24,24,12,20,16],
 "Supplementary Table S17 | Independent re-annotation of the 294 annotated finger millet NLR loci with NLR-Annotator v2.1b (gene span plus 20 kb each side). 87% of models lacking an LRR gain one, confirming that CN and N classes are largely truncated gene models.")



# ---- S18 statistical support ----
import statistics as _st
pd_=json.load(open(f"{ROOT}/nlr/robust_pairdist.json"))
kg=json.load(open(f"{ROOT}/nlr/kaks_groups.json"))
d=[["ID-NLRs enriched among genes in head-to-head pairs","26/210 (12.4%) vs 147/3131 (4.7%)","Fisher exact","OR 2.87","2.8e-05","supported"],
   ["ID-NLRs enriched in physical clusters","105/2069 (5.1%) vs 68/1272 (5.3%)","Fisher exact","OR 0.95","0.75","no effect"],
   ["E. coracana head-to-head pairs free of integrated domains","0 of 14 paired genes; 0.19 expected","one-sided binomial","—","0.83","NOT supported"],
   ["E. coracana ID-NLR rate lower than the other genomes","4/294 (1.4%) vs 169/3047 (5.5%)","Fisher exact","OR 0.23","8.0e-04","supported"],
   ["Ka/Ks core+soft-core vs shell orthogroups",f'{_st.median(kg["core"]):.3f} vs {_st.median(kg["shell"]):.3f}',"Mann-Whitney U","—","2.2e-05","supported"],
   ["Ka/Ks shell vs species-specific orthogroups",f'{_st.median(kg["shell"]):.3f} vs {_st.median(kg["spec"]):.3f}',"Mann-Whitney U","—","2.5e-06","supported"],
   ["Ka/Ks ID-NLRs vs NLRs without an ID","0.538 vs 0.405 (n = 85 vs 1,974)","Mann-Whitney U","—","1.5e-10","supported"],
   ["Ka/Ks sensor vs helper within pairs","0.595 vs 0.419 (n = 12 pairs)","Mann-Whitney U","—","0.10","underpowered"],
   ["Ka/Ks sensor vs helper, matched within pair","median difference +0.072 (n = 9)","Wilcoxon signed-rank","—","0.82","underpowered"]]
sheet("S18_Statistical_support",["Claim tested","Observation","Test","Effect size","P value","Outcome"],d,
 [56,40,20,13,12,16],
 "Supplementary Table S18 | Statistical tests of the main claims, including two that the data do not support.")

d=[[f'{int(k)//1000} kb',v["pairs"],v["idpairs"],round(100*v["idpairs"]/max(v["pairs"],1),1),
    round(v["od"],2),f'{v["p"]:.1e}',v["ecor"]] for k,v in sorted(pd_.items(),key=lambda kv:int(kv[0]))]
sheet("S19_Robustness_pairs",
 ["Intergenic distance cut-off","Head-to-head pairs (n)","Pairs with an ID (n)","Pairs with an ID (%)",
  "ID enrichment (odds ratio)","P value","E. coracana pairs (n)"],d,
 [20,18,17,17,18,12,18],
 "Supplementary Table S19 | Head-to-head pair detection and integrated-domain enrichment across intergenic distance cut-offs.")

d=[["40","70",633,13,15,206,399],["50","70",1055,6,10,319,720],["50","80",1271,5,8,324,934],
   ["50","90",1577,4,6,307,1260],["60","80",1744,3,6,320,1415],["70","80",2151,2,5,274,1870]]
sheet("S20_Robustness_orthogroups",
 ["Identity cut-off (%)","Coverage cut-off (%)","Orthogroups (n)","Core","Soft-core","Shell","Species-specific"],d,
 [17,17,16,9,11,9,16],
 "Supplementary Table S20 | Orthogroup composition across sequence identity and coverage cut-offs. The core stays between 2 and 13 orthogroups and species-specific orthogroups stay the majority throughout.")

wb.save(f"{ROOT}/figures/Supplementary_Tables.xlsx")
print("workbook saved; sheets:",len(wb.sheetnames))
