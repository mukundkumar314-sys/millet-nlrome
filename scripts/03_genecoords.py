import gzip, re, os, json
import os
ROOT = os.environ.get("PROJECT_ROOT", "/home/claude")
GFF=f"{ROOT}/gff"
SP={"Sitalica":"Sitalica.gff3.gz","Sviridis":"Sviridis.gff3.gz","Phallii":"Phallii.gff3.gz",
    "Dexilis":"Dexilis.gff3.gz","Ecrusgalli":"Ecrusgalli.gff3.gz","Sbicolor":"Sbicolor.gff3.gz",
    "Osativa":"Osativa.gff3.gz","Zmays":"Zmays.gff3.gz","Bdistachyon":"Bdistachyon.gff3.gz",
    "Ecoracana":"Ecoracana.gff.gz"}
out={}
for sp,fn in SP.items():
    genes=[]
    with gzip.open(f"{GFF}/{fn}",'rt') as fh:
        for ln in fh:
            if ln.startswith('#'): continue
            f=ln.rstrip('\n').split('\t')
            if len(f)<9 or f[2]!='gene': continue
            a=f[8]
            if sp=="Ecoracana":
                m=re.search(r'locus_tag=([^;]+)',a) or re.search(r'ID=gene-([^;]+)',a)
            else:
                m=re.search(r'ID=gene:([^;]+)',a) or re.search(r'ID=([^;]+)',a)
            if not m: continue
            genes.append((m.group(1), f[0], int(f[3]), int(f[4])))
    out[sp]=genes
    chroms=len({g[1] for g in genes})
    print(f"{sp:<13} genes={len(genes):>7}  seqs={chroms}")
json.dump(out,open(f"{ROOT}/nlr/all_genes.json","w"))
print("saved all_genes.json")
