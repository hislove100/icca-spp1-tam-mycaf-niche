import sys, os, re, gc, glob, time
import anndata as ad, pandas as pd, numpy as np
from scipy import sparse
RAW="/sessions/nice-pensive-brown/mnt/outputs/GSE138709/raw"
PER="/sessions/nice-pensive-brown/mnt/outputs/GSE138709/per_sample"
os.makedirs(PER,exist_ok=True)
files=sorted(glob.glob(f"{RAW}/*_UMI.csv.gz"))
i=int(sys.argv[1]); f=files[i]; base=os.path.basename(f)
out=f"{PER}/s{i}.h5ad"
if os.path.exists(out):
    print("already done",base); sys.exit(0)
t=time.time()
m=re.match(r"(GSM\d+)_ICC_(\d+)_(\w+?)_UMI",base)
gsm,patient,tissuefull=m.group(1),m.group(2),m.group(3)
tissue="Tumor" if "Tumor" in tissuefull else "Adjacent"
df=pd.read_csv(f,index_col=0)
genes=df.index.astype(str).to_numpy(); cells=df.columns.astype(str).to_numpy()
X=sparse.csr_matrix(df.values.T.astype(np.float32))
del df; gc.collect()
a=ad.AnnData(X); a.obs_names=cells; a.var_names=genes; a.var_names_make_unique()
a.obs["sample"]=f"ICC_{patient}_{tissuefull}"; a.obs["patient"]=f"ICC_{patient}"
a.obs["tissue"]=tissue; a.obs["gsm"]=gsm
a.write(out)
print(f"OK {base} {a.shape} {time.time()-t:.1f}s")
