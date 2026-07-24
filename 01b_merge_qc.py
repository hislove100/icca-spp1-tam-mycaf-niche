import scanpy as sc, anndata as ad, pandas as pd, numpy as np, glob, gc
from scipy import sparse
sc.settings.verbosity=0
OUT="/sessions/nice-pensive-brown/mnt/outputs/GSE138709"
parts=[ad.read_h5ad(f) for f in sorted(glob.glob(f"{OUT}/per_sample/s*.h5ad"))]
adata=ad.concat(parts,join="outer",index_unique="-")
del parts; gc.collect()
adata.X=sparse.csr_matrix(adata.X)
print("MERGED",adata.shape)

adata.var["mt"]=adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"]=adata.var_names.str.upper().str.startswith(("RPS","RPL"))
sc.pp.calculate_qc_metrics(adata,qc_vars=["mt","ribo"],inplace=True,percent_top=None)
print("mito genes:",int(adata.var['mt'].sum()))

qc=adata.obs.groupby("sample").agg(
    n_cells=("n_genes_by_counts","size"),
    med_genes=("n_genes_by_counts","median"),
    med_counts=("total_counts","median"),
    med_pct_mt=("pct_counts_mt","median")).round(1)
tis=adata.obs.groupby("sample")["tissue"].first()
qc.insert(0,"tissue",tis)
qc.to_csv(f"{OUT}/qc_by_sample_prefilter.csv")
print(qc.to_string())
print("\nTOTAL cells (raw):",adata.n_obs,"| genes:",adata.n_vars)
print("tissue split:\n",adata.obs['tissue'].value_counts().to_string())
adata.write(f"{OUT}/adata_raw.h5ad")
print("saved adata_raw.h5ad")
