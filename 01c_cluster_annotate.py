# GSE138709 single-cell: filtering, normalization, HVG, PCA, Harmony, clustering, UMAP, annotation
# (run after 01b_merge_qc.py which produced the merged raw AnnData)
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, harmonypy
sc.settings.verbosity = 0

a = ad.read_h5ad('adata_raw.h5ad')                 # merged 8-sample raw counts (33,991 cells)

# ---- QC filtering ----
a.var['mt'] = a.var_names.str.upper().str.startswith('MT-')
sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True, percent_top=None)
sc.pp.filter_cells(a, min_genes=200)
a = a[(a.obs.n_genes_by_counts < 6000) & (a.obs.pct_counts_mt < 20)].copy()
sc.pp.filter_genes(a, min_cells=3)                 # -> 32,626 cells x 19,813 genes

# ---- Normalization + HVG ----
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor='seurat', batch_key='patient')

# ---- PCA + Harmony (batch = patient) ----
ah = a[:, a.var.highly_variable].copy()
sc.pp.scale(ah, max_value=10); sc.tl.pca(ah, n_comps=30)
ho = harmonypy.run_harmony(ah.obsm['X_pca'], ah.obs, ['patient'], max_iter_harmony=10)
Zc = np.asarray(ho.Z_corr)
a.obsm['X_pca_harmony'] = (Zc if Zc.shape[0] == ah.n_obs else Zc.T).astype(np.float32)

# ---- Graph, Leiden, UMAP ----
sc.pp.neighbors(a, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=30)
sc.tl.leiden(a, resolution=1.0, flavor='igraph', n_iterations=2, directed=False)   # 25 clusters
sc.tl.umap(a)

# ---- Lineage annotation by background-corrected scoring ----
markers = {
 'Malignant/Epithelial':['EPCAM','KRT19','KRT7','KRT8','KRT18','SOX9','ANXA4','MUC1'],
 'Hepatocyte':['ALB','APOA1','APOA2','TTR','TF','CYP2E1'],
 'T':['CD3D','CD3E','CD2','TRAC','CD8A','CD4','IL7R','CD7'],
 'NK':['NKG7','GNLY','KLRD1','KLRF1','NCAM1'],
 'B':['CD79A','MS4A1','CD19','CD79B','BANK1'],
 'Plasma':['MZB1','IGHG1','JCHAIN','DERL3'],
 'Myeloid/Macrophage':['CD68','CD14','LYZ','C1QA','C1QB','C1QC','AIF1','FCGR3A','ITGAM'],
 'Fibroblast/CAF':['COL1A1','COL1A2','DCN','ACTA2','PDGFRB','FAP','LUM','COL3A1'],
 'Endothelial':['PECAM1','VWF','CDH5','ENG','CLDN5'],
 'Mast':['TPSAB1','TPSB2','CPA3','MS4A2','KIT'],
 'DC':['LILRA4','CLEC9A','CD1C','FCER1A','CLEC10A'],
}
for ct, gs in markers.items():
    sc.tl.score_genes(a, [g for g in gs if g in a.var_names], score_name='sc_'+ct, ctrl_size=50)
scores = a.obs.groupby('leiden')[['sc_'+c for c in markers]].mean()
assign = scores.idxmax(axis=1).str.replace('sc_', '', regex=False)
a.obs['cell_type'] = a.obs['leiden'].map(assign).astype('category')
a.write('annotated.h5ad')
