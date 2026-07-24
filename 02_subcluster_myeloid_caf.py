# Independent re-clustering of the myeloid and stromal compartments + subtype annotation
import scanpy as sc, anndata as ad, numpy as np, harmonypy
sc.settings.verbosity = 0
a = ad.read_h5ad('annotated.h5ad')

# ---------- Myeloid (4,660 cells) ----------
mye = a[a.obs['cell_type'] == 'Myeloid/Macrophage'].copy()
sc.pp.highly_variable_genes(mye, n_top_genes=1500, flavor='seurat', batch_key='patient')
mh = mye[:, mye.var.highly_variable].copy()
sc.pp.scale(mh, max_value=10); sc.tl.pca(mh, n_comps=20)
ho = harmonypy.run_harmony(mh.obsm['X_pca'], mh.obs, ['patient'], max_iter_harmony=10)
Zc = np.asarray(ho.Z_corr); mye.obsm['X_harmony'] = (Zc if Zc.shape[0]==mh.n_obs else Zc.T).astype(np.float32)
sc.pp.neighbors(mye, use_rep='X_harmony', n_neighbors=15, n_pcs=20)
sc.tl.leiden(mye, resolution=0.6, flavor='igraph', n_iterations=2, directed=False, key_added='mye_clust')
sc.tl.umap(mye)
mye_labels = {'0':'Mono-ISG15','1':'SPP1+CXCL9+ TAM','2':'Cycling Mac','3':'SPP1+TREM2+ TAM',
 '4':'cDC2 (CD1C+)','5':'FCN1+ Monocyte','6':'FCN1+ Mono-Mac','7':'Transitional Mac','8':'C1QC+MARCO+ Kupffer-like'}
mye.obs['mye_subtype'] = mye.obs['mye_clust'].map(mye_labels).astype('category')
mye.obs['is_SPP1_TAM'] = mye.obs['mye_subtype'].isin(['SPP1+CXCL9+ TAM','SPP1+TREM2+ TAM'])
mye.write('myeloid.h5ad')

# ---------- Stroma / CAF (521 cells) ----------
caf = a[a.obs['cell_type'] == 'Fibroblast/CAF'].copy()
sc.pp.highly_variable_genes(caf, n_top_genes=1000, flavor='seurat')
ch = caf[:, caf.var.highly_variable].copy()
sc.pp.scale(ch, max_value=10); sc.tl.pca(ch, n_comps=15)
sc.pp.neighbors(ch, use_rep='X_pca', n_neighbors=12, n_pcs=15)
sc.tl.leiden(ch, resolution=0.5, flavor='igraph', n_iterations=2, directed=False, key_added='caf_clust')
sc.tl.umap(ch)
caf.obs['caf_clust'] = ch.obs['caf_clust'].values; caf.obsm['X_umap'] = ch.obsm['X_umap']
# consolidate 11 fine clusters into biological subtypes by marker means (see manuscript)
grp = {'0':'Pericyte','2':'Pericyte','4':'Pericyte','6':'Pericyte','1':'vSMC',
       '3':'myCAF','5':'myCAF','9':'myCAF','7':'iCAF','8':'iCAF','10':'Matrix-CAF'}
caf.obs['caf_subtype'] = caf.obs['caf_clust'].map(grp).astype('category')
caf.write('caf.h5ad')
