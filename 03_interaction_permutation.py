# SPP1 ligand-receptor interaction with a CellPhoneDB-style label-permutation test
import scanpy as sc, anndata as ad, pandas as pd, numpy as np
sc.settings.verbosity = 0
a = ad.read_h5ad('annotated_fine.h5ad')      # cell_type with myeloid & CAF subtypes merged into 'label'
genes = ['SPP1','CD44','ITGA5','ITGB1','ITGAV','ITGB5']
E = sc.get.obs_df(a, keys=genes).values.astype(np.float32)
codes, cats = pd.factorize(a.obs['label'].values); K = len(cats)
gi = {g:i for i,g in enumerate(genes)}
def group_means(code):
    s = np.zeros((K, len(genes))); cnt = np.zeros(K)
    np.add.at(s, code, E); np.add.at(cnt, code, 1.0)
    return s / cnt[:,None]
pairs = {'SPP1->CD44':['CD44'], 'SPP1->ITGA5_ITGB1':['ITGA5','ITGB1'],
         'SPP1->ITGAV_ITGB1':['ITGAV','ITGB1'], 'SPP1->ITGAV_ITGB5':['ITGAV','ITGB5']}
def score_matrix(mean):                        # sender x receiver; complex = min subunit
    L = mean[:, gi['SPP1']]; out = {}
    for name, rec in pairs.items():
        R = np.min(np.stack([mean[:, gi[r]] for r in rec], 1), axis=1)
        out[name] = (L[:,None] + R[None,:]) / 2.0
    return out
obs = score_matrix(group_means(codes))
NP = 1000; rng = np.random.default_rng(0)
pval = {k: np.zeros_like(v) for k,v in obs.items()}
for _ in range(NP):
    ps = score_matrix(group_means(rng.permutation(codes)))
    for k in obs: pval[k] += (ps[k] >= obs[k])
for k in pval: pval[k] /= NP                    # empirical p-value
