import scanpy as sc, anndata as ad, pandas as pd, numpy as np, glob, os, re, gzip
from scipy.io import mmread
from scipy import sparse
from scipy.stats import pearsonr
sc.settings.verbosity=0
RAW='raw'
samples=sorted(set(re.match(r'(GSM\d+_sample_\d+_(?:pre|post))',os.path.basename(f)).group(1)
          for f in glob.glob(f'{RAW}/*_matrix.mtx.gz')))
sigs={
 'SPP1_TAM':['SPP1','TREM2','C1QC','C1QA','C1QB','APOE','CD68'],
 'CXCL9_TAM':['CXCL9','CXCL10','CXCL11','GBP1','STAT1'],
 'TREM2_LATAM':['TREM2','GPNMB','C1QC','FCGR3A','APOC1','CD68'],
 'myCAF':['FAP','POSTN','ITGA5','ACTA2','COL1A1','MMP2','DCN','LUM','PDGFRB'],
}
def load(s):
    feat=pd.read_csv(f'{RAW}/{s}_features.tsv.gz',sep='\t',header=None)
    bc=pd.read_csv(f'{RAW}/{s}_barcodes.tsv.gz',header=None)[0].values
    M=mmread(f'{RAW}/{s}_matrix.mtx.gz').tocsr()   # genes x spots
    a=ad.AnnData(sparse.csr_matrix(M.T.astype(np.float32)))
    a.var_names=feat[1].astype(str).values; a.var_names_make_unique(); a.obs_names=bc
    pos=pd.read_csv(f'{RAW}/{s}_tissue_positions.csv.gz')
    pos.columns=[c.strip('"') for c in pos.columns]
    pos=pos.set_index('barcode')
    pos=pos.reindex(a.obs_names)
    a.obs['in_tissue']=pos['in_tissue'].values
    a.obs['row']=pos['array_row'].values; a.obs['col']=pos['array_col'].values
    a.obs['px']=pos['pxl_row_in_fullres'].values; a.obs['py']=pos['pxl_col_in_fullres'].values
    a=a[a.obs['in_tissue']==1].copy()
    return a

rows=[]; pooled={k:[] for k in ['SPP1_TAM','CXCL9_TAM','TREM2_LATAM']}; pooled_caf=[]; pooled_samp=[]
for s in samples:
    a=load(s)
    if a.n_obs<50: 
        rows.append((s,a.n_obs,np.nan,np.nan,np.nan,'too_few_spots')); continue
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    for k,g in sigs.items():
        sc.tl.score_genes(a,[x for x in g if x in a.var_names],score_name=k,ctrl_size=50)
    caf=a.obs['myCAF'].values
    r_spp1=pearsonr(a.obs['SPP1_TAM'],caf)[0]
    r_cxcl9=pearsonr(a.obs['CXCL9_TAM'],caf)[0]
    r_trem2=pearsonr(a.obs['TREM2_LATAM'],caf)[0]
    rows.append((s,a.n_obs,round(r_spp1,3),round(r_cxcl9,3),round(r_trem2,3),'ok'))
    # pooled (z within sample)
    z=lambda v:(v-v.mean())/(v.std()+1e-9)
    for k in pooled: pooled[k].append(z(a.obs[k].values))
    pooled_caf.append(z(caf)); pooled_samp+=[s]*a.n_obs
    a.obs[['row','col','px','py','SPP1_TAM','CXCL9_TAM','TREM2_LATAM','myCAF']].to_csv(f'spots_{s}.csv')
res=pd.DataFrame(rows,columns=['sample','n_spots','r_SPP1TAM_myCAF','r_CXCL9TAM_myCAF','r_TREM2TAM_myCAF','status'])
res.to_csv('spatial_colocalization_by_sample.csv',index=False)
print(res.to_string(index=False))
# pooled correlations
for k in pooled:
    x=np.concatenate(pooled[k]); y=np.concatenate(pooled_caf)
    r,p=pearsonr(x,y); print(f'POOLED {k} vs myCAF: r={r:.3f} p={p:.2e} (n_spots={len(x)})')
