import pandas as pd, numpy as np, glob, re
from scipy.spatial import cKDTree
rng=np.random.default_rng(0)
def hex_neighbors(df):
    # Visium hex: neighbors within array distance; use pixel KDTree (6 nearest ~ adjacency)
    xy=df[['px','py']].values.astype(float)
    tree=cKDTree(xy); d,idx=tree.query(xy,k=7)  # self+6
    pairs=[]
    for i in range(len(df)):
        for j in idx[i,1:]:
            pairs.append((i,j))
    return np.array(pairs)
def bivar_moran(x,y,pairs):
    zx=(x-x.mean())/(x.std()+1e-9); zy=(y-y.mean())/(y.std()+1e-9)
    w=len(pairs)
    I=np.sum(zx[pairs[:,0]]*zy[pairs[:,1]])/w
    # permutation: shuffle y
    null=np.empty(300)
    for k in range(300):
        p=rng.permutation(len(y)); zyp=zy[p]
        null[k]=np.sum(zx[pairs[:,0]]*zyp[pairs[:,1]])/w
    p_emp=(np.sum(null>=I)+1)/(len(null)+1)
    return I,p_emp
rows=[]
for f in sorted(glob.glob('spots_*.csv')):
    s=re.match(r'spots_(.+)\.csv',f.split('/')[-1]).group(1)
    df=pd.read_csv(f)
    if len(df)<80: continue
    cond='pre' if s.endswith('pre') else 'post'
    pairs=hex_neighbors(df)
    for tam in ['SPP1_TAM','TREM2_LATAM','CXCL9_TAM']:
        I,p=bivar_moran(df[tam].values,df['myCAF'].values,pairs)
        rows.append((s,cond,len(df),tam,round(I,3),round(p,3)))
res=pd.DataFrame(rows,columns=['sample','cond','n_spots','TAM_sig','bivar_moranI','perm_p'])
res.to_csv('spatial_neighborhood_moran.csv',index=False)
# summary by condition & signature (n-weighted mean I, fraction significant)
print('=== per-sample bivariate Moran I (TAM x neighbor myCAF) ===')
print(res.to_string(index=False))
print('\n=== summary: mean I by signature x condition ===')
for tam in ['SPP1_TAM','TREM2_LATAM','CXCL9_TAM']:
    for cond in ['pre','post']:
        sub=res[(res.TAM_sig==tam)&(res.cond==cond)]
        wI=np.average(sub.bivar_moranI,weights=sub.n_spots)
        nsig=(sub.perm_p<0.05).sum()
        print(f'{tam:12s} {cond:4s}: weighted I={wI:.3f}  significant {nsig}/{len(sub)} samples')
