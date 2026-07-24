# Bulk prognostic scoring and survival (TCGA-CHOL & GSE107943) + 3-cohort meta-analysis
import pandas as pd, numpy as np, gzip
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from scipy.stats import norm

# refined, macrophage-coherent signatures (validated by within-signature correlation)
SIGS = {'CXCL9_IFN_TAM':['CXCL9','CXCL10','CXCL11','GBP1','STAT1'],
        'TREM2_LATAM':['TREM2','GPNMB','C1QC','FCGR3A'],
        'myCAF':['FAP','POSTN','ITGA5','ACTA2','COL1A1','MMP2','DCN','LUM','PDGFRB']}

def score_expr(logexpr):                          # gene x sample (log scale) -> signature scores
    z = logexpr.sub(logexpr.mean(1), axis=0).div(logexpr.std(1)+1e-9, axis=0)
    s = pd.DataFrame({k: z.loc[[g for g in v if g in z.index]].mean() for k,v in SIGS.items()})
    s['niche_LATAM_myCAF'] = (s['TREM2_LATAM'] + s['myCAF'])/2
    return s

def cox(df, col, tcol='os', ecol='ev'):
    d = df[[tcol, ecol, col]].dropna(); d = d[d[tcol] > 0].copy()
    d[col] = (d[col]-d[col].mean())/d[col].std()
    r = CoxPHFitter().fit(d, tcol, ecol).summary.loc[col]
    return r['coef'], r['se(coef)'], int(len(d)), int(d[ecol].sum())

# TCGA-CHOL: STAR TPM (log2 already); map signature genes by Ensembl ID; keep -01 tumours,
# merge with survival (OS_time, OS_status) and histology (intrahepatic subset).
# GSE107943: RPKM matrix (log2(x+1)); survival(mo)+death from the series-matrix characteristics.
# (Full extraction code omitted for brevity; both produce a 'df' with columns os, ev + signature scores.)

def meta_DL(betas, ses):                          # DerSimonian-Laird random-effects
    b, s = np.array(betas), np.array(ses); w = 1/s**2
    fe = (w*b).sum()/w.sum(); Q = (w*(b-fe)**2).sum(); dfr = len(b)-1
    c = w.sum()-(w**2).sum()/w.sum(); tau2 = max(0,(Q-dfr)/c) if c>0 else 0
    wr = 1/(s**2+tau2); re = (wr*b).sum()/wr.sum(); sre = np.sqrt(1/wr.sum())
    I2 = max(0,(Q-dfr)/Q)*100 if Q>0 else 0
    return re, sre, I2
