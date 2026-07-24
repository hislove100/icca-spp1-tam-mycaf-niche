#!/usr/bin/env python3
"""
E-MTAB-6389 iCCA - immunosuppressive SPP1/TREM2 TAM-myCAF axis prognostic validation
+ 3-cohort random-effects meta-analysis (GSE107943 + E-MTAB-6389 + TCGA-iCCA).

Fixes the real data issue: data_exp_icc.txt is indexed by Affymetrix HTA-2.0 EXON-level
PSR probe IDs (e.g. PSR01001649.hg.1), not gene symbols. We map PSR -> gene using
psr2genes.json (built from the HTA-2.0 na36 probeset annotation via transcript_cluster_id),
collapse probes to gene level (mean of log2 intensities), then run the identical scoring /
Cox / KM / meta pipeline.
"""
import os, sys, csv, re, json, argparse
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--folder", default=os.path.dirname(os.path.abspath(__file__)))
args = ap.parse_args()
HERE = args.folder
EXPR  = os.path.join(HERE, "data_exp_icc.txt")
SDRF  = os.path.join(HERE, "E-MTAB-6389.sdrf.txt")
PRIOR = os.path.join(HERE, "prior_cohorts.csv")
MAP   = os.path.join(HERE, "psr2genes.json")

SIGS = {
    "CXCL9_IFN_TAM": ["CXCL9", "CXCL10", "CXCL11", "GBP1", "STAT1"],
    "TREM2_LATAM":   ["TREM2", "GPNMB", "C1QC", "FCGR3A"],
    "myCAF":         ["FAP", "POSTN", "ITGA5", "ACTA2", "COL1A1", "MMP2", "DCN", "LUM", "PDGFRB"],
    "SPP1_TAM":      ["SPP1", "TREM2", "C1QC", "C1QA", "C1QB", "APOE", "CD68"],
}
def die(m): print("\n[STOP] " + m); sys.exit(1)
for p in (EXPR, SDRF, PRIOR, MAP):
    if not os.path.exists(p): die(f"missing required file: {p}")

psr2genes = json.load(open(MAP))
need = set(psr2genes.keys())
genes_all = sorted({g for v in psr2genes.values() for g in v})
print(f"[map] {len(need)} PSR probes -> {len(genes_all)} genes")

with open(EXPR) as f:
    header = f.readline().rstrip("\n").split("\t")
samples = [h.strip().strip('"') for h in header]
print(f"[expr] samples in matrix = {len(samples)}")

gene_sum = {g: np.zeros(len(samples)) for g in genes_all}
gene_cnt = {g: np.zeros(len(samples)) for g in genes_all}
grabbed = 0
with open(EXPR) as f:
    f.readline()
    for line in f:
        tab = line.find("\t")
        pid = line[:tab].strip().strip('"')
        if pid not in need: continue
        vals = line.rstrip("\n").split("\t")[1:]
        arr = np.array([float(x) if x not in ("", "NA", "NaN") else np.nan for x in vals], float)
        for g in psr2genes[pid]:
            m = ~np.isnan(arr); gene_sum[g][m] += arr[m]; gene_cnt[g][m] += 1
        grabbed += 1
print(f"[expr] signature probe rows captured = {grabbed}/{len(need)}")

genemat = {g: gene_sum[g] / np.where(gene_cnt[g] == 0, np.nan, gene_cnt[g]) for g in genes_all}
expr = pd.DataFrame(genemat, index=samples).T
if np.nanmax(expr.values) > 50:
    expr = np.log2(expr + 1); print("[expr] applied log2(x+1)")
print(f"[expr] gene-level matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

sd = pd.read_csv(SDRF, sep="\t", dtype=str)
def findcol(pat):
    for c in sd.columns:
        if re.search(pat, c, re.I): return c
    return None
os_col = findcol(r"overall survival"); ev_col = findcol(r"event death|vital|death")
src = findcol(r"^source name") or "Source Name"
print(f"[sdrf] os={os_col!r} event={ev_col!r} name={src!r}")
sd["_k"] = sd[src].str.strip()
clin = pd.DataFrame({"os": pd.to_numeric(sd[os_col], errors="coerce").values,
                     "ev": pd.to_numeric(sd[ev_col], errors="coerce").values}, index=sd["_k"].values).dropna()
clin = clin[~clin.index.duplicated()]

tum = [s for s in expr.columns if s.endswith("G") and not s.endswith("NG") and s in clin.index]
X = expr[tum]
z = X.sub(X.mean(1), axis=0).div(X.std(1) + 1e-9, axis=0)
score = pd.DataFrame({k: z.loc[[g for g in v if g in z.index]].mean() for k, v in SIGS.items()})
score["niche_LATAM_myCAF"] = (score["TREM2_LATAM"] + score["myCAF"]) / 2
df = score.copy()
df["os"] = clin.loc[df.index, "os"].values
df["ev"] = clin.loc[df.index, "ev"].values
df = df.dropna(subset=["os", "ev"]); df = df[df["os"] > 0]
df.to_csv(os.path.join(HERE, "emtab6389_scored.csv"))
print(f"[surv] analyzable tumors = {len(df)}, events = {int(df.ev.sum())}")

from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
out = []
for col in ["CXCL9_IFN_TAM", "TREM2_LATAM", "myCAF", "niche_LATAM_myCAF", "SPP1_TAM"]:
    d = df[["os", "ev", col]].copy(); d[col] = (d[col] - d[col].mean()) / d[col].std()
    c = CoxPHFitter().fit(d, "os", "ev"); r = c.summary.loc[col]
    hi = df[col] > df[col].median()
    lr = logrank_test(df.os[hi], df.os[~hi], df.ev[hi], df.ev[~hi])
    out.append([col, r["coef"], r["se(coef)"], np.exp(r["coef"]),
                np.exp(r["coef"]-1.96*r["se(coef)"]), np.exp(r["coef"]+1.96*r["se(coef)"]),
                r["p"], lr.p_value, len(df), int(df.ev.sum())])
res = pd.DataFrame(out, columns=["sig","logHR","se","HR","CI_low","CI_high","cox_p","logrank_p","n","events"])
res.to_csv(os.path.join(HERE, "emtab6389_cox.csv"), index=False)
print("\n=== E-MTAB-6389 Cox (HR per SD) ===")
print(res[["sig","HR","CI_low","CI_high","cox_p","logrank_p"]].round(3).to_string(index=False))

prior = pd.read_csv(PRIOR)
def meta(bs, ses):
    bs, ses = np.array(bs), np.array(ses); w = 1/ses**2
    fe = np.sum(w*bs)/np.sum(w); Q = np.sum(w*(bs-fe)**2); dfr = len(bs)-1
    cc = np.sum(w)-np.sum(w**2)/np.sum(w); tau2 = max(0,(Q-dfr)/cc) if cc>0 else 0
    wr = 1/(ses**2+tau2); re = np.sum(wr*bs)/np.sum(wr); sre = np.sqrt(1/np.sum(wr))
    I2 = max(0,(Q-dfr)/Q)*100 if Q>0 else 0
    return re, sre, I2
from scipy.stats import norm as NORM
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
meta_rows = []
fig, axes = plt.subplots(2, 1, figsize=(8, 5.4))
for ax, sig, title in [(axes[0], "TREM2_LATAM", "TREM2+ LA-TAM signature"),
                       (axes[1], "niche_LATAM_myCAF", "Immunosuppressive niche (TREM2+ LA-TAM x myCAF)")]:
    g107 = prior[(prior.sig==sig) & (prior.cohort=="GSE107943(iCCA)")].iloc[0]
    tica = prior[(prior.sig==sig) & (prior.cohort=="TCGA-iCCA")].iloc[0]
    em = res[res.sig==sig].iloc[0]
    entries = [("GSE107943 iCCA (n=30)", g107.logHR, g107.se),
               ("E-MTAB-6389 iCCA (n=%d)" % em.n, em.logHR, em.se),
               ("TCGA-iCCA (n=27)", tica.logHR, tica.se)]
    b = [e[1] for e in entries]; s = [e[2] for e in entries]
    re, sre, I2 = meta(b, s)
    entries.append(("Meta (3 iCCA cohorts, RE)", re, sre))
    for nm, b_, s_ in entries:
        meta_rows.append([sig, nm, b_, s_, np.exp(b_), np.exp(b_-1.96*s_), np.exp(b_+1.96*s_),
                          2*(1-NORM.cdf(abs(b_/s_))), I2 if "Meta" in nm else np.nan])
    y = 0; yt=[]; yl=[]
    for nm, b_, s_ in entries:
        ismeta = "Meta" in nm; col = "#c0392b" if ismeta else "#2c3e50"
        hr, lo, hi = np.exp(b_), np.exp(b_-1.96*s_), np.exp(b_+1.96*s_)
        ax.plot([lo,hi],[y,y], color=col, lw=2.5 if ismeta else 1.8)
        ax.plot(hr, y, "D" if ismeta else "o", color=col, ms=9 if ismeta else 7)
        ax.text(5.2, y, f"{hr:.2f} ({lo:.2f}-{hi:.2f})", va="center", fontsize=8.5)
        yt.append(y); yl.append(nm); y -= 1
    ax.axvline(1, ls="--", color="grey", lw=1); ax.set_xscale("log"); ax.set_xlim(0.4, 9)
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=9); ax.set_ylim(min(yt)-0.6, max(yt)+0.6)
    ax.set_xticks([0.5,1,2,4]); ax.set_xticklabels(["0.5","1","2","4"])
    ax.set_xlabel("Hazard ratio per SD (95% CI)", fontsize=9)
    ax.set_title(f"{title}   (I2={I2:.0f}%)", fontsize=10)
fig.suptitle("3-cohort iCCA prognostic meta-analysis (higher HR = worse survival)", fontsize=11)
plt.tight_layout(); plt.savefig(os.path.join(HERE, "meta3_forest.png"), dpi=140, bbox_inches="tight")
pd.DataFrame(meta_rows, columns=["sig","cohort","logHR","se","HR","CI_low","CI_high","p","I2"]).to_csv(
    os.path.join(HERE, "meta3_results.csv"), index=False)
print("\n=== meta written ===\nDone.")
