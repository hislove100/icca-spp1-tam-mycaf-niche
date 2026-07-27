#!/usr/bin/env python3
"""
07_cellchat_network.py
Aggregate cell-cell communication network + SPP1 signalling-pathway network.

Produces the matrices behind:
  * Figure 4A  - SPP1 signalling-pathway network at fine cell-state resolution
  * Figure S5  - global intercellular communication network (interaction number /
                 interaction strength / SPP1 pathway) at lineage resolution

Inputs
  annotated_fine.h5ad   (from 01c/02): log1p-normalised AnnData with
                        obs['cell_type'] = 10 lineages, obs['label'] = 22 fine states.

Outputs (CSV, sender x receiver)
  cellchat_count.csv          lineage x lineage : number of expressed L-R pairs
  cellchat_weight.csv         lineage x lineage : summed interaction strength
  spp1_pathway_weight.csv     lineage x lineage : SPP1-pathway strength
  spp1_finestate_weight.csv   22 x 22 fine state: SPP1-pathway strength

Method
  * The SPP1 pathway score reuses the ligand-receptor scoring of
    03_interaction_permutation.py: for each SPP1 receptor pair the score is the
    mean of the (log-normalised) ligand and receptor group means, complexes taken
    as the minimum subunit; the pathway weight is the sum over SPP1 receptor pairs.
  * The global network is a transparent CellChat-style aggregation over the LIANA
    'consensus' ligand-receptor resource, using a Hill function on linear group-mean
    expression:  P = (L*R) / (Kh + L*R), Kh = 0.5  (receptor complex = geometric
    mean of subunit means). 'count' = number of pairs expressed in both partners;
    'weight' = sum of P.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sc.settings.verbosity = 0
HERE = os.path.dirname(os.path.abspath(__file__))
H5   = os.environ.get("ANNOTATED_FINE", "annotated_fine.h5ad")
OUT  = os.environ.get("NETWORK_OUT", ".")
KH   = 0.5   # Hill half-saturation (CellChat default)

a = ad.read_h5ad(H5)

# ----------------------------------------------------------------------------- #
# 1. SPP1 pathway network (same L-R scoring as 03_interaction_permutation.py)     #
# ----------------------------------------------------------------------------- #
SPP1_PAIRS = {
    "SPP1->CD44":        ["CD44"],
    "SPP1->ITGA5_ITGB1": ["ITGA5", "ITGB1"],
    "SPP1->ITGAV_ITGB1": ["ITGAV", "ITGB1"],
    "SPP1->ITGAV_ITGB5": ["ITGAV", "ITGB5"],
}
SPP1_GENES = ["SPP1", "CD44", "ITGA5", "ITGB1", "ITGAV", "ITGB5"]


def spp1_network(adata, group_key):
    """sender x receiver SPP1-pathway weight for the given grouping."""
    E = sc.get.obs_df(adata, keys=SPP1_GENES).values.astype(np.float64)
    cats = pd.Categorical(adata.obs[group_key].astype(str))
    codes = cats.codes
    labels = list(cats.categories)
    K = len(labels)
    gi = {g: i for i, g in enumerate(SPP1_GENES)}
    s = np.zeros((K, len(SPP1_GENES))); n = np.zeros(K)
    np.add.at(s, codes, E); np.add.at(n, codes, 1.0)
    mean = s / n[:, None]                       # group mean (log-normalised)
    L = mean[:, gi["SPP1"]]
    W = np.zeros((K, K))
    for rec in SPP1_PAIRS.values():
        R = np.min(np.stack([mean[:, gi[r]] for r in rec], 1), axis=1)
        W += (L[:, None] + R[None, :]) / 2.0    # sender x receiver
    return pd.DataFrame(W, index=labels, columns=labels)


spp1_lineage = spp1_network(a, "cell_type")
spp1_fine    = spp1_network(a, "label")
spp1_lineage.to_csv(os.path.join(OUT, "spp1_pathway_weight.csv"))
spp1_fine.to_csv(os.path.join(OUT, "spp1_finestate_weight.csv"))
print(f"[spp1] lineage {spp1_lineage.shape}, fine {spp1_fine.shape}")

# ----------------------------------------------------------------------------- #
# 2. Global CellChat-style aggregate network over the LIANA consensus resource    #
# ----------------------------------------------------------------------------- #
try:
    import liana as li
    res = li.rs.select_resource("consensus")        # ligand / receptor columns
except Exception as e:                              # pragma: no cover
    raise SystemExit(f"[stop] liana consensus resource unavailable: {e}")

lin = pd.Categorical(a.obs["cell_type"].astype(str))
lin_labels = list(lin.categories)
K = len(lin_labels)

var_names = set(a.var_names)
needed = sorted({g for col in ("ligand", "receptor") for pair in res[col]
                 for g in str(pair).split("_") if g in var_names})
Elin = np.expm1(sc.get.obs_df(a, keys=needed).values.astype(np.float64))
cats = lin.codes
gmean = np.zeros((K, len(needed))); n = np.zeros(K)
np.add.at(gmean, cats, Elin); np.add.at(n, cats, 1.0)
gmean = gmean / n[:, None]
col = {g: i for i, g in enumerate(needed)}


def subunit_mean(complex_str):
    subs = [s for s in str(complex_str).split("_") if s in col]
    if not subs:
        return None
    vals = np.stack([gmean[:, col[s]] for s in subs], 1)
    return np.exp(np.log(vals + 1e-12).mean(1))     # geometric mean of subunits


count = np.zeros((K, K)); weight = np.zeros((K, K))
for _, row in res.iterrows():
    L = subunit_mean(row["ligand"]); R = subunit_mean(row["receptor"])
    if L is None or R is None:
        continue
    LR = L[:, None] * R[None, :]
    P = LR / (KH + LR)
    expressed = (L[:, None] > 0) & (R[None, :] > 0)
    count += expressed.astype(float)
    weight += P
pd.DataFrame(count, index=lin_labels, columns=lin_labels).to_csv(os.path.join(OUT, "cellchat_count.csv"))
pd.DataFrame(weight, index=lin_labels, columns=lin_labels).to_csv(os.path.join(OUT, "cellchat_weight.csv"))
print(f"[global] count/weight ({K}x{K}) over {len(res)} L-R pairs")
print("Done.")
