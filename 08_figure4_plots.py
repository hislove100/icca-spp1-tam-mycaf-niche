#!/usr/bin/env python3
"""
08_figure4_plots.py
Render the cell-cell communication figures from the network matrices (07) and the
SPP1 permutation results (03).

  Figure 4A : SPP1 signalling-pathway network, 22 fine cell states (circle plot)
  Figure 4B : SPP1 ligand-receptor dot-plot (senders = SPP1+ TAM states)
  Figure S5 : global network, lineage resolution
              (A) number of interactions, (B) interaction strength,
              (C) SPP1 signalling-pathway network

Inputs (produced by 07 and 03)
  spp1_finestate_weight.csv, spp1_pathway_weight.csv,
  cellchat_count.csv, cellchat_weight.csv, LR_SPP1_TAM_to_all.csv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Arial if available (see 00_figstyle.py); mathtext regular for superscripts
for f in ("Arial", "Liberation Sans", "DejaVu Sans"):
    try:
        matplotlib.rcParams["font.family"] = f; break
    except Exception:
        continue
matplotlib.rcParams["mathtext.default"] = "regular"
matplotlib.rcParams["svg.fonttype"] = "none"

IN  = os.environ.get("NETWORK_OUT", ".")
OUT = os.environ.get("FIG_OUT", ".")
sup = lambda s: s.replace("+", r"$^{+}$")


# --------------------------------------------------------------------------- #
# circle plot (CellChat netVisual_circle style)                                #
# --------------------------------------------------------------------------- #
def circle(ax, M, labels, title, cmap="tab20", fs=8):
    M = np.asarray(M, float); K = len(labels)
    colors = [matplotlib.colors.to_hex(matplotlib.cm.get_cmap(cmap)(i % 20)) for i in range(K)]
    ang = np.array([np.pi / 2 - 2 * np.pi * i / K for i in range(K)])
    pos = np.c_[np.cos(ang), np.sin(ang)]
    vmax = M.max() if M.max() > 0 else 1.0
    tot = M.sum(0) + M.sum(1); tmax = tot.max() if tot.max() > 0 else 1.0
    for s in range(K):
        for r in range(K):
            v = M[s, r]
            if v <= 0:
                continue
            lw = 0.2 + 3.6 * (v / vmax); al = 0.22 + 0.55 * (v / vmax)
            if s == r:
                ax.add_patch(Circle(pos[s] * 1.09, 0.055, fill=False, lw=lw,
                                    ec=colors[s], alpha=al))
            else:
                rad = 0.2 if (r - s) % K <= K / 2 else -0.2
                ax.add_patch(FancyArrowPatch(pos[s], pos[r],
                             connectionstyle=f"arc3,rad={rad}", arrowstyle="-|>",
                             mutation_scale=7, lw=lw, color=colors[s], alpha=al,
                             shrinkA=9, shrinkB=11, zorder=1))
    for i in range(K):
        ax.scatter(*pos[i], s=120 + 900 * (tot[i] / tmax), color=colors[i],
                   edgecolors="white", linewidths=0.8, zorder=3)
        an = ang[i]
        ha = "left" if np.cos(an) > 0.1 else ("right" if np.cos(an) < -0.1 else "center")
        va = "bottom" if np.sin(an) > 0.1 else ("top" if np.sin(an) < -0.1 else "center")
        ax.text(pos[i, 0] * 1.20, pos[i, 1] * 1.20, sup(labels[i]), fontsize=fs,
                ha=ha, va=va, zorder=4)
    ax.set_xlim(-1.65, 1.65); ax.set_ylim(-1.6, 1.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=12)


# --------------------------------------------------------------------------- #
# Figure 4A : SPP1 pathway network, fine cell states                           #
# --------------------------------------------------------------------------- #
fine = pd.read_csv(os.path.join(IN, "spp1_finestate_weight.csv"), index_col=0)
fig, ax = plt.subplots(figsize=(8.6, 8.6))
circle(ax, fine.values, list(fine.index),
       "SPP1 signaling pathway network (fine cell states)")
fig.savefig(os.path.join(OUT, "Figure_4A.svg"), bbox_inches="tight"); plt.close(fig)

# --------------------------------------------------------------------------- #
# Figure S5 : global network (lineage)                                         #
# --------------------------------------------------------------------------- #
cc = pd.read_csv(os.path.join(IN, "cellchat_count.csv"), index_col=0)
cw = pd.read_csv(os.path.join(IN, "cellchat_weight.csv"), index_col=0)
sp = pd.read_csv(os.path.join(IN, "spp1_pathway_weight.csv"), index_col=0)
fig, axes = plt.subplots(1, 3, figsize=(19, 6.6))
circle(axes[0], cc.values, list(cc.index), "Number of interactions", fs=9)
circle(axes[1], cw.values, list(cw.index), "Interaction weights / strength", fs=9)
circle(axes[2], sp.values, list(sp.index), "SPP1 signaling pathway network", fs=9)
fig.savefig(os.path.join(OUT, "Figure_S5.svg"), bbox_inches="tight"); plt.close(fig)

# --------------------------------------------------------------------------- #
# Figure 4B : SPP1 ligand-receptor dot-plot                                    #
# --------------------------------------------------------------------------- #
lr = pd.read_csv(os.path.join(IN, "LR_SPP1_TAM_to_all.csv"))
pairs = ["SPP1->ITGA5_ITGB1", "SPP1->CD44", "SPP1->ITGAV_ITGB1",
         "SPP1->ITGAV_ITGB5", "SPP1->ITGB1"]
recv = ["myCAF", "Matrix-CAF", "iCAF", "vSMC", "Pericyte", "Endothelial",
        "Malignant/Epithelial", "T", "NK"]
senders = ["SPP1+TREM2+ TAM", "SPP1+CXCL9+ TAM"]
sub = lr[lr.pair.isin(pairs) & lr.receiver.isin(recv)].copy()
sub["z"] = (sub.score - sub.score.mean()) / sub.score.std()
cmap = matplotlib.cm.get_cmap("RdBu_r"); norm = Normalize(-2, 2)
lp = lambda p: -np.log10(p + 1e-3)
sz = lambda l: 10 + l * 80.0
prlabel = [p.replace("->", "_") for p in pairs]
FS = 9
fig = plt.figure(figsize=(9.6, 4.8))
p1 = fig.add_axes([0.30, 0.30, 0.27, 0.60]); p2 = fig.add_axes([0.59, 0.30, 0.27, 0.60])
axL = fig.add_axes([0.90, 0.55, 0.09, 0.34]); axL.axis("off")
axC = fig.add_axes([0.915, 0.13, 0.020, 0.32])


def panel(ax, snd, left):
    d = sub[sub.sender == snd]
    pz = d.pivot(index="pair", columns="receiver", values="z").reindex(index=pairs, columns=recv)
    pp = d.pivot(index="pair", columns="receiver", values="pval").reindex(index=pairs, columns=recv)
    for yi, pr in enumerate(pairs):
        for xi, rc in enumerate(recv):
            z = pz.loc[pr, rc]; pv = pp.loc[pr, rc]
            if pd.isna(pv):
                continue
            l = lp(pv)
            if pv < 0.05:
                ax.scatter(xi, yi, s=sz(l), facecolor=cmap(norm(z)), edgecolor="0.35", linewidths=0.4, zorder=3)
            else:
                ax.scatter(xi, yi, s=max(sz(l), 10), facecolor="none", edgecolor="0.7", linewidths=0.5, zorder=2)
    ax.set_xticks(range(len(recv))); ax.set_xticklabels(recv, rotation=90, fontsize=FS)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels(prlabel if left else [], fontsize=FS, style="italic")
    ax.set_ylim(-0.6, len(pairs) - 0.4); ax.set_xlim(-0.6, len(recv) - 0.4); ax.invert_yaxis()
    for s in ax.spines.values():
        s.set_color("0.8")
    ax.tick_params(length=0)
    ax.set_title("Sender: " + sup(snd), fontsize=FS, pad=4)


panel(p1, senders[0], True); panel(p2, senders[1], False)
axL.set_title("-log10(p-value)", fontsize=FS, loc="left", pad=2)
for l, y in zip([0, 1, 2, 3], np.linspace(0.8, 0.15, 4)):
    axL.scatter(0.28, y, s=sz(l), facecolor="0.6", edgecolor="0.35", linewidths=0.4, transform=axL.transAxes)
    axL.text(0.6, y, str(l), va="center", fontsize=FS, transform=axL.transAxes)
axL.set_xlim(0, 1); axL.set_ylim(0, 1)
sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, cax=axC, ticks=[-2, -1, 0, 1, 2]); cb.ax.tick_params(labelsize=FS)
cb.set_label("Scaled interaction score", fontsize=FS)
fig.suptitle("SPP1 ligand-receptor interaction (dot size = -log10 p; * significant filled)", fontsize=FS, x=0.5, y=0.99)
fig.savefig(os.path.join(OUT, "Figure_4B.svg"), bbox_inches="tight"); plt.close(fig)
print("Saved Figure_4A.svg, Figure_4B.svg, Figure_S5.svg")
