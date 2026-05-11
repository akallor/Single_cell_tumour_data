#Visualizations of the scRNA-seq data

#!/usr/bin/env python3
"""
scRNA-seq Reporting Script — Python / matplotlib + scanpy
Plots: UMAP (cluster, condition, marker gene), dot plot, composition bar chart
Dataset: NF1 neurofibromatosis — Sus scrofa (Tumor1, Tumor2, NormalAdjacent)
Dependencies: scanpy, anndata, matplotlib, seaborn, pandas, numpy
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")                         # non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import seaborn as sns
import scanpy as sc

# ── 0. Configuration ──────────────────────────────────────────────────────────

ANNDATA_H5AD  = "integrated_adata.h5ad"      # path to the AnnData object
OUT_DIR       = "figures/"                    # output directory
MARKER_GENE   = "NF1"                         # gene of interest
CONDITION_KEY = "condition"                   # obs column: "tumour" / "normal"
SAMPLE_KEY    = "sample"                      # obs column: sample name
CLUSTER_KEY   = "leiden"                      # obs column: cluster labels
N_TOP_MARKERS = 5                             # top markers per cluster for dot plot
DPI           = 200                           # raster output resolution

CONDITION_PAL = {"normal": "#4393C3", "tumour": "#D6604D"}

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.verbosity = 1

# ── 1. Load data ──────────────────────────────────────────────────────────────

print("Loading AnnData object...")
adata = sc.read_h5ad(ANNDATA_H5AD)

# Validate required columns
for col in [CONDITION_KEY, CLUSTER_KEY, SAMPLE_KEY]:
    assert col in adata.obs.columns, f"Missing obs column: {col}"

# Ensure condition order
adata.obs[CONDITION_KEY] = pd.Categorical(
    adata.obs[CONDITION_KEY], categories=["normal", "tumour"], ordered=True
)

clusters = sorted(adata.obs[CLUSTER_KEY].unique(), key=lambda x: int(x))
n_clusters = len(clusters)

# Cluster colour palette
cluster_colors = dict(zip(
    clusters,
    plt.cm.tab20.colors[:n_clusters] if n_clusters <= 20
    else plt.cm.nipy_spectral(np.linspace(0, 1, n_clusters))
))

# ── 2. Shared plot style ──────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family"      : "DejaVu Sans",
    "axes.spines.top"  : False,
    "axes.spines.right": False,
    "figure.facecolor" : "white",
    "axes.facecolor"   : "#F7F7F7",
    "axes.grid"        : False,
    "savefig.dpi"      : DPI,
    "savefig.bbox"     : "tight",
    "savefig.facecolor": "white",
})

def style_umap_ax(ax, title, subtitle=""):
    """Remove ticks, add labels and titles for UMAP axes."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1", fontsize=9, color="grey")
    ax.set_ylabel("UMAP 2", fontsize=9, color="grey")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
    if subtitle:
        ax.text(0.5, 1.01, subtitle, ha="center", va="bottom",
                transform=ax.transAxes, fontsize=8, color="grey")

# ── 3. Extract UMAP coordinates ───────────────────────────────────────────────

print("Extracting UMAP embedding...")
umap = adata.obsm["X_umap"]
obs  = adata.obs.copy()
obs["UMAP_1"] = umap[:, 0]
obs["UMAP_2"] = umap[:, 1]

# Marker gene expression
if MARKER_GENE in adata.var_names:
    gene_idx = list(adata.var_names).index(MARKER_GENE)
    import scipy.sparse as sp
    raw_expr = adata.X[:, gene_idx]
    obs[MARKER_GENE] = np.asarray(raw_expr.todense()).flatten() \
        if sp.issparse(adata.X) else raw_expr
else:
    warnings.warn(f"{MARKER_GENE} not found in var_names. Expression UMAP will be blank.")
    obs[MARKER_GENE] = np.nan

# Sort so high-expressing cells plot on top
obs_expr = obs.sort_values(MARKER_GENE, na_position="first")

# ── 4. UMAP plots ─────────────────────────────────────────────────────────────

print("Generating UMAP plots...")

fig, axes = plt.subplots(1, 3, figsize=(21, 7))
fig.suptitle(
    "scRNA-seq · Sus scrofa NF1 neurofibromatosis",
    fontsize=14, fontweight="bold", y=1.01
)

# 4a. UMAP by cluster
ax = axes[0]
for cl in clusters:
    mask = obs[CLUSTER_KEY] == cl
    ax.scatter(obs.loc[mask, "UMAP_1"], obs.loc[mask, "UMAP_2"],
               c=[cluster_colors[cl]], s=1.5, alpha=0.6, linewidths=0, label=cl)
legend_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=cluster_colors[cl],
           markersize=6, label=cl)
    for cl in clusters
]
ax.legend(handles=legend_handles, title="Cluster", fontsize=7,
          title_fontsize=8, markerscale=1, ncol=max(1, n_clusters // 12),
          loc="lower right", framealpha=0.7, edgecolor="none")
style_umap_ax(ax, "Cell clusters",
              f"{len(obs):,} cells · Leiden clustering")

# 4b. UMAP by condition
ax = axes[1]
for cond, colour in CONDITION_PAL.items():
    mask = obs[CONDITION_KEY] == cond
    ax.scatter(obs.loc[mask, "UMAP_1"], obs.loc[mask, "UMAP_2"],
               c=colour, s=1.5, alpha=0.5, linewidths=0, label=cond)
ax.legend(title="Condition", fontsize=9, title_fontsize=9,
          markerscale=3, loc="lower right", framealpha=0.7, edgecolor="none")
style_umap_ax(ax, "Condition",
              "Blue = NormalAdjacent · Red = Tumor1 + Tumor2")

# 4c. UMAP by marker gene expression
ax = axes[2]
expr_vals = obs_expr[MARKER_GENE].values
vmax = np.nanpercentile(expr_vals, 99) or 1.0
sc_plot = ax.scatter(
    obs_expr["UMAP_1"], obs_expr["UMAP_2"],
    c=expr_vals, cmap="magma", s=1.5, alpha=0.8,
    linewidths=0, vmin=0, vmax=vmax
)
cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.04, pad=0.02, aspect=30)
cbar.set_label(f"{MARKER_GENE}\n(log-norm)", fontsize=8)
cbar.ax.tick_params(labelsize=7)
style_umap_ax(ax, f"{MARKER_GENE} expression",
              "Normalised log expression · high expressors on top")

plt.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT_DIR, f"umap_combined.{ext}"))
plt.close(fig)
print(f"Saved umap_combined.pdf/.png")

# ── 5. Find top marker genes ──────────────────────────────────────────────────

print("Computing cluster marker genes...")
sc.tl.rank_genes_groups(
    adata,
    groupby   = CLUSTER_KEY,
    method    = "wilcoxon",
    use_raw   = False,
    pts       = True,       # compute percent-expressed
    n_genes   = N_TOP_MARKERS * 3,
)

# Extract to DataFrame
marker_dfs = []
for cl in clusters:
    df = sc.get.rank_genes_groups_df(adata, group=cl, pval_cutoff=0.05)
    df.insert(0, "cluster", cl)
    marker_dfs.append(df.head(N_TOP_MARKERS))

top_markers = pd.concat(marker_dfs, ignore_index=True)
top_markers.to_csv(os.path.join(OUT_DIR, "top_cluster_markers.csv"), index=False)
print("Saved top_cluster_markers.csv")

genes_ordered = top_markers.drop_duplicates("names")["names"].tolist()

# ── 6. Dot plot ───────────────────────────────────────────────────────────────

print("Generating dot plot...")

# Compute per-cluster statistics for each gene
import scipy.sparse as sp

def get_expr_matrix(adata, genes):
    """Return dense expression matrix for selected genes."""
    idx = [list(adata.var_names).index(g) for g in genes if g in adata.var_names]
    valid_genes = [g for g in genes if g in adata.var_names]
    X = adata.X[:, idx]
    if sp.issparse(X):
        X = np.asarray(X.todense())
    return valid_genes, X

valid_genes, expr_mat = get_expr_matrix(adata, genes_ordered)

dot_rows = []
for cl in clusters:
    mask = (adata.obs[CLUSTER_KEY] == cl).values
    sub  = expr_mat[mask, :]
    for gi, gene in enumerate(valid_genes):
        col = sub[:, gi]
        dot_rows.append({
            "cluster"     : cl,
            "gene"        : gene,
            "avg_expr"    : np.mean(np.expm1(col)),
            "pct_express" : np.mean(col > 0) * 100,
        })

dot_df = pd.DataFrame(dot_rows)

# Scale average expression within each gene (0–1)
dot_df["avg_scaled"] = dot_df.groupby("gene")["avg_expr"].transform(
    lambda x: (x - x.min()) / (x.max() - x.min() + 1e-9)
)

dot_df["cluster"] = pd.Categorical(dot_df["cluster"], categories=clusters, ordered=True)
dot_df["gene"]    = pd.Categorical(dot_df["gene"], categories=list(reversed(valid_genes)), ordered=True)

fig_w = max(8, n_clusters * 0.85)
fig_h = max(5, len(valid_genes) * 0.32)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# Scatter: size = pct expressing, color = scaled mean expression
size_scale = 200
sc_dot = ax.scatter(
    x           = dot_df["cluster"].cat.codes,
    y           = dot_df["gene"].cat.codes,
    s           = dot_df["pct_express"] * size_scale / 100,
    c           = dot_df["avg_scaled"],
    cmap        = "Blues",
    vmin        = 0, vmax = 1,
    edgecolors  = "grey",
    linewidths  = 0.3,
    alpha       = 0.9,
)

# Axes labels
ax.set_xticks(range(n_clusters))
ax.set_xticklabels(clusters, fontsize=9)
ax.set_yticks(range(len(valid_genes)))
ax.set_yticklabels(list(reversed(valid_genes)), fontsize=8, fontstyle="italic")
ax.set_xlabel("Cluster", fontsize=10)
ax.set_ylabel("Gene", fontsize=10)
ax.set_title(f"Top {N_TOP_MARKERS} marker genes per cluster",
             fontsize=12, fontweight="bold")
ax.text(0.5, 1.01,
        "Dot size = % cells expressing · Colour = scaled mean expression",
        ha="center", va="bottom", transform=ax.transAxes,
        fontsize=8, color="grey")

# Colour bar
cbar = plt.colorbar(sc_dot, ax=ax, fraction=0.03, pad=0.02, aspect=25)
cbar.set_label("Scaled\nexpression", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# Size legend
for pct in [10, 25, 50, 75]:
    ax.scatter([], [], c="grey", alpha=0.6,
               s=pct * size_scale / 100,
               label=f"{pct}%", edgecolors="grey", linewidths=0.3)
ax.legend(title="% expressing", fontsize=8, title_fontsize=8,
          loc="lower right", framealpha=0.8, edgecolor="none")

ax.set_facecolor("#F7F7F7")
ax.grid(axis="x", color="white", linewidth=0.8)
ax.grid(axis="y", color="white", linewidth=0.8)

plt.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT_DIR, f"dotplot_top_markers.{ext}"))
plt.close(fig)
print("Saved dotplot_top_markers.pdf/.png")

# ── 7. Cluster composition bar chart ─────────────────────────────────────────

print("Generating cluster composition bar chart...")

comp = (
    obs.groupby([CLUSTER_KEY, CONDITION_KEY], observed=True)
    .size()
    .reset_index(name="n_cells")
)
comp["total"] = comp.groupby(CLUSTER_KEY)["n_cells"].transform("sum")
comp["proportion"] = comp["n_cells"] / comp["total"]

# Pivot to wide for stacking
pivot = comp.pivot_table(
    index   = CLUSTER_KEY,
    columns = CONDITION_KEY,
    values  = "proportion",
    fill_value = 0
)
pivot = pivot.loc[clusters]                              # enforce cluster order
totals = comp.groupby(CLUSTER_KEY)["n_cells"].sum()

fig_w = max(8, n_clusters * 0.65)
fig, ax = plt.subplots(figsize=(fig_w, 6))

bar_width = 0.7
x = np.arange(n_clusters)

# Stack bars
bottom = np.zeros(n_clusters)
for cond in ["normal", "tumour"]:
    if cond not in pivot.columns:
        continue
    vals = pivot[cond].values
    bars = ax.bar(x, vals, bar_width,
                  bottom    = bottom,
                  color     = CONDITION_PAL[cond],
                  label     = "Normal adjacent" if cond == "normal" else "Tumour (T1 + T2)",
                  edgecolor = "white",
                  linewidth = 0.5)
    # Proportion labels inside bars (if segment > 8%)
    for xi, (v, b) in enumerate(zip(vals, bottom)):
        if v > 0.08:
            ax.text(xi, b + v / 2, f"{v:.0%}",
                    ha="center", va="center",
                    fontsize=7, color="white", fontweight="bold")
    bottom += vals

# Dashed 50% reference line
ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8, zorder=3)

# Total cell count above each bar
for xi, cl in enumerate(clusters):
    total = int(totals[cl])
    ax.text(xi, 1.025, f"{total:,}",
            ha="center", va="bottom", fontsize=7.5, color="#444444")

ax.set_xticks(x)
ax.set_xticklabels(clusters, fontsize=9)
ax.set_xlabel("Cluster", fontsize=10)
ax.set_ylabel("Proportion of cells", fontsize=10)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax.set_ylim(0, 1.08)
ax.set_title("Cluster composition — tumour vs normal",
             fontsize=12, fontweight="bold")
ax.text(0.5, 1.04,
        "Stacked proportions · numbers above bars = total cell count",
        ha="center", va="bottom", transform=ax.transAxes,
        fontsize=8, color="grey")
ax.legend(title="Condition", fontsize=9, title_fontsize=9,
          loc="upper right", framealpha=0.8, edgecolor="none")
ax.set_facecolor("#F7F7F7")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(color="white", linewidth=0.8, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT_DIR, f"cluster_composition.{ext}"))
plt.close(fig)
print("Saved cluster_composition.pdf/.png")

# ── 8. Full report figure ─────────────────────────────────────────────────────

print("Assembling full report figure...")

fig = plt.figure(figsize=(30, 20))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

def _umap_scatter(ax, colors, label_col=None, cmap=None, vmin=None, vmax=None,
                  title="", subtitle="", legend_elements=None):
    if cmap:
        sc_obj = ax.scatter(obs_expr["UMAP_1"], obs_expr["UMAP_2"],
                            c=obs_expr[MARKER_GENE], cmap=cmap, s=1.5,
                            alpha=0.8, linewidths=0, vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(sc_obj, ax=ax, fraction=0.04, pad=0.02, aspect=30)
        cbar.set_label(f"{MARKER_GENE}\n(log-norm)", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
    else:
        for key, col in colors.items():
            mask = obs[label_col] == key
            ax.scatter(obs.loc[mask, "UMAP_1"], obs.loc[mask, "UMAP_2"],
                       c=col, s=1.5, alpha=0.6, linewidths=0)
        if legend_elements:
            ax.legend(handles=legend_elements, fontsize=7, title_fontsize=8,
                      markerscale=1, loc="lower right",
                      framealpha=0.7, edgecolor="none")
    style_umap_ax(ax, title, subtitle)

# Row 0
ax0 = fig.add_subplot(gs[0, 0])
_umap_scatter(ax0, cluster_colors, CLUSTER_KEY,
              title="Cell clusters",
              subtitle=f"{len(obs):,} cells",
              legend_elements=legend_handles[:min(n_clusters, 20)])

ax1 = fig.add_subplot(gs[0, 1])
_umap_scatter(ax1, CONDITION_PAL, CONDITION_KEY,
              title="Condition",
              subtitle="Blue = normal · Red = tumour",
              legend_elements=[
                  Line2D([0], [0], marker="o", color="w",
                         markerfacecolor=c, markersize=6, label=k)
                  for k, c in CONDITION_PAL.items()
              ])

ax2 = fig.add_subplot(gs[0, 2])
_umap_scatter(ax2, {}, cmap="magma",
              vmin=0, vmax=vmax,
              title=f"{MARKER_GENE} expression",
              subtitle="log-norm, 99th pct clipped")

# Row 1: dot plot
ax3 = fig.add_subplot(gs[1, 0:2])
ax3.scatter(
    x          = dot_df["cluster"].cat.codes,
    y          = dot_df["gene"].cat.codes,
    s          = dot_df["pct_express"] * size_scale / 100,
    c          = dot_df["avg_scaled"],
    cmap       = "Blues",
    vmin=0, vmax=1,
    edgecolors = "grey", linewidths=0.3, alpha=0.9
)
ax3.set_xticks(range(n_clusters))
ax3.set_xticklabels(clusters, fontsize=8)
ax3.set_yticks(range(len(valid_genes)))
ax3.set_yticklabels(list(reversed(valid_genes)), fontsize=7, fontstyle="italic")
ax3.set_title(f"Top {N_TOP_MARKERS} marker genes per cluster",
              fontsize=11, fontweight="bold")
ax3.set_facecolor("#F7F7F7")

# Row 1: composition bar
ax4 = fig.add_subplot(gs[1, 2])
bottom = np.zeros(n_clusters)
for cond in ["normal", "tumour"]:
    if cond not in pivot.columns:
        continue
    vals = pivot[cond].values
    ax4.bar(x, vals, bar_width, bottom=bottom,
            color=CONDITION_PAL[cond], edgecolor="white", linewidth=0.5,
            label="Normal adjacent" if cond == "normal" else "Tumour")
    bottom += vals
ax4.axhline(0.5, color="grey", linestyle="--", linewidth=0.8)
ax4.set_xticks(x)
ax4.set_xticklabels(clusters, fontsize=8)
ax4.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
ax4.set_ylim(0, 1.1)
ax4.set_title("Cluster composition", fontsize=11, fontweight="bold")
ax4.legend(fontsize=8, framealpha=0.8, edgecolor="none")
ax4.set_facecolor("#F7F7F7")

fig.suptitle(
    f"scRNA-seq report · Sus scrofa NF1 neurofibromatosis\n"
    f"Samples: NormalAdjacent, Tumor1, Tumor2 · Marker: {MARKER_GENE}",
    fontsize=14, fontweight="bold", y=1.01
)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT_DIR, f"full_report.{ext}"),
                bbox_inches="tight", facecolor="white")
plt.close(fig)
print("Saved full_report.pdf/.png")

# ── Done ──────────────────────────────────────────────────────────────────────

print("\nDone. All figures written to:", OUT_DIR)
print("Files produced:")
for f in [
    "umap_combined.pdf/.png",
    "dotplot_top_markers.pdf/.png",
    "cluster_composition.pdf/.png",
    "full_report.pdf/.png",
    "top_cluster_markers.csv",
]:
    print(f"  {f}")
