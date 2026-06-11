#!/usr/bin/env python3
"""
==============================================================================
Glycine × Heat Stress 16S rRNA Microbiome Analysis Pipeline
Part 1: Data Loading, Preprocessing, Alpha/Beta Diversity
==============================================================================
"""
import os, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from skbio import TreeNode, DistanceMatrix
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim
import statsmodels.stats.multitest as smm
from itertools import combinations
import json, pickle

# ==== Paths ====
FEAT_FP = "/mnt/user-data/uploads/feature-table.tsv"
TAX_FP  = "/mnt/user-data/uploads/taxonomy.tsv"
META_FP = "/mnt/user-data/uploads/sample-metadata.txt"
TREE_FP = "/mnt/user-data/uploads/tree.nwk"
FIG_DIR = "/home/claude/results/figures"
TAB_DIR = "/home/claude/results/tables"

# ==== Style constants ====
TREAT_ORDER  = ["NC", "PC", "G25", "G50"]
TREAT_COLORS = {"NC": "#1976D2", "PC": "#388E3C", "G25": "#F57C00", "G50": "#D32F2F"}
FONT_FAMILY  = "Helvetica"

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.transparent': False,
})

# ===========================================================================
# 1) Load Data
# ===========================================================================
print("=" * 60)
print("[1] Loading data...")

# --- Feature table ---
feat = pd.read_csv(FEAT_FP, sep="\t", skiprows=1, index_col=0)
if "taxonomy" in feat.columns:
    feat = feat.drop(columns=["taxonomy"])
feat.columns = [c.replace("-", ".") for c in feat.columns]
feat = feat.astype(float).fillna(0).astype(int)

# --- Taxonomy ---
tax_raw = pd.read_csv(TAX_FP, sep="\t")
tax_raw.columns = ["FeatureID", "Taxon", "Confidence"]
tax_levels = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
tax_split = tax_raw["Taxon"].str.split(";", expand=True)
tax_split.columns = tax_levels[:tax_split.shape[1]]
for col in tax_levels:
    if col not in tax_split.columns:
        tax_split[col] = np.nan
tax_split = tax_split.apply(lambda x: x.str.strip())
tax_split = tax_split.replace({"": np.nan, "None": np.nan})
# Remove prefix like d__, p__, etc.
for col in tax_split.columns:
    tax_split[col] = tax_split[col].str.replace(r'^[a-z]__', '', regex=True)
    tax_split[col] = tax_split[col].replace({"": np.nan})
tax_split.index = tax_raw["FeatureID"]
tax_df = tax_split.copy()

# --- Metadata ---
meta = pd.read_csv(META_FP, sep="\t")
meta.columns = meta.columns.str.strip()
meta = meta.dropna(subset=["SampleID"])
meta["SampleID"] = meta["SampleID"].str.strip()
meta["Treat"] = meta["Treat"].str.strip()
meta = meta[meta["Treat"].isin(TREAT_ORDER)]
meta["Treat"] = pd.Categorical(meta["Treat"], categories=TREAT_ORDER, ordered=True)
meta = meta.set_index("SampleID")

# --- Tree ---
tree = TreeNode.read(TREE_FP)

# ===========================================================================
# 2) Match & Filter
# ===========================================================================
common_samples = sorted(set(feat.columns) & set(meta.index))
common_taxa = sorted(set(feat.index) & set(tax_df.index))

feat = feat.loc[common_taxa, common_samples]
tax_df = tax_df.loc[common_taxa]
meta = meta.loc[common_samples]

# Remove zero-sum taxa
keep = feat.sum(axis=1) > 0
feat = feat.loc[keep]
tax_df = tax_df.loc[feat.index]

# Prune tree to keep only matching tips
feat_ids_set = set(feat.index)
tree = tree.shear(feat_ids_set)

print(f"  Samples : {feat.shape[1]}")
print(f"  ASVs    : {feat.shape[0]}")
print(f"  Groups  : {meta['Treat'].value_counts().to_dict()}")

# ===========================================================================
# 3) Relative Abundance
# ===========================================================================
feat_rel = feat.div(feat.sum(axis=0), axis=1)

# Save processed data for later parts
pickle.dump({
    'feat': feat, 'feat_rel': feat_rel,
    'tax_df': tax_df, 'meta': meta, 'tree': tree
}, open("/home/claude/results/processed_data.pkl", "wb"))

# ===========================================================================
# 4) Alpha Diversity
# ===========================================================================
print("\n" + "=" * 60)
print("[2] Alpha Diversity...")

otu_mat = feat.T.values  # samples × ASVs
sample_ids = list(feat.columns)

# Calculate metrics
alpha_obs = alpha_diversity('sobs', otu_mat, ids=sample_ids)
alpha_sha = alpha_diversity('shannon', otu_mat, ids=sample_ids)
alpha_sim = alpha_diversity('simpson', otu_mat, ids=sample_ids)

# Faith's PD
from skbio.diversity import alpha as alpha_mod
# faith_pd requires tree
try:
    alpha_fpd = alpha_diversity('faith_pd', otu_mat, ids=sample_ids, tree=tree, taxa=list(feat.index))
except Exception:
    alpha_fpd = alpha_diversity('faith_pd', otu_mat, ids=sample_ids, tree=tree, otu_ids=list(feat.index))

alpha_df = pd.DataFrame({
    "SampleID": sample_ids,
    "Observed": alpha_obs.values,
    "Shannon": alpha_sha.values,
    "Simpson": alpha_sim.values,
    "Faith_PD": alpha_fpd.values
})
alpha_df = alpha_df.merge(meta.reset_index(), on="SampleID", how="left")
alpha_df.to_csv(os.path.join(TAB_DIR, "alpha_diversity.csv"), index=False)

# --- Kruskal-Wallis + Dunn's post-hoc ---
from scipy.stats import kruskal, mannwhitneyu

alpha_stats = []
metrics = ["Observed", "Shannon", "Simpson", "Faith_PD"]
metric_labels = {"Observed": "Observed ASVs", "Shannon": "Shannon Index",
                 "Simpson": "Simpson Index", "Faith_PD": "Faith's PD"}

for m in metrics:
    groups = [alpha_df.loc[alpha_df["Treat"] == t, m].values for t in TREAT_ORDER]
    h_stat, kw_p = kruskal(*groups)
    alpha_stats.append({"Metric": m, "H_statistic": h_stat, "KW_p_value": kw_p})
    # Pairwise Mann-Whitney
    pairs = list(combinations(TREAT_ORDER, 2))
    pvals = []
    for t1, t2 in pairs:
        v1 = alpha_df.loc[alpha_df["Treat"] == t1, m].values
        v2 = alpha_df.loc[alpha_df["Treat"] == t2, m].values
        if len(v1) > 1 and len(v2) > 1:
            _, p = mannwhitneyu(v1, v2, alternative='two-sided')
        else:
            p = np.nan
        pvals.append(p)
    # BH correction
    _, pvals_adj, _, _ = smm.multipletests(pvals, method='fdr_bh')
    for (t1, t2), raw_p, adj_p in zip(pairs, pvals, pvals_adj):
        alpha_stats.append({
            "Metric": m, "Comparison": f"{t1} vs {t2}",
            "MannWhitney_p": raw_p, "BH_adjusted_p": adj_p
        })

alpha_stats_df = pd.DataFrame(alpha_stats)
alpha_stats_df.to_csv(os.path.join(TAB_DIR, "alpha_diversity_stats.csv"), index=False)
print("  Alpha diversity stats saved.")

# Extract KW p-values for annotation
kw_pvals = {}
for m in metrics:
    row = alpha_stats_df[(alpha_stats_df["Metric"] == m) & (alpha_stats_df["H_statistic"].notna())]
    if len(row) > 0:
        kw_pvals[m] = row.iloc[0]["KW_p_value"]

# --- Alpha diversity figure ---
fig, axes = plt.subplots(1, 4, figsize=(16, 4.5))

for ax, m in zip(axes, metrics):
    data_plot = alpha_df[["Treat", m]].copy()
    
    # Boxplot
    bp = ax.boxplot(
        [data_plot.loc[data_plot["Treat"] == t, m].values for t in TREAT_ORDER],
        positions=range(len(TREAT_ORDER)),
        widths=0.55, patch_artist=True,
        showfliers=False, medianprops=dict(color='black', linewidth=1.5)
    )
    for patch, treat in zip(bp['boxes'], TREAT_ORDER):
        patch.set_facecolor(TREAT_COLORS[treat])
        patch.set_alpha(0.6)
    
    # Jitter points
    for i, t in enumerate(TREAT_ORDER):
        vals = data_plot.loc[data_plot["Treat"] == t, m].values
        jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(i + jitter, vals, color=TREAT_COLORS[t],
                   s=40, zorder=5, edgecolor='white', linewidth=0.5, alpha=0.9)
    
    # KW p-value annotation
    kw_p = kw_pvals.get(m, np.nan)
    if not np.isnan(kw_p):
        p_txt = f"KW p = {kw_p:.3f}" if kw_p >= 0.001 else f"KW p < 0.001"
        ax.set_title(f"{metric_labels[m]}\n{p_txt}", fontsize=11, fontweight='bold')
    else:
        ax.set_title(metric_labels[m], fontsize=11, fontweight='bold')
    
    ax.set_xticks(range(len(TREAT_ORDER)))
    ax.set_xticklabels(TREAT_ORDER, fontweight='bold')
    ax.set_ylabel(metric_labels[m])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.suptitle("Alpha Diversity", fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig1_AlphaDiversity.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig1_AlphaDiversity.pdf"))
plt.close()
print("  Alpha diversity figure saved.")

# ===========================================================================
# 5) Beta Diversity
# ===========================================================================
print("\n" + "=" * 60)
print("[3] Beta Diversity...")

dist_methods = {
    "Bray-Curtis": None,
    "Unweighted UniFrac": None,
    "Weighted UniFrac": None,
}

# Calculate distance matrices
dm_bc = beta_diversity('braycurtis', otu_mat, ids=sample_ids)
dm_uu = beta_diversity('unweighted_unifrac', otu_mat, ids=sample_ids,
                        tree=tree, taxa=list(feat.index))
dm_wu = beta_diversity('weighted_unifrac', otu_mat, ids=sample_ids,
                        tree=tree, taxa=list(feat.index))

dist_methods["Bray-Curtis"] = dm_bc
dist_methods["Unweighted UniFrac"] = dm_uu
dist_methods["Weighted UniFrac"] = dm_wu

# --- PERMANOVA ---
grouping = meta.loc[sample_ids, "Treat"].astype(str)

permanova_results = []
for name, dm in dist_methods.items():
    res = permanova(dm, grouping, permutations=9999)
    permanova_results.append({
        "Distance": name,
        "pseudo_F": res["test statistic"],
        "p_value": res["p-value"],
        "n_permutations": 9999
    })
    print(f"  PERMANOVA {name}: F={res['test statistic']:.3f}, p={res['p-value']:.4f}")

# Pairwise PERMANOVA
pairwise_permanova = []
for name, dm in dist_methods.items():
    pairs = list(combinations(TREAT_ORDER, 2))
    pvals_raw = []
    for t1, t2 in pairs:
        idx = [s for s in sample_ids if meta.loc[s, "Treat"] in [t1, t2]]
        dm_sub = dm.filter(idx)
        grp_sub = meta.loc[idx, "Treat"].astype(str)
        res = permanova(dm_sub, grp_sub, permutations=9999)
        pvals_raw.append(res["p-value"])
        pairwise_permanova.append({
            "Distance": name, "Comparison": f"{t1} vs {t2}",
            "pseudo_F": res["test statistic"], "p_value": res["p-value"]
        })
    # BH correction
    _, adj_pvals, _, _ = smm.multipletests(pvals_raw, method='fdr_bh')
    for i, (t1, t2) in enumerate(pairs):
        pairwise_permanova[-len(pairs)+i]["BH_adjusted_p"] = adj_pvals[i]

pd.DataFrame(permanova_results).to_csv(os.path.join(TAB_DIR, "PERMANOVA_global.csv"), index=False)
pd.DataFrame(pairwise_permanova).to_csv(os.path.join(TAB_DIR, "PERMANOVA_pairwise.csv"), index=False)

# --- ANOSIM ---
anosim_results = []
for name, dm in dist_methods.items():
    res = anosim(dm, grouping, permutations=9999)
    anosim_results.append({
        "Distance": name, "R_statistic": res["test statistic"], "p_value": res["p-value"]
    })
pd.DataFrame(anosim_results).to_csv(os.path.join(TAB_DIR, "ANOSIM.csv"), index=False)

# --- PCoA Figure ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

for ax, (name, dm) in zip(axes, dist_methods.items()):
    pc = pcoa(dm)
    prop_exp = pc.proportion_explained
    coords = pc.samples[["PC1", "PC2"]]
    coords.index = sample_ids
    coords["Treat"] = meta.loc[sample_ids, "Treat"].values
    
    # Get PERMANOVA p-value
    perm_p = [r for r in permanova_results if r["Distance"] == name][0]["p_value"]
    perm_f = [r for r in permanova_results if r["Distance"] == name][0]["pseudo_F"]
    
    for treat in TREAT_ORDER:
        sub = coords[coords["Treat"] == treat]
        color = TREAT_COLORS[treat]
        
        # Confidence ellipse with fill
        if len(sub) >= 3:
            from matplotlib.patches import Ellipse
            import matplotlib.transforms as transforms
            
            mean_x, mean_y = sub["PC1"].mean(), sub["PC2"].mean()
            cov_mat = np.cov(sub["PC1"].values, sub["PC2"].values)
            eigenvalues, eigenvectors = np.linalg.eigh(cov_mat)
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
            
            # 95% CI ellipse (chi-squared 2df, alpha=0.05 → 5.991)
            chi2_val = 5.991
            width = 2 * np.sqrt(eigenvalues[0] * chi2_val)
            height = 2 * np.sqrt(eigenvalues[1] * chi2_val)
            
            ellipse = Ellipse(xy=(mean_x, mean_y), width=width, height=height,
                             angle=angle, facecolor=color, alpha=0.15,
                             edgecolor=color, linewidth=1.5, linestyle='--')
            ax.add_patch(ellipse)
        
        # Points
        ax.scatter(sub["PC1"], sub["PC2"], c=color, s=70, label=treat,
                  edgecolor='white', linewidth=0.5, zorder=5, alpha=0.9)
    
    ax.set_xlabel(f"PCoA1 ({prop_exp.iloc[0]*100:.1f}%)", fontweight='bold')
    ax.set_ylabel(f"PCoA2 ({prop_exp.iloc[1]*100:.1f}%)", fontweight='bold')
    p_str = f"p = {perm_p:.4f}" if perm_p >= 0.001 else "p < 0.001"
    ax.set_title(f"{name}\nPERMANOVA: F = {perm_f:.2f}, {p_str}",
                fontsize=11, fontweight='bold')
    ax.legend(title="Treatment", frameon=True, fontsize=9, title_fontsize=10)
    ax.axhline(0, color='grey', linewidth=0.5, alpha=0.3)
    ax.axvline(0, color='grey', linewidth=0.5, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig2_BetaDiversity_PCoA.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig2_BetaDiversity_PCoA.pdf"))
plt.close()
print("  Beta diversity figures saved.")

# --- Save distance matrices ---
for name, dm in dist_methods.items():
    fn = name.replace(" ", "_").replace("-", "")
    dm_df = dm.to_data_frame()
    dm_df.to_csv(os.path.join(TAB_DIR, f"distance_matrix_{fn}.csv"))

# Save PCoA coordinates
pcoa_all = {}
for name, dm in dist_methods.items():
    pc = pcoa(dm)
    coords = pc.samples.copy()
    coords["Treat"] = meta.loc[sample_ids, "Treat"].values
    pcoa_all[name] = coords
    fn = name.replace(" ", "_").replace("-", "")
    coords.to_csv(os.path.join(TAB_DIR, f"PCoA_coords_{fn}.csv"))

print("\n[DONE] Part 1 complete.")
