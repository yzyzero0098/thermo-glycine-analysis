#!/usr/bin/env python3
"""
==============================================================================
Part 2: Taxonomy Composition, Differential Abundance, Centroid Distance,
        Glycine Mimicking Analysis
==============================================================================
"""
import os, warnings, pickle
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import braycurtis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
import statsmodels.stats.multitest as smm
from itertools import combinations

# ==== Load processed data ====
data = pickle.load(open("/home/claude/results/processed_data.pkl", "rb"))
feat, feat_rel, tax_df, meta, tree = data['feat'], data['feat_rel'], data['tax_df'], data['meta'], data['tree']

FIG_DIR = "/home/claude/results/figures"
TAB_DIR = "/home/claude/results/tables"
TREAT_ORDER  = ["NC", "PC", "G25", "G50"]
TREAT_COLORS = {"NC": "#1976D2", "PC": "#388E3C", "G25": "#F57C00", "G50": "#D32F2F"}
sample_ids = list(feat.columns)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titleweight': 'bold',
    'axes.labelweight': 'bold',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ===========================================================================
# 6) Taxonomy Composition Barplot — Phylum & Genus level
# ===========================================================================
print("=" * 60)
print("[4] Taxonomy Composition...")

def aggregate_taxa(feat_rel, tax_df, rank, top_n=10):
    """Aggregate relative abundance to a given rank."""
    tax_rank = tax_df[rank].fillna("Unclassified").copy()
    # Build aggregation df
    agg = feat_rel.copy()
    agg['taxon'] = tax_rank.loc[agg.index].values
    agg = agg.groupby('taxon').sum()
    # Top N
    total = agg.sum(axis=1).sort_values(ascending=False)
    top = total.head(top_n).index.tolist()
    # Merge Others
    agg_top = agg.loc[agg.index.isin(top)].copy()
    others = agg.loc[~agg.index.isin(top)].sum(axis=0)
    agg_top.loc["Others"] = others
    return agg_top

# --- Phylum ---
phylum_agg = aggregate_taxa(feat_rel, tax_df, "Phylum", top_n=8)
# --- Genus ---
genus_agg = aggregate_taxa(feat_rel, tax_df, "Genus", top_n=15)

# Color palette for taxa
def get_taxa_palette(taxa_list):
    # Use a qualitative colormap
    base_colors = [
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
        "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
        "#E377C2", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78",
        "#98DF8A"
    ]
    palette = {}
    for i, t in enumerate(taxa_list):
        if t == "Others":
            palette[t] = "#CCCCCC"
        elif t == "Unclassified":
            palette[t] = "#999999"
        else:
            palette[t] = base_colors[i % len(base_colors)]
    return palette

# Sort samples by treatment
sample_order = []
for t in TREAT_ORDER:
    samps = meta[meta["Treat"] == t].index.tolist()
    sample_order.extend(sorted(samps))

# --- Combined Phylum + Genus figure ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))

for ax, agg, title, top_n in [(ax1, phylum_agg, "Phylum", 8), (ax2, genus_agg, "Genus", 15)]:
    agg_ordered = agg[sample_order]
    taxa_order = agg_ordered.sum(axis=1).sort_values(ascending=True).index.tolist()
    if "Others" in taxa_order:
        taxa_order.remove("Others")
        taxa_order = ["Others"] + taxa_order
    
    palette = get_taxa_palette([t for t in taxa_order if t != "Others"] + ["Others"])
    
    bottom = np.zeros(len(sample_order))
    for taxon in taxa_order:
        vals = agg_ordered.loc[taxon].values
        ax.bar(range(len(sample_order)), vals, bottom=bottom,
               color=palette.get(taxon, "#CCCCCC"), label=taxon, width=0.85, edgecolor='white', linewidth=0.3)
        bottom += vals
    
    # Treatment separators
    cumsum = 0
    for t in TREAT_ORDER:
        n = (meta["Treat"] == t).sum()
        mid = cumsum + n / 2 - 0.5
        ax.text(mid, 1.02, t, ha='center', fontweight='bold', fontsize=11,
                color=TREAT_COLORS[t], transform=ax.get_xaxis_transform())
        cumsum += n
        if t != TREAT_ORDER[-1]:
            ax.axvline(cumsum - 0.5, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
    
    ax.set_ylabel("Relative Abundance", fontweight='bold')
    ax.set_title(f"{title}-level Composition", fontweight='bold', fontsize=12)
    ax.set_xlim(-0.5, len(sample_order) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks(range(len(sample_order)))
    ax.set_xticklabels([s.split('.')[-2] for s in sample_order], rotation=45, ha='right', fontsize=7)
    
    handles, labels = ax.get_legend_handles_labels()
    # Reverse to match stacking order
    ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=8, frameon=False, title=title, title_fontsize=9)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig3_Taxonomy_Composition.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig3_Taxonomy_Composition.pdf"))
plt.close()
print("  Taxonomy composition figure saved.")

# Save taxonomy tables
phylum_agg[sample_order].T.to_csv(os.path.join(TAB_DIR, "taxonomy_phylum_relabund.csv"))
genus_agg[sample_order].T.to_csv(os.path.join(TAB_DIR, "taxonomy_genus_relabund.csv"))

# ===========================================================================
# 7) Genus-level mean ± SE per group (heatmap)
# ===========================================================================
genus_mean = pd.DataFrame()
genus_se = pd.DataFrame()
for t in TREAT_ORDER:
    samps = meta[meta["Treat"] == t].index.tolist()
    genus_mean[t] = genus_agg[samps].mean(axis=1)
    genus_se[t] = genus_agg[samps].std(axis=1) / np.sqrt(len(samps))

genus_mean.to_csv(os.path.join(TAB_DIR, "genus_mean_relabund.csv"))
genus_se.to_csv(os.path.join(TAB_DIR, "genus_se_relabund.csv"))

# Heatmap of top genera (excluding Others)
hm_data = genus_mean.drop("Others", errors='ignore')
hm_data = hm_data.loc[hm_data.max(axis=1) > 0.005]  # >0.5% at least
hm_data = hm_data.sort_values(TREAT_ORDER[0], ascending=False)

fig, ax = plt.subplots(figsize=(6, max(4, len(hm_data)*0.4)))
sns.heatmap(hm_data, cmap='YlOrRd', annot=True, fmt='.3f',
            linewidths=0.5, linecolor='white', ax=ax,
            cbar_kws={'label': 'Relative Abundance'})
ax.set_title("Mean Relative Abundance by Treatment (Genus)", fontweight='bold')
ax.set_ylabel("")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig4_Genus_Heatmap.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig4_Genus_Heatmap.pdf"))
plt.close()
print("  Genus heatmap saved.")

# ===========================================================================
# 8) Differential Abundance — Kruskal-Wallis + pairwise at Genus level
# ===========================================================================
print("\n" + "=" * 60)
print("[5] Differential Abundance (Genus level)...")

genus_full = aggregate_taxa(feat_rel, tax_df, "Genus", top_n=999)
genus_full = genus_full.drop("Others", errors='ignore')
genus_full = genus_full.loc[genus_full.mean(axis=1) >= 0.001]  # ≥0.1% prevalence

da_results = []
for taxon in genus_full.index:
    groups = {}
    for t in TREAT_ORDER:
        samps = meta[meta["Treat"] == t].index.tolist()
        groups[t] = genus_full.loc[taxon, samps].values
    
    # Kruskal-Wallis
    try:
        h, kw_p = stats.kruskal(*groups.values())
    except:
        h, kw_p = np.nan, np.nan
    
    row = {"Genus": taxon, "KW_H": h, "KW_p": kw_p}
    for t in TREAT_ORDER:
        row[f"mean_{t}"] = np.mean(groups[t])
        row[f"se_{t}"] = np.std(groups[t]) / np.sqrt(len(groups[t]))
    da_results.append(row)

da_df = pd.DataFrame(da_results)
_, da_df["KW_BH_p"], _, _ = smm.multipletests(da_df["KW_p"].fillna(1), method='fdr_bh')
da_df = da_df.sort_values("KW_p")

# Pairwise for significant genera
pairs = list(combinations(TREAT_ORDER, 2))
pairwise_list = []
for _, row in da_df[da_df["KW_p"] < 0.1].iterrows():
    taxon = row["Genus"]
    for t1, t2 in pairs:
        s1 = meta[meta["Treat"] == t1].index.tolist()
        s2 = meta[meta["Treat"] == t2].index.tolist()
        v1 = genus_full.loc[taxon, s1].values
        v2 = genus_full.loc[taxon, s2].values
        try:
            _, p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
        except:
            p = np.nan
        pairwise_list.append({
            "Genus": taxon, "Comparison": f"{t1} vs {t2}", "MWU_p": p
        })

pairwise_da = pd.DataFrame(pairwise_list)
if len(pairwise_da) > 0:
    _, pairwise_da["BH_p"], _, _ = smm.multipletests(pairwise_da["MWU_p"].fillna(1), method='fdr_bh')

da_df.to_csv(os.path.join(TAB_DIR, "DA_KruskalWallis_Genus.csv"), index=False)
pairwise_da.to_csv(os.path.join(TAB_DIR, "DA_Pairwise_Genus.csv"), index=False)

sig_genera = da_df[da_df["KW_p"] < 0.1]["Genus"].tolist()
print(f"  Genera with KW p < 0.1: {len(sig_genera)}")
for g in sig_genera[:10]:
    p = da_df[da_df["Genus"] == g]["KW_p"].values[0]
    print(f"    {g}: p = {p:.4f}")

# ===========================================================================
# 9) Centroid Distance Analysis — Glycine mimicking score
# ===========================================================================
print("\n" + "=" * 60)
print("[6] Centroid Distance Analysis...")

otu_mat = feat.T.values
dm_bc = beta_diversity('braycurtis', otu_mat, ids=sample_ids)
dm_wu = beta_diversity('weighted_unifrac', otu_mat, ids=sample_ids,
                        tree=tree, taxa=list(feat.index))

centroid_results = []
for dm_name, dm in [("Bray-Curtis", dm_bc), ("Weighted_UniFrac", dm_wu)]:
    dm_df = dm.to_data_frame()
    
    for t1, t2 in combinations(TREAT_ORDER, 2):
        s1 = meta[meta["Treat"] == t1].index.tolist()
        s2 = meta[meta["Treat"] == t2].index.tolist()
        
        # Inter-group distances
        dists = []
        for a in s1:
            for b in s2:
                dists.append(dm_df.loc[a, b])
        
        centroid_results.append({
            "Distance_metric": dm_name,
            "Group1": t1, "Group2": t2,
            "Mean_distance": np.mean(dists),
            "SD_distance": np.std(dists),
            "n_comparisons": len(dists)
        })

centroid_df = pd.DataFrame(centroid_results)
centroid_df.to_csv(os.path.join(TAB_DIR, "centroid_intergroup_distances.csv"), index=False)

# Distance to PC (key metric)
print("\n  Distance to PC group (mean ± SD):")
for dm_name in ["Bray-Curtis", "Weighted_UniFrac"]:
    print(f"\n  [{dm_name}]")
    sub = centroid_df[(centroid_df["Distance_metric"] == dm_name) & (centroid_df["Group2"] == "PC")]
    sub2 = centroid_df[(centroid_df["Distance_metric"] == dm_name) & (centroid_df["Group1"] == "PC")]
    dist_to_pc = pd.concat([sub, sub2])
    for _, row in dist_to_pc.iterrows():
        g1, g2 = row["Group1"], row["Group2"]
        other = g2 if g1 == "PC" else g1
        print(f"    {other} ↔ PC: {row['Mean_distance']:.4f} ± {row['SD_distance']:.4f}")

# --- Centroid distance barplot ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, dm_name in zip(axes, ["Bray-Curtis", "Weighted_UniFrac"]):
    sub = centroid_df[centroid_df["Distance_metric"] == dm_name]
    # Filter: distance from NC, G25, G50 to PC
    pc_dists = []
    for grp in ["NC", "G25", "G50"]:
        row = sub[((sub["Group1"] == grp) & (sub["Group2"] == "PC")) |
                  ((sub["Group1"] == "PC") & (sub["Group2"] == grp))]
        if len(row) > 0:
            pc_dists.append({"Group": grp, "Distance": row.iloc[0]["Mean_distance"],
                           "SD": row.iloc[0]["SD_distance"]})
    
    pc_df = pd.DataFrame(pc_dists)
    bars = ax.bar(pc_df["Group"], pc_df["Distance"],
                  color=[TREAT_COLORS[g] for g in pc_df["Group"]],
                  edgecolor='black', linewidth=0.5, alpha=0.8,
                  yerr=pc_df["SD"], capsize=5)
    
    # Highlight minimum
    min_idx = pc_df["Distance"].idxmin()
    min_grp = pc_df.loc[min_idx, "Group"]
    bars[min_idx].set_edgecolor('#D32F2F')
    bars[min_idx].set_linewidth(2.5)
    
    ax.set_ylabel("Mean Distance to PC", fontweight='bold')
    ax.set_title(f"{dm_name}\nDistance to Energy Control (PC)", fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Annotate minimum
    ax.annotate(f"Min: {min_grp}", xy=(min_idx, pc_df.loc[min_idx, "Distance"]),
               xytext=(0, 15), textcoords='offset points', ha='center',
               fontweight='bold', color='#D32F2F', fontsize=10,
               arrowprops=dict(arrowstyle='->', color='#D32F2F'))

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig5_Centroid_Distance_to_PC.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig5_Centroid_Distance_to_PC.pdf"))
plt.close()
print("  Centroid distance figure saved.")

# ===========================================================================
# 10) Glycine Mimicking Taxa Analysis
# ===========================================================================
print("\n" + "=" * 60)
print("[7] Glycine Mimicking Taxa...")

# Log fold change calculation (using mean relative abundance + pseudocount)
pseudo = 1e-6

def calc_lfc(agg_df, g1_samples, g2_samples):
    """log2 fold change: g2/g1"""
    mean1 = agg_df[g1_samples].mean(axis=1) + pseudo
    mean2 = agg_df[g2_samples].mean(axis=1) + pseudo
    return np.log2(mean2 / mean1)

nc_samps = meta[meta["Treat"] == "NC"].index.tolist()
pc_samps = meta[meta["Treat"] == "PC"].index.tolist()
g25_samps = meta[meta["Treat"] == "G25"].index.tolist()
g50_samps = meta[meta["Treat"] == "G50"].index.tolist()

lfc_pc_nc = calc_lfc(genus_full, nc_samps, pc_samps)
lfc_g25_nc = calc_lfc(genus_full, nc_samps, g25_samps)
lfc_g50_nc = calc_lfc(genus_full, nc_samps, g50_samps)

mimick_df = pd.DataFrame({
    "Genus": genus_full.index,
    "LFC_PC_vs_NC": lfc_pc_nc.values,
    "LFC_G25_vs_NC": lfc_g25_nc.values,
    "LFC_G50_vs_NC": lfc_g50_nc.values,
})
mimick_df["Delta_G25"] = np.abs(mimick_df["LFC_G25_vs_NC"] - mimick_df["LFC_PC_vs_NC"])
mimick_df["Delta_G50"] = np.abs(mimick_df["LFC_G50_vs_NC"] - mimick_df["LFC_PC_vs_NC"])
mimick_df["Min_delta"] = mimick_df[["Delta_G25", "Delta_G50"]].min(axis=1)
mimick_df["Best_dose"] = mimick_df[["Delta_G25", "Delta_G50"]].idxmin(axis=1).str.replace("Delta_", "")
mimick_df = mimick_df.sort_values("Min_delta")
mimick_df.to_csv(os.path.join(TAB_DIR, "glycine_mimicking_taxa.csv"), index=False)

# --- Mimicking heatmap (top 20) ---
top_mimic = mimick_df.head(20).set_index("Genus")
hm_cols = ["LFC_PC_vs_NC", "LFC_G25_vs_NC", "LFC_G50_vs_NC"]
hm_labels = ["PC vs NC", "G25 vs NC", "G50 vs NC"]

fig, ax = plt.subplots(figsize=(7, max(5, len(top_mimic)*0.35)))
vmax = max(abs(top_mimic[hm_cols].values.min()), abs(top_mimic[hm_cols].values.max()))
sns.heatmap(top_mimic[hm_cols], cmap='RdBu_r', center=0, vmin=-vmax, vmax=vmax,
            annot=True, fmt='.2f', linewidths=0.5, linecolor='white',
            xticklabels=hm_labels, ax=ax,
            cbar_kws={'label': 'Log₂ Fold Change'})
ax.set_title("Top 20 Glycine-Mimicking Taxa\n(Smallest |ΔlogFC| to PC)", fontweight='bold')
ax.set_ylabel("")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig6_Mimicking_Taxa_Heatmap.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig6_Mimicking_Taxa_Heatmap.pdf"))
plt.close()
print("  Mimicking taxa heatmap saved.")

# ===========================================================================
# 11) Dose-Response Trend (Jonckheere-Terpstra style with Spearman)
# ===========================================================================
print("\n" + "=" * 60)
print("[8] Dose-Response Trend Analysis...")

# Glycine dose axis: NC=0, G25=0.341, G50=0.683
dose_map = {"NC": 0, "G25": 0.341, "G50": 0.683}
dose_samples = meta[meta["Treat"].isin(["NC", "G25", "G50"])].copy()
dose_samples["Dose"] = dose_samples["Treat"].map(dose_map)

trend_results = []
for taxon in genus_full.index:
    vals = genus_full.loc[taxon, dose_samples.index].values
    doses = dose_samples["Dose"].values
    rho, p = stats.spearmanr(doses, vals)
    trend_results.append({"Genus": taxon, "Spearman_rho": rho, "Spearman_p": p})

trend_df = pd.DataFrame(trend_results)
_, trend_df["BH_p"], _, _ = smm.multipletests(trend_df["Spearman_p"].fillna(1), method='fdr_bh')
trend_df = trend_df.sort_values("Spearman_p")
trend_df.to_csv(os.path.join(TAB_DIR, "dose_response_trend.csv"), index=False)

sig_trend = trend_df[trend_df["Spearman_p"] < 0.1]
print(f"  Taxa with dose-response trend (p < 0.1): {len(sig_trend)}")
for _, r in sig_trend.head(10).iterrows():
    print(f"    {r['Genus']}: rho={r['Spearman_rho']:.3f}, p={r['Spearman_p']:.4f}")

# --- Trend figure for top taxa ---
trend_top = trend_df.head(12)
fig, axes = plt.subplots(3, 4, figsize=(16, 10))
axes = axes.flatten()

for i, (_, row) in enumerate(trend_top.iterrows()):
    if i >= 12:
        break
    ax = axes[i]
    taxon = row["Genus"]
    
    for t in ["NC", "G25", "G50"]:
        samps = meta[meta["Treat"] == t].index.tolist()
        vals = genus_full.loc[taxon, samps].values
        dose = dose_map[t]
        jitter = np.random.uniform(-0.015, 0.015, len(vals))
        ax.scatter([dose + j for j in jitter], vals, color=TREAT_COLORS[t],
                  s=50, zorder=5, alpha=0.8, edgecolor='white', linewidth=0.5)
    
    # Mean trend line
    means = [genus_full.loc[taxon, meta[meta["Treat"] == t].index].mean() for t in ["NC", "G25", "G50"]]
    ax.plot([0, 0.341, 0.683], means, 'k--', linewidth=1, alpha=0.6)
    
    p_str = f"p={row['Spearman_p']:.3f}" if row['Spearman_p'] >= 0.001 else "p<0.001"
    ax.set_title(f"{taxon}\nρ={row['Spearman_rho']:.2f}, {p_str}", fontsize=9, fontweight='bold')
    ax.set_xlabel("Gly dose (%)")
    ax.set_ylabel("Rel. Abund.")
    ax.set_xticks([0, 0.341, 0.683])
    ax.set_xticklabels(["NC", "G25", "G50"], fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle("Dose-Response Trend (Spearman Correlation, Glycine Axis)", fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig7_DoseResponse_Trend.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig7_DoseResponse_Trend.pdf"))
plt.close()
print("  Dose-response figure saved.")

# ===========================================================================
# 12) Phylum-level differential abundance boxplot
# ===========================================================================
print("\n" + "=" * 60)
print("[9] Phylum-level Boxplots...")

phylum_full = aggregate_taxa(feat_rel, tax_df, "Phylum", top_n=999)
phylum_full = phylum_full.drop("Others", errors='ignore')
top_phyla = phylum_full.mean(axis=1).sort_values(ascending=False).head(6).index.tolist()

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()

for i, phylum in enumerate(top_phyla):
    ax = axes[i]
    plot_data = []
    for t in TREAT_ORDER:
        samps = meta[meta["Treat"] == t].index.tolist()
        vals = phylum_full.loc[phylum, samps].values
        for v in vals:
            plot_data.append({"Treat": t, "RelAbund": v})
    pdf = pd.DataFrame(plot_data)
    
    bp = ax.boxplot(
        [pdf[pdf["Treat"] == t]["RelAbund"].values for t in TREAT_ORDER],
        positions=range(len(TREAT_ORDER)), widths=0.55, patch_artist=True,
        showfliers=False, medianprops=dict(color='black', linewidth=1.5)
    )
    for patch, t in zip(bp['boxes'], TREAT_ORDER):
        patch.set_facecolor(TREAT_COLORS[t])
        patch.set_alpha(0.6)
    for j, t in enumerate(TREAT_ORDER):
        vals = pdf[pdf["Treat"] == t]["RelAbund"].values
        jitter = np.random.uniform(-0.1, 0.1, len(vals))
        ax.scatter(j + jitter, vals, color=TREAT_COLORS[t], s=35, zorder=5,
                  edgecolor='white', linewidth=0.5, alpha=0.9)
    
    # KW test
    groups = [pdf[pdf["Treat"] == t]["RelAbund"].values for t in TREAT_ORDER]
    try:
        _, kw_p = stats.kruskal(*groups)
        p_str = f"KW p={kw_p:.3f}" if kw_p >= 0.001 else "KW p<0.001"
    except:
        p_str = ""
    
    ax.set_title(f"p. {phylum}\n{p_str}", fontsize=10, fontweight='bold', style='italic')
    ax.set_xticks(range(len(TREAT_ORDER)))
    ax.set_xticklabels(TREAT_ORDER, fontweight='bold')
    ax.set_ylabel("Relative Abundance")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.suptitle("Phylum-Level Relative Abundance", fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig8_Phylum_Boxplots.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig8_Phylum_Boxplots.pdf"))
plt.close()
print("  Phylum boxplot saved.")

# ===========================================================================
# 13) Within-group vs Between-group beta-dispersion
# ===========================================================================
print("\n" + "=" * 60)
print("[10] Beta Dispersion...")

dm_bc_df = dm_bc.to_data_frame()
disp_results = []

for t in TREAT_ORDER:
    samps = meta[meta["Treat"] == t].index.tolist()
    # Within-group distances
    within = []
    for s1, s2 in combinations(samps, 2):
        within.append(dm_bc_df.loc[s1, s2])
    disp_results.append({
        "Group": t,
        "Mean_within_dist": np.mean(within) if within else np.nan,
        "SD_within_dist": np.std(within) if within else np.nan,
        "n_pairs": len(within)
    })

pd.DataFrame(disp_results).to_csv(os.path.join(TAB_DIR, "beta_dispersion.csv"), index=False)

print("\n[DONE] Part 2 complete.")
print(f"\nAll figures saved in: {FIG_DIR}")
print(f"All tables saved in:  {TAB_DIR}")
