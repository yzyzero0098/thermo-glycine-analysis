#!/usr/bin/env python3
"""
==============================================================================
Glycine × Heat Stress: Phenotype–Microbiome Integration
==============================================================================
논문(Nam et al., 2023 Poultry Science)의 phenotype 데이터를 
jejunal 16S microbiome 데이터와 통합하여 정리 및 분석
==============================================================================
"""
import os, warnings, pickle
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova
import statsmodels.stats.multitest as smm

# ==== Load processed data ====
data = pickle.load(open("/home/claude/results/processed_data.pkl", "rb"))
feat = data['feat']
feat_rel = data['feat_rel']
tax_df = data['tax_df']
meta = data['meta']
tree = data['tree']

FIG_DIR = "/home/claude/results/figures"
TAB_DIR = "/home/claude/results/tables"
TREAT_ORDER = ["NC", "PC", "G25", "G50"]
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
# 1) Sample ID 구조 파싱 — 개체(replicate) 매칭 정보 추출
# ===========================================================================
print("=" * 70)
print("[1] Sample ID 구조 파싱 및 Replicate 매칭...")

meta_extended = meta.copy()
meta_extended["SampleID"] = meta_extended.index

# Parse: P1.X.JC.Treatment → Replicate = X
meta_extended["Replicate"] = meta_extended["SampleID"].apply(
    lambda x: x.split(".")[1]  # 두 번째 element가 replicate 번호
)
meta_extended["Tissue"] = "JC"  # Jejunum Content

print("\n  Sample-Replicate 매칭 현황:")
for t in TREAT_ORDER:
    reps = sorted(meta_extended[meta_extended["Treat"] == t]["Replicate"].tolist())
    print(f"    {t}: replicates {reps} (n={len(reps)})")

print(f"\n  * 논문: 7 replicates/group × 4 groups = 28 total")
print(f"  * 시퀀싱 후 QC 통과: {len(meta_extended)} samples")
print(f"  * Missing samples은 sequencing depth 또는 QC failure로 추정")

# ===========================================================================
# 2) 논문 Phenotype 데이터 구조화
# ===========================================================================
print("\n" + "=" * 70)
print("[2] 논문 Phenotype 데이터 구조화...")

# ---- A. 개체 수준 측정 (Individual-level: 도축 시 1수/replicate) ----
# 논문에서 "1 hen with BW close to average in each replicate was selected"
# → 이 개체들이 바로 microbiome sample을 제공한 개체!

individual_phenotypes = {
    # Serum Measurements (Table 4)
    "AST_UL": {"NC": 306.0, "PC": 232.1, "G25": 200.7, "G50": 269.3,
               "SEM": 37.70, "unit": "U/L", "source": "Table 4"},
    "ALT_UL": {"NC": 4.00, "PC": 2.40, "G25": 1.87, "G50": 3.65,
               "SEM": 0.783, "unit": "U/L", "source": "Table 4"},
    "Uric_acid_mg_dL": {"NC": 5.30, "PC": 3.67, "G25": 3.92, "G50": 5.96,
                         "SEM": 0.905, "unit": "mg/dL", "source": "Table 4"},
    "Creatinine_mg_dL": {"NC": 0.25, "PC": 0.25, "G25": 0.23, "G50": 0.27,
                          "SEM": 0.016, "unit": "mg/dL", "source": "Table 4"},
    
    # Stress Indicators (Table 5)
    "Blood_HL_ratio": {"NC": 0.25, "PC": 0.18, "G25": 0.25, "G50": 0.16,
                        "SEM": 0.027, "unit": "ratio", "source": "Table 5"},
    "Feather_CORT_ng_g": {"NC": 7.9, "PC": 11.5, "G25": 6.8, "G50": 6.1,
                           "SEM": 2.50, "unit": "ng/g", "source": "Table 5"},
    
    # Fatty Liver Incidence (Table 6) — 같은 개체의 간
    "Liver_color_score": {"NC": 3.04, "PC": 1.50, "G25": 1.43, "G50": 1.33,
                           "SEM": 0.242, "unit": "score 1-5", "source": "Table 6"},
    "Liver_L_star": {"NC": 32.0, "PC": 25.5, "G25": 28.6, "G50": 25.7,
                      "SEM": 1.05, "unit": "CIE L*", "source": "Table 6"},
    "Liver_a_star": {"NC": 20.7, "PC": 14.4, "G25": 18.5, "G50": 16.2,
                      "SEM": 1.19, "unit": "CIE a*", "source": "Table 6"},
    "Liver_b_star": {"NC": 14.8, "PC": 7.3, "G25": 10.5, "G50": 7.4,
                      "SEM": 1.28, "unit": "CIE b*", "source": "Table 6"},
    "Hemorrhagic_score": {"NC": 1.39, "PC": 0.68, "G25": 1.04, "G50": 0.58,
                           "SEM": 0.218, "unit": "score 0-5", "source": "Table 6"},
    "Liver_fat_pct_DM": {"NC": 29.5, "PC": 15.1, "G25": 22.1, "G50": 19.3,
                          "SEM": 3.23, "unit": "% DM", "source": "Table 6"},
}

# ---- B. Replicate 수준 측정 (12주 평균) ----
replicate_phenotypes = {
    "FI_g_hen_d": {"NC": 96, "PC": 91, "G25": 95, "G50": 92,
                    "SEM": 1.3, "unit": "g/hen/d", "source": "Table 2"},
    "FCR_g_g": {"NC": 1.84, "PC": 1.79, "G25": 1.86, "G50": 1.79,
                 "SEM": 0.021, "unit": "g/g", "source": "Table 2"},
    "HD_pct": {"NC": 95.1, "PC": 93.4, "G25": 93.6, "G50": 95.4,
                "SEM": 0.85, "unit": "%", "source": "Table 2"},
    "EW_g": {"NC": 54.7, "PC": 54.6, "G25": 54.7, "G50": 53.9,
              "SEM": 0.58, "unit": "g", "source": "Table 2"},
    "EM_g": {"NC": 52.1, "PC": 51.0, "G25": 51.2, "G50": 51.4,
              "SEM": 0.88, "unit": "g", "source": "Table 2"},
    
    # Egg Quality (Table 3) — replicate 당 10개 계란 평균
    "Eggshell_L_star": {"NC": 50.1, "PC": 50.2, "G25": 50.6, "G50": 51.3,
                         "SEM": 0.40, "unit": "CIE L*", "source": "Table 3"},
    "Egg_yolk_color": {"NC": 8.5, "PC": 8.4, "G25": 8.5, "G50": 8.3,
                        "SEM": 0.07, "unit": "Roche fan", "source": "Table 3"},
    "Eggshell_strength": {"NC": 4.70, "PC": 4.68, "G25": 4.66, "G50": 4.67,
                           "SEM": 0.122, "unit": "kg/cm²", "source": "Table 3"},
    "Haugh_unit": {"NC": 91.0, "PC": 91.1, "G25": 91.0, "G50": 90.5,
                    "SEM": 0.57, "unit": "HU", "source": "Table 3"},
}

# ---- C. Diet 정보 (고정값: treatment별) ----
diet_info = {
    "SID_GlySer_pct": {"NC": 1.365, "PC": 1.365, "G25": 1.706, "G50": 2.048},
    "AMEn_kcal_kg": {"NC": 2725, "PC": 2825, "G25": 2731, "G50": 2736},
    "Analyzed_Gly_pct": {"NC": 0.645, "PC": 0.635, "G25": 1.003, "G50": 1.313},
    "Gly_supplement_pct": {"NC": 0.0, "PC": 0.0, "G25": 0.341, "G50": 0.683},
    "CP_analyzed_pct": {"NC": 16.62, "PC": 16.53, "G25": 17.09, "G50": 17.80},
}

# ===========================================================================
# 3) 통합 Metadata 구축
# ===========================================================================
print("\n" + "=" * 70)
print("[3] 통합 Metadata 구축...")

rows = []
for sid in meta_extended.index:
    treat = str(meta_extended.loc[sid, "Treat"])
    rep = meta_extended.loc[sid, "Replicate"]
    
    row = {
        "SampleID": sid,
        "Treat": treat,
        "Replicate": rep,
        "Tissue": "Jejunum_Content",
    }
    
    # Diet info
    for k, v in diet_info.items():
        row[k] = v[treat]
    
    # Individual phenotypes (group mean assigned)
    for k, v in individual_phenotypes.items():
        row[k] = v[treat]
        row[f"{k}_SEM"] = v["SEM"]
    
    # Replicate phenotypes (group mean assigned)
    for k, v in replicate_phenotypes.items():
        row[k] = v[treat]
        row[f"{k}_SEM"] = v["SEM"]
    
    rows.append(row)

integrated_meta = pd.DataFrame(rows)
integrated_meta = integrated_meta.set_index("SampleID")

# Column 분류 태그 추가용 별도 테이블
col_categories = {}
for k in diet_info:
    col_categories[k] = "Diet"
for k in individual_phenotypes:
    col_categories[k] = "Individual_phenotype"
for k in replicate_phenotypes:
    col_categories[k] = "Replicate_phenotype"

print(f"  통합 metadata: {integrated_meta.shape[0]} samples × {integrated_meta.shape[1]} columns")
print(f"  Diet variables: {len(diet_info)}")
print(f"  Individual phenotypes: {len(individual_phenotypes)}")
print(f"  Replicate phenotypes: {len(replicate_phenotypes)}")

# ===========================================================================
# 4) Phenotype Summary Table 생성
# ===========================================================================
print("\n" + "=" * 70)
print("[4] Phenotype Summary 정리...")

# A. 전체 phenotype summary
summary_rows = []
all_pheno = {**individual_phenotypes, **replicate_phenotypes}
for pname, pdata in all_pheno.items():
    row = {
        "Phenotype": pname,
        "Category": "Individual" if pname in individual_phenotypes else "Replicate",
        "Unit": pdata["unit"],
        "Source": pdata["source"],
    }
    for t in TREAT_ORDER:
        row[f"Mean_{t}"] = pdata[t]
    row["SEM"] = pdata["SEM"]
    
    # Direction: NC vs others
    row["PC_vs_NC_direction"] = "↑" if pdata["PC"] > pdata["NC"] else "↓" if pdata["PC"] < pdata["NC"] else "="
    row["G50_vs_NC_direction"] = "↑" if pdata["G50"] > pdata["NC"] else "↓" if pdata["G50"] < pdata["NC"] else "="
    row["G50_mimics_PC"] = "Yes" if (row["PC_vs_NC_direction"] == row["G50_vs_NC_direction"]) else "No"
    
    summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(TAB_DIR, "phenotype_summary_all.csv"), index=False)

print("\n  Phenotype Direction Analysis (G50 mimics PC?):")
for _, r in summary_df.iterrows():
    mimic_flag = "✓" if r["G50_mimics_PC"] == "Yes" else "✗"
    print(f"    {mimic_flag} {r['Phenotype']}: PC{r['PC_vs_NC_direction']}NC, G50{r['G50_vs_NC_direction']}NC")

# ===========================================================================
# 5) Phenotype Heatmap (z-score normalized)
# ===========================================================================
print("\n" + "=" * 70)
print("[5] Phenotype Heatmap 생성...")

# Select key phenotypes for visualization
key_phenos_individual = [
    "AST_UL", "ALT_UL", "Uric_acid_mg_dL", "Creatinine_mg_dL",
    "Blood_HL_ratio", "Feather_CORT_ng_g",
    "Liver_color_score", "Liver_L_star", "Liver_b_star",
    "Hemorrhagic_score", "Liver_fat_pct_DM"
]

key_phenos_replicate = [
    "FI_g_hen_d", "FCR_g_g", "HD_pct",
    "Eggshell_L_star", "Egg_yolk_color"
]

# Build group mean matrix
pheno_mean = pd.DataFrame()
for pname in key_phenos_individual + key_phenos_replicate:
    all_p = {**individual_phenotypes, **replicate_phenotypes}
    for t in TREAT_ORDER:
        pheno_mean.loc[pname, t] = all_p[pname][t]

# Z-score normalize rows
pheno_z = pheno_mean.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# Pretty labels
label_map = {
    "AST_UL": "AST (U/L)", "ALT_UL": "ALT (U/L)",
    "Uric_acid_mg_dL": "Uric acid (mg/dL)", "Creatinine_mg_dL": "Creatinine (mg/dL)",
    "Blood_HL_ratio": "Blood H:L ratio", "Feather_CORT_ng_g": "Feather CORT (ng/g)",
    "Liver_color_score": "Liver color score", "Liver_L_star": "Liver L*",
    "Liver_b_star": "Liver b*", "Hemorrhagic_score": "Hemorrhagic score",
    "Liver_fat_pct_DM": "Liver fat (% DM)",
    "FI_g_hen_d": "Feed intake (g/d)", "FCR_g_g": "FCR (g/g)",
    "HD_pct": "Hen-day prod (%)",
    "Eggshell_L_star": "Eggshell L*", "Egg_yolk_color": "Yolk color"
}

pheno_z.index = [label_map.get(x, x) for x in pheno_z.index]

# Category color bar
category_colors = []
for p in key_phenos_individual + key_phenos_replicate:
    if p in individual_phenotypes:
        category_colors.append("#E64B35")  # Individual
    else:
        category_colors.append("#4DBBD5")  # Replicate

fig, ax = plt.subplots(figsize=(7, 9))
sns.heatmap(pheno_z, cmap='RdBu_r', center=0, annot=pheno_mean.values,
            fmt='.2f', linewidths=0.8, linecolor='white',
            xticklabels=TREAT_ORDER, ax=ax,
            cbar_kws={'label': 'Z-score', 'shrink': 0.7})

# Add category annotation on left
for i, c in enumerate(category_colors):
    ax.add_patch(plt.Rectangle((-0.35, i), 0.25, 1, 
                               facecolor=c, edgecolor='none', 
                               clip_on=False, transform=ax.get_yaxis_transform()))

ax.set_title("Phenotype Profile by Treatment\n(Color: z-score, Values: group mean)",
             fontweight='bold', fontsize=12, pad=20)
ax.set_ylabel("")

# Legend for category
legend_patches = [
    mpatches.Patch(facecolor="#E64B35", label="Individual-level"),
    mpatches.Patch(facecolor="#4DBBD5", label="Replicate-level"),
]
ax.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(0, -0.05),
          ncol=2, fontsize=9, frameon=False)

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig9_Phenotype_Heatmap.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig9_Phenotype_Heatmap.pdf"))
plt.close()
print("  Phenotype heatmap saved.")

# ===========================================================================
# 6) Microbiome–Phenotype Correlation (Spearman, group-level)
# ===========================================================================
print("\n" + "=" * 70)
print("[6] Microbiome–Phenotype Correlation...")

# Aggregate genus-level relative abundance
def aggregate_genus(feat_rel, tax_df):
    tax_genus = tax_df["Genus"].fillna("Unclassified").copy()
    agg = feat_rel.copy()
    agg['taxon'] = tax_genus.loc[agg.index].values
    agg = agg.groupby('taxon').sum()
    # Filter: >0.5% in at least one sample
    agg = agg.loc[agg.max(axis=1) >= 0.005]
    return agg

genus_agg = aggregate_genus(feat_rel, tax_df)

# Phenotype values per sample (group mean assigned)
pheno_cols_for_corr = list(individual_phenotypes.keys()) + list(replicate_phenotypes.keys())

corr_results = []
for taxon in genus_agg.index:
    if taxon in ["Others", "Unclassified"]:
        continue
    for pheno in pheno_cols_for_corr:
        tax_vals = genus_agg.loc[taxon, integrated_meta.index].values
        pheno_vals = integrated_meta.loc[genus_agg.columns, pheno].values
        
        # Spearman (note: within-group variation in phenotype is 0 since group means)
        # So we're effectively correlating taxon abundance with treatment-level phenotype
        rho, p = stats.spearmanr(tax_vals, pheno_vals)
        
        corr_results.append({
            "Genus": taxon, "Phenotype": pheno,
            "Spearman_rho": rho, "Spearman_p": p
        })

corr_df = pd.DataFrame(corr_results)
_, corr_df["BH_p"], _, _ = smm.multipletests(corr_df["Spearman_p"].fillna(1), method='fdr_bh')
corr_df = corr_df.sort_values("Spearman_p")
corr_df.to_csv(os.path.join(TAB_DIR, "microbiome_phenotype_correlation.csv"), index=False)

# Significant correlations
sig_corr = corr_df[corr_df["Spearman_p"] < 0.05]
print(f"  Significant correlations (raw p < 0.05): {len(sig_corr)}")
for _, r in sig_corr.head(15).iterrows():
    print(f"    {r['Genus']} × {r['Phenotype']}: ρ={r['Spearman_rho']:.3f}, p={r['Spearman_p']:.4f}")

# ===========================================================================
# 7) Correlation Heatmap — key taxa × key phenotypes
# ===========================================================================
print("\n" + "=" * 70)
print("[7] Correlation Heatmap 생성...")

# Select taxa with at least one significant correlation
sig_taxa = sig_corr["Genus"].unique()[:15]

# Select key phenotypes (biologically meaningful for gut-liver axis)
key_phenos = [
    "AST_UL", "ALT_UL", "Uric_acid_mg_dL",
    "Blood_HL_ratio", "Feather_CORT_ng_g",
    "Liver_color_score", "Liver_fat_pct_DM", "Hemorrhagic_score",
    "FI_g_hen_d", "FCR_g_g"
]

key_labels = {
    "AST_UL": "AST", "ALT_UL": "ALT",
    "Uric_acid_mg_dL": "Uric acid",
    "Blood_HL_ratio": "H:L ratio", "Feather_CORT_ng_g": "CORT",
    "Liver_color_score": "Liver score", "Liver_fat_pct_DM": "Liver fat",
    "Hemorrhagic_score": "Hemorrhage",
    "FI_g_hen_d": "Feed intake", "FCR_g_g": "FCR"
}

# Build correlation matrix
rho_mat = pd.DataFrame(index=sig_taxa, columns=key_phenos, dtype=float)
p_mat = pd.DataFrame(index=sig_taxa, columns=key_phenos, dtype=float)

for taxon in sig_taxa:
    for pheno in key_phenos:
        row = corr_df[(corr_df["Genus"] == taxon) & (corr_df["Phenotype"] == pheno)]
        if len(row) > 0:
            rho_mat.loc[taxon, pheno] = row.iloc[0]["Spearman_rho"]
            p_mat.loc[taxon, pheno] = row.iloc[0]["Spearman_p"]

rho_mat = rho_mat.astype(float)
p_mat = p_mat.astype(float)

# Annotation: rho + significance stars
annot = rho_mat.copy().astype(str)
for i in rho_mat.index:
    for j in rho_mat.columns:
        r = rho_mat.loc[i, j]
        p = p_mat.loc[i, j]
        stars = ""
        if p < 0.001:
            stars = "***"
        elif p < 0.01:
            stars = "**"
        elif p < 0.05:
            stars = "*"
        annot.loc[i, j] = f"{r:.2f}{stars}"

fig, ax = plt.subplots(figsize=(12, max(6, len(sig_taxa) * 0.5)))
vmax = max(abs(rho_mat.values.min()), abs(rho_mat.values.max()))

sns.heatmap(rho_mat, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
            annot=annot, fmt='', linewidths=0.6, linecolor='white',
            xticklabels=[key_labels.get(k, k) for k in key_phenos],
            ax=ax, cbar_kws={'label': 'Spearman ρ', 'shrink': 0.7})

ax.set_title("Microbiome–Phenotype Correlation\n(Spearman ρ, * p<0.05, ** p<0.01, *** p<0.001)",
             fontweight='bold', fontsize=12)
ax.set_ylabel("Genus")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig10_Microbiome_Phenotype_Correlation.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig10_Microbiome_Phenotype_Correlation.pdf"))
plt.close()
print("  Correlation heatmap saved.")

# ===========================================================================
# 8) Mantel Test — Microbiome distance vs Phenotype distance
# ===========================================================================
print("\n" + "=" * 70)
print("[8] Mantel Test: Microbiome ↔ Phenotype...")

from skbio.stats.distance import mantel

# Microbiome distance matrices
otu_mat = feat.T.values
dm_bc = beta_diversity('braycurtis', otu_mat, ids=sample_ids)
dm_wu = beta_diversity('weighted_unifrac', otu_mat, ids=sample_ids,
                        tree=tree, taxa=list(feat.index))

# Phenotype distance matrix (Euclidean on z-scored phenotypes)
pheno_for_mantel = integrated_meta.loc[sample_ids, list(individual_phenotypes.keys())].astype(float)
# Z-score
pheno_z_samples = (pheno_for_mantel - pheno_for_mantel.mean()) / pheno_for_mantel.std()
pheno_dist_vals = pdist(pheno_z_samples.values, metric='euclidean')
dm_pheno = DistanceMatrix(squareform(pheno_dist_vals), ids=sample_ids)

# Subsets of phenotypes
pheno_subsets = {
    "All_individual": list(individual_phenotypes.keys()),
    "Liver_only": ["Liver_color_score", "Liver_L_star", "Liver_b_star",
                    "Hemorrhagic_score", "Liver_fat_pct_DM"],
    "Stress_only": ["Blood_HL_ratio", "Feather_CORT_ng_g"],
    "Serum_only": ["AST_UL", "ALT_UL", "Uric_acid_mg_dL", "Creatinine_mg_dL"],
}

mantel_results = []
for micro_name, dm_micro in [("Bray-Curtis", dm_bc), ("Weighted_UniFrac", dm_wu)]:
    for pheno_name, pheno_cols in pheno_subsets.items():
        pheno_sub = integrated_meta.loc[sample_ids, pheno_cols].astype(float)
        pheno_sub_z = (pheno_sub - pheno_sub.mean()) / pheno_sub.std()
        # Handle constant columns
        pheno_sub_z = pheno_sub_z.fillna(0)
        
        pheno_dist = pdist(pheno_sub_z.values, metric='euclidean')
        dm_pheno_sub = DistanceMatrix(squareform(pheno_dist), ids=sample_ids)
        
        r, p, n = mantel(dm_micro, dm_pheno_sub, method='spearman', permutations=9999)
        
        mantel_results.append({
            "Microbiome_distance": micro_name,
            "Phenotype_set": pheno_name,
            "Mantel_r": r, "p_value": p, "n_permutations": n
        })
        
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"  {micro_name} × {pheno_name}: r={r:.3f}, p={p:.4f} {sig}")

mantel_df = pd.DataFrame(mantel_results)
mantel_df.to_csv(os.path.join(TAB_DIR, "mantel_test_results.csv"), index=False)

# --- Mantel test visualization ---
fig, ax = plt.subplots(figsize=(8, 4.5))

mantel_pivot = mantel_df.pivot(index="Phenotype_set", columns="Microbiome_distance", values="Mantel_r")
mantel_p_pivot = mantel_df.pivot(index="Phenotype_set", columns="Microbiome_distance", values="p_value")

# Annotation with p-values
annot_mantel = mantel_pivot.copy().astype(str)
for i in mantel_pivot.index:
    for j in mantel_pivot.columns:
        r = mantel_pivot.loc[i, j]
        p = mantel_p_pivot.loc[i, j]
        stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        annot_mantel.loc[i, j] = f"r={r:.3f}\np={p:.3f}{stars}"

sns.heatmap(mantel_pivot, cmap='YlOrRd', annot=annot_mantel, fmt='',
            linewidths=0.8, linecolor='white', vmin=0, vmax=0.5,
            ax=ax, cbar_kws={'label': 'Mantel r', 'shrink': 0.8})
ax.set_title("Mantel Test: Microbiome Distance vs Phenotype Distance\n(* p<0.05, ** p<0.01, *** p<0.001)",
             fontweight='bold')
ax.set_ylabel("Phenotype Set")
ax.set_xlabel("Microbiome Distance")
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "Fig11_Mantel_Test.png"), dpi=300)
fig.savefig(os.path.join(FIG_DIR, "Fig11_Mantel_Test.pdf"))
plt.close()
print("  Mantel test figure saved.")

# ===========================================================================
# 9) Key Taxa × Phenotype Scatter (Top associations)
# ===========================================================================
print("\n" + "=" * 70)
print("[9] Key association scatter plots...")

# Top 6 significant correlations (biologically meaningful)
top_assoc = sig_corr[
    sig_corr["Phenotype"].isin(["Liver_fat_pct_DM", "Blood_HL_ratio", 
                                 "Liver_color_score", "AST_UL", "ALT_UL",
                                 "FCR_g_g", "FI_g_hen_d", "Hemorrhagic_score"])
].head(6)

if len(top_assoc) > 0:
    n_plots = min(6, len(top_assoc))
    ncols = 3
    nrows = (n_plots + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    if nrows == 1:
        axes = [axes] if n_plots == 1 else axes
    axes = np.array(axes).flatten()
    
    for idx, (_, row) in enumerate(top_assoc.iterrows()):
        if idx >= n_plots:
            break
        ax = axes[idx]
        taxon = row["Genus"]
        pheno = row["Phenotype"]
        
        for t in TREAT_ORDER:
            samps = integrated_meta[integrated_meta["Treat"] == t].index
            samps = [s for s in samps if s in genus_agg.columns]
            tax_vals = genus_agg.loc[taxon, samps].values
            pheno_vals = integrated_meta.loc[samps, pheno].values
            
            ax.scatter(tax_vals, pheno_vals, color=TREAT_COLORS[t],
                      s=60, label=t, edgecolor='white', linewidth=0.5, alpha=0.9, zorder=5)
        
        # Trend line
        all_tax = genus_agg.loc[taxon, integrated_meta.index].values
        all_pheno = integrated_meta[pheno].values
        z = np.polyfit(all_tax, all_pheno, 1)
        xrange = np.linspace(all_tax.min(), all_tax.max(), 100)
        ax.plot(xrange, np.polyval(z, xrange), 'k--', alpha=0.4, linewidth=1)
        
        p_str = f"p={row['Spearman_p']:.3f}" if row['Spearman_p'] >= 0.001 else "p<0.001"
        ax.set_title(f"{taxon} × {label_map.get(pheno, pheno)}\nρ={row['Spearman_rho']:.2f}, {p_str}",
                    fontsize=10, fontweight='bold')
        ax.set_xlabel(f"{taxon} (Rel. Abund.)")
        ax.set_ylabel(label_map.get(pheno, pheno))
        ax.legend(fontsize=8, frameon=True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    for j in range(idx + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle("Microbiome–Phenotype Associations (Top Correlations)",
                 fontweight='bold', y=1.02)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "Fig12_Taxa_Phenotype_Scatter.png"), dpi=300)
    fig.savefig(os.path.join(FIG_DIR, "Fig12_Taxa_Phenotype_Scatter.pdf"))
    plt.close()
    print("  Scatter plots saved.")

# ===========================================================================
# 10) Integrated Metadata Export (Excel)
# ===========================================================================
print("\n" + "=" * 70)
print("[10] 통합 Metadata Excel 저장...")

integrated_meta.to_csv(os.path.join(TAB_DIR, "integrated_metadata_full.csv"))

# Separate clean version (without SEM columns)
clean_cols = [c for c in integrated_meta.columns if not c.endswith("_SEM")]
integrated_clean = integrated_meta[clean_cols]
integrated_clean.to_csv(os.path.join(TAB_DIR, "integrated_metadata_clean.csv"))

# ===========================================================================
# 11) Data Integration Summary — 분석 가능 범위 정리
# ===========================================================================
print("\n" + "=" * 70)
print("[11] 통합 분석 가능 범위 요약\n")

print("┌─────────────────────────────────────────────────────────────┐")
print("│  Data Integration Summary                                    │")
print("├─────────────────────────────────────────────────────────────┤")
print("│                                                               │")
print("│  [Level 1] 확실히 연결 가능 (Same Individual)                │")
print("│  ────────────────────────────────────────────────             │")
print("│  • Jejunal microbiome (16S) ←→ Serum (AST,ALT,UA,Cr)        │")
print("│  • Jejunal microbiome (16S) ←→ Blood H:L ratio              │")
print("│  • Jejunal microbiome (16S) ←→ Liver scores/fat             │")
print("│  • Jejunal microbiome (16S) ←→ Feather CORT                 │")
print("│  → 도축 시 동일 개체에서 채취                                 │")
print("│  → 현재: group mean 사용 (raw data 확보 시 individual 가능)  │")
print("│                                                               │")
print("│  [Level 2] 간접 연결 (Same Replicate)                        │")
print("│  ────────────────────────────────────────────────             │")
print("│  • Jejunal microbiome ←→ FI, FCR, HD, EM (12주 평균)        │")
print("│  • Jejunal microbiome ←→ Egg quality (실험 종료 시)          │")
print("│  → Replicate 단위 생산성, 도축 개체 = rep 대표               │")
print("│  → 통계적 독립성 주의 필요                                    │")
print("│                                                               │")
print("│  [Level 3] Treatment-level only                               │")
print("│  ────────────────────────────────────────────────             │")
print("│  • Diet composition (AMEn, Gly+Ser, CP)                      │")
print("│  → 처리구 고정값, dose variable로 사용                       │")
print("│                                                               │")
print("│  [추천] Raw individual data 확보 시 upgrade 가능:            │")
print("│  • Replicate 번호로 phenotype-microbiome 1:1 매칭            │")
print("│  • Within-group variation 활용한 진정한 correlation           │")
print("│  • Mixed-effects model 적용 가능                              │")
print("└─────────────────────────────────────────────────────────────┘")

print("\n[DONE] Phenotype-Microbiome Integration complete.")
