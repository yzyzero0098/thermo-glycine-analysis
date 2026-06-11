# Additional run package

This folder contains the **new follow-up outputs** generated after the original lab-meeting package was assembled.

## Included contents

### 1. Targeted host-pathway × taxa integration
Source directory:
- `experiments/glycine-targeted-integration/`

Copied files:
- `tables/module_taxa_correlations.csv`
- `tables/module_pc_convergence_summary.csv`
- `tables/gsva_scores_matched17.csv`
- `tables/mean_zscore_module_scores_matched17.csv`
- `tables/shortlist_taxa_relative_abundance_matched17.csv`
- `tables/module_gene_membership.csv`
- `tables/module_gene_sets.json`
- `tables/analysis_summary.json`
- `figures/gsva_module_scores_heatmap.png`
- `figures/module_taxa_correlation_heatmap.png`
- `figures/module_pc_convergence_ratios.png`
- `docs/targeted-integration-README.md`

### 2. PICRUSt2 pathway / EC follow-up
Source directory:
- `experiments/glycine-picrust-followup/`

Copied files:
- `tables/pathway_selected_features.csv`
- `tables/ec_selected_features.csv`
- `tables/picrust_analysis_summary.json`
- `figures/pathway_selected_heatmap.png`
- `figures/ec_selected_heatmap.png`
- `docs/pathway_summary.md`
- `docs/ec_summary.md`
- `docs/picrust-followup-README.md`

## What these files are for
- `module_taxa_correlations.csv`: main targeted integration result linking selected taxa with GSVA-based host modules
- `module_pc_convergence_summary.csv`: host-module convergence toward the PC reference state
- `pathway_selected_features.csv`: prioritized PICRUSt2 pathway features from saved pathway tables
- `ec_selected_features.csv`: prioritized PICRUSt2 EC features from saved EC tables

## Best files to open first
1. `tables/module_taxa_correlations.csv`
2. `tables/module_pc_convergence_summary.csv`
3. `tables/pathway_selected_features.csv`
4. `figures/module_taxa_correlation_heatmap.png`
5. `figures/module_pc_convergence_ratios.png`
6. `figures/pathway_selected_heatmap.png`

## Interpretation reminders
- GSVA-based host-module integration is a **candidate adaptive-axis analysis**, not causal proof.
- PICRUSt2 outputs are **predicted functions from 16S**, not directly measured microbial metabolism.
- These files extend the earlier lab-meeting package by moving beyond global concordance screening into targeted integration and predicted-function summaries.
