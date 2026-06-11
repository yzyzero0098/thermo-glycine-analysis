# Glycine lab-meeting asset folder guide

This folder is a curated copy of the most reusable glycine project assets for lab meeting, slide drafting, and manuscript planning.

## Folder structure
- `docs/` — previously written markdown summaries and interpretation notes
- `figures/16s/` — current microbiome figures
- `figures/rnaseq/` — current host transcriptome figures
- `figures/integration/` — current integration summary figures
- `tables/16s/` — current microbiome result tables
- `tables/rnaseq/` — current host transcriptome result tables
- `tables/integration/` — current matching, concordance, and phenotype-limitation tables

## Recommended first-open files
### Narrative / context
1. `docs/glycine-intro-animal-nutrition.md`
2. `docs/glycine-omics-reanalysis.md`
3. `docs/glycine-log-interpretation.md`
4. `docs/glycine-procrustes-integration.md`
5. `docs/glycine-results-inventory.md`

### Fast figure set for slides
1. `figures/16s/Fig3_Taxonomy_Composition_Combined.png`
2. `figures/16s/Fig5_Centroid_Distance_to_PC_Canonical.png`
3. `figures/rnaseq/RNAseq_MDS_edgeR_panels.png`
4. `figures/rnaseq/NC_vs_G50_GO_BP_dotplot.png`
5. `figures/integration/glycine-procrustes-integration-pc-distance-ratios.png`

### Fast number tables for discussion
1. `tables/16s/BetaDiversity_global_summary.csv`
2. `tables/rnaseq/DEG_method_comparison.csv`
3. `tables/rnaseq/deseq2_threshold_sensitivity_counts.csv`
4. `tables/rnaseq/G50_vs_PC_PyDESeq2_summary.csv`
5. `tables/integration/integration_summary.csv`
6. `tables/integration/centroid_procrustes_summary.csv`
7. `tables/integration/direct_pc_closeness_summary.csv`

## Use notes by folder

## `docs/`
- `glycine-intro-animal-nutrition.md` — current introduction framing packet for Animal Nutrition / Poultry Science style
- `glycine-results-inventory.md` — figure-by-figure 16S inventory and rebuild notes
- `glycine-omics-reanalysis.md` — RNA-seq + integration feasibility reanalysis memo
- `glycine-gsea-followup.md` — earlier pathway-centered follow-up logic
- `glycine-gsea-rerun.md` — rerun fgsea/GSEA summary with new plots and direct `G50_vs_PC` interpretation
- `glycine-gsea-slide-notes.md` — recommended GSEA figures plus slide-ready English captions and interpretations
- `glycine-log-interpretation.md` — current biological interpretation memo
- `glycine-procrustes-integration.md` — strongest current integration writeup
- `glycine-procrustes-integration.provenance.md` — provenance sidecar for the integration memo

## `figures/16s/`
- `Fig1_AlphaDiversity.png` — alpha diversity overview
- `Fig2_BetaDiversity_PCoA_PERMDISP_Jaccard.png` — extended beta-diversity summary
- `Fig3_Taxonomy_Composition_Combined.png` — best current composition figure
- `Fig4_Treatment_Response_Heatmap.png` — response-focused taxa heatmap
- `Fig5_Centroid_Distance_to_PC_Canonical.png` — PC-reference proximity summary
- `Fig6_Glycine_Mimicking_Taxa_Canonical.png` — partial PC-like taxa figure
- `Fig7_Dose_Response_Canonical.png` — dose-responsive taxa figure
- `Fig9_PICRUSt2_Pathway_Heatmap.png` — predicted pathway heatmap
- `Fig10_PICRUSt2_Pathway_Dose_Mimicry.png` — predicted pathway dose/PC-reference figure
- `ANCOMBC_genus_ref_NC_top20.png` / `ANCOMBC_genus_ref_PC_top20.png` / `MaAsLin2_treatment_top20.png` — DA support figures

## `figures/rnaseq/`
- `Fig1_RNAseq_PCA.png` — saved PCA view
- `RNAseq_MDS_edgeR_panels.png` — most useful group-structure figure
- `NC_vs_PC_volcano.png`, `NC_vs_G25_volcano.png`, `NC_vs_G50_volcano.png` — contrast volcanoes
- `NC_vs_G50_heatmap_top40.png` — strong host-response heatmap
- `NC_vs_G50_GO_BP_dotplot.png` / `NC_vs_G50_KEGG_dotplot.png` — older pathway summary plots
- `GSEA_summary_dotplot.png` — rerun cross-contrast fgsea overview
- `NC_vs_G50_GSEA_top_pathways.png` — strongest archived contrast in the rerun
- `G50_vs_PC_GSEA_top_pathways.png` — direct non-equivalence contrast in the rerun
- `RNAseq_MDS_outlier_candidates.png`, `DEG_stability_heatmap_top20.png`, `MDS_comparison_top_scenarios.png` — outlier/sensitivity backup

## `figures/integration/`
- `glycine-procrustes-integration-pc-distance-ratios.png` — best current integration summary figure
- `integration_procrustes.png` — exploratory whole-ordination summary figure
- `phenotype_summary.png` — treatment-profile phenotype overview
- `microbiome_beta_diversity_bray_pcoa.png` — Bray PCoA used in old multiomics pipeline
- `rnaseq_mds.png` — RNA MDS used in old multiomics pipeline

## `tables/16s/`
Main discussion tables:
- `alpha_diversity.csv`
- `BetaDiversity_global_summary.csv`
- `PERMANOVA_global_extended.csv`
- `PERMDISP_global.csv`
- `taxonomy_composition_phylum.csv`
- `taxonomy_composition_genus.csv`
- `Fig4_taxon_kw_ranked.csv`
- `Fig5_centroid_distance_to_PC.csv`
- `Fig6_glycine_mimic_taxa_ranked.csv`
- `Fig7_dose_response_ranked.csv`
- `ANCOMBC_genus_global.csv`
- `ANCOMBC_genus_glydose_linear.csv`
- `Pathway_kw_annotated.csv`
- `Pathway_dose_trend_annotated.csv`
- `Pathway_pc_mimicry_annotated.csv`
- `NSTI_summary.csv`

Input-like tables retained because they help figure interpretation:
- `genus_top15_plus_others_relative_abundance.tsv`
- `genus_table_DA_primary_prev2_maxRA0.1pct.tsv`
- `taxon_level6_or_lowest_DA_primary_prev2_maxRA0.1pct.tsv`

## `tables/rnaseq/`
Main discussion tables:
- `DEG_method_comparison.csv`
- `DEG_summary.csv`
- `deseq2_threshold_sensitivity_counts.csv`
- `deseq2_near_threshold_genes.csv`
- `mds_group_dispersion_summary.csv`
- `mds_centroid_distances.csv`
- `fgsea_go_bp_summary.csv`
- `NC_vs_PC_fgsea_go_bp_top30_sig.csv`
- `NC_vs_G25_fgsea_go_bp_top30_sig.csv`
- `NC_vs_G50_fgsea_go_bp_top30_sig.csv`
- `G50_vs_PC_PyDESeq2_summary.csv`
- `G50_vs_PC_PyDESeq2_sig_fdr005_lfc1.csv`

Structure/QC support:
- `sample_metadata.csv`
- `RNAseq_MDS_edgeR_panel_coordinates.csv`
- `MDS_coordinates_and_outlier_candidates.csv`
- `outlier_key_scenarios.csv`

## `tables/integration/`
Main discussion tables:
- `integration_summary.csv` — exploratory Mantel/Procrustes summary
- `centroid_procrustes_summary.csv` — stronger centroid-level reanalysis summary
- `direct_pc_closeness_summary.csv` — best current convergence table
- `matched_group_counts.csv` — matched-17 group counts
- `matched_host_microbiome_retained16_counts.csv` — retained exact-pair count summary for later Procrustes framing
- `phenotype_repeat_pattern_summary.csv` — shows why phenotype is not a defensible sample-level block
- `phenotype_pairwise_standardized_distances.csv` — treatment-profile phenotype distance summary

Supporting coordinate/distance tables:
- `host_pca_full25_variance.csv`
- `host_pca_full25_coordinates.csv`
- `bray_pcoa_coordinates.csv`
- `weighted_unifrac_pcoa_coordinates.csv`
- `host_full25_centroid_pairwise_distances.csv`
- `host_matched16_centroid_pairwise_distances.csv`
- `bray_centroid_pairwise_distances.csv`
- `weighted_unifrac_centroid_pairwise_distances.csv`
- `sample_procrustes_summary.csv`
- `matched_17_samples_full_metadata.csv`
- `matched_17_samples_integration_ready.csv`
- `matched_host_microbiome_retained16.csv`

## Interpretation rules to remember
1. Do **not** sell Mantel/Procrustes as proof of strong global concordance.
2. Do **not** treat phenotype as an independent sample-level block unless raw per-bird values are recovered.
3. Do **not** describe PICRUSt2 as direct metagenomics/metabolomics.
4. Do use this folder to support a **shared-heat-stress adaptive remodeling** story.
