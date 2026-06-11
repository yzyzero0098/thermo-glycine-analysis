# Glycine 16S results inventory and figure-4–8 rebuild plan

## Scope
Inventory currently usable results under `/home/yzyzero/glycine_16s_analysis` and assess whether Figures 4–8 still fit the current filtering policy and manuscript flow.

## High-confidence usable result classes

### 1) Diversity / ordination
- **Alpha diversity**: `16s_preprocessing/core-metrics-results/*`, `Functional_analysis/Fig1_AlphaDiversity.png`, `Functional_analysis/alpha_diversity.csv`
- **Beta diversity**: original and extended PCoA outputs in `Functional_analysis/` including Jaccard and PERMDISP
- **Filtering-aware diversity input**: `16S_filtering/feature-table_filtered_offtarget_removed.tsv`

### 2) Taxonomy composition
- **Canonical figure input**: `16S_filtering/genus_top15_plus_others_relative_abundance.tsv`
- **Canonical off-target-cleaned ASV table** for phylum aggregation: `16S_filtering/feature-table_filtered_offtarget_removed.tsv`
- **Updated Figure 3**: `Functional_analysis/Fig3_Taxonomy_Composition_Combined.png`

### 3) Differential abundance / treatment-response taxa
- **Filtered DA tables**: `16S_filtering/taxon_level6_or_lowest_DA_primary_prev2_maxRA0.1pct.tsv`, `genus_table_DA_primary_prev2_maxRA0.1pct.tsv`, `genus_table_DA_sensitivity_prev4.tsv`
- **Cross-method DA summary**: `da_functional_phenotype_20260323/DA_cross_method_summary.md`
- **ANCOMBC tables**: `da_functional_phenotype_20260323/tables/ANCOMBC_*`
- **MaAsLin2 tables**: `da_functional_phenotype_20260323/maaslin2_*`

### 4) Functional prediction (PICRUSt2-derived)
- **Usable summary-level report**: `functional_analysis_20260324/functional_analysis_report.md`
- **Pathway / KO / EC tables**: `functional_analysis_20260324/tables/*`
- **Revised functional figures already exist**: `functional_analysis_20260324/figures_revised/*`

### 5) Integrative interpretation artifacts
- `glycine_16S_interpretation_report.md`
- `da_functional_phenotype_20260323/integration_status_report.md`

## Current figure-by-figure assessment

### Figure 4 — current: mean genus abundance heatmap
**Current file:** `Functional_analysis/Fig4_Genus_Heatmap.png`

**Assessment:** only partly appropriate now.
- Strength: easy to read and shows group-level shifts.
- Weakness 1: it is partly redundant with Figure 3, which already shows composition at sample level.
- Weakness 2: it comes from the older `Functional_analysis` workflow, not explicitly from the new locked filtering policy.
- Weakness 3: a mean-abundance heatmap is weaker than a treatment-response / DA-oriented heatmap if the paper's story is glycine mimicry and dose response.

**Recommendation:** rebuild.
- Preferred role: make Figure 4 a **response-focused genus heatmap** using filtered DA candidates rather than just mean abundance.
- Input: primary DA table (`taxon_level6_or_lowest...`) or a genus-level filtered summary derived from the same policy.

### Figure 5 — current: centroid distance to PC
**Current file:** `Functional_analysis/Fig5_Centroid_Distance_to_PC.png`

**Assessment:** conceptually appropriate and still useful.
- Strength: this is one of the clearest figures for the claim that **G50 is closest to PC**.
- Strength: aligns with the phenotype-mimicry narrative.
- Weakness: it should be recomputed from the current canonical diversity table / current beta-diversity workflow so all figures share the same preprocessing policy.
- Caveat: it is descriptive, not a formal significance test by itself.

**Recommendation:** keep the concept, rebuild the figure from the current off-target-cleaned diversity input and current Bray/weighted UniFrac outputs.

### Figure 6 — current: glycine-mimicking taxa heatmap
**Current file:** `Functional_analysis/Fig6_Mimicking_Taxa_Heatmap.png`

**Assessment:** not safe to keep as-is.
- Serious issue: the current figure includes **Chloroplast**, which conflicts with the now-locked off-target removal policy.
- It is likely based on older unfiltered taxon summaries and is therefore inconsistent with the current cleaned workflow.
- Some extreme log-fold changes appear to be driven by near-zero abundances, so ranking by smallest |ΔlogFC| alone may overstate unstable taxa.

**Recommendation:** rebuild.
- Input should come from the off-target-cleaned and DA-policy-consistent table.
- Add a stability rule: only consider taxa that pass the DA primary filter or an explicit minimum abundance/prevalence threshold.
- This figure is still useful conceptually, but only after policy-consistent recomputation.

### Figure 7 — current: dose-response trend grid
**Current file:** `Functional_analysis/Fig7_DoseResponse_Trend.png`

**Assessment:** useful concept, but current version is not fully acceptable.
- Strength: dose-response is an important biological angle.
- Problem: current grid includes **Chloroplast**, again indicating the old workflow leaked off-target taxa into trend selection.
- Problem: several displayed taxa are weak or noisy and may not survive the now-locked filtering rules.

**Recommendation:** rebuild.
- Restrict candidates to taxa passing the primary DA filter or a dose-trend candidate list derived from the cleaned table.
- Keep only the strongest trend candidates; fewer, cleaner panels will be better.

### Figure 8 — current: phylum boxplots
**Current file:** `Functional_analysis/Fig8_Phylum_Boxplots.png`

**Assessment:** lowest priority and probably not a strong main-text figure.
- Strength: simple and easy to interpret.
- Weakness 1: mostly nonsignificant phylum-level comparisons.
- Weakness 2: redundant with Figure 3 phylum panel.
- Weakness 3: phylum-level patterns are not the most mechanistically informative part of this dataset.

**Recommendation:** demote or replace.
- If kept, remake from the off-target-cleaned canonical table.
- Better option: move to supplement, or replace this slot with a stronger functional/PICRUSt2 summary figure.

## Does the overall current figure flow make sense?

## Current flow (roughly)
1. Alpha diversity
2. Beta diversity
3. Taxonomy composition
4. Mean genus heatmap
5. PC-closeness centroid distance
6. Glycine-mimicking taxa
7. Dose-response taxa
8. Phylum boxplots

## Assessment
**Partly, but not fully.**

### What works
- Figures 1–3 form a coherent introduction: diversity -> ordination -> composition.
- Figure 5 is strong for the manuscript's central biological claim that **G50 approximates PC**.
- Figure 7 is conceptually valuable because glycine dose is a biologically interpretable axis.

### What does not work well now
- Figure 4 is too close to Figure 3 in function.
- Figures 6 and 7 are currently contaminated by older filtering choices (including off-target taxa such as chloroplast).
- Figure 8 is weak and likely redundant.

## Better revised flow
1. **Fig 1** Alpha diversity
2. **Fig 2** Beta diversity (extended version with Jaccard/PERMDISP as needed)
3. **Fig 3** Taxonomy composition (already rebuilt)
4. **Fig 4** Differential-abundance / treatment-response heatmap (rebuild)
5. **Fig 5** Distance to PC / mimicry summary (rebuild, same concept)
6. **Fig 6** Glycine-mimicking taxa heatmap (rebuild)
7. **Fig 7** Dose-response taxa (rebuild)
8. **Fig 8** Either supplement-level phylum boxplots **or** replace with a stronger PICRUSt2 functional summary

## Recommended rebuild plan for Figures 4–8

### Figure 4 — rebuild as DA/treatment-response heatmap
- **Use:** primary DA table + cross-method candidate taxa
- **Goal:** show the most interpretable taxa across NC, PC, G25, G50
- **Why:** stronger than a simple mean-abundance heatmap and less redundant with Figure 3

### Figure 5 — rebuild as PC-closeness figure
- **Use:** off-target-cleaned diversity table and current Bray/weighted UniFrac distance outputs
- **Goal:** preserve the strongest mimicry story

### Figure 6 — rebuild mimicking taxa from cleaned table
- **Use:** off-target-cleaned + DA-policy-consistent taxon table
- **Goal:** show taxa whose direction is closest to PC, but exclude unstable low-information artifacts

### Figure 7 — rebuild dose-response candidates
- **Use:** primary filtered genus/taxon table on NC/G25/G50 only
- **Goal:** cleaner, smaller panel set of dose-responsive taxa

### Figure 8 — rethink
- **Option A:** remake phylum boxplots and move to supplement
- **Option B:** replace with a PICRUSt2 functional summary figure using `functional_analysis_20260324/figures_revised/*`
- **Preferred:** Option B if the manuscript wants one more mechanistic layer

## Best currently usable non-figure results to carry forward
- Weighted UniFrac significance / abundance-weighted community shift
- G50 minimal distance to PC in centroid-style comparisons
- Recurrent taxa across methods: `Romboutsia`, `Turicibacter`, `Clostridia UCG-014`, `Enterococcus`, `Muribaculaceae`, `Bacteroides`, `Lachnospiraceae NK4A136 group`
- PICRUSt2-derived pathway/KO/EC tables as predicted-function hypothesis support (not direct metagenomics)

## Main caveats
- `glycine_16S_interpretation_report.md` contains useful narrative synthesis but is partly based on older figure generations and should not be treated as fully policy-updated without cross-checking against the rebuilt figures.
- `analysis_part2_taxa_DA.py` is still labeled as exploratory/legacy relative to the canonical filtering outputs in `16S_filtering/`.
- Functional prediction results should be framed explicitly as **PICRUSt2-predicted** functions, not measured metagenomes.

## Next recommended step
Rebuild Figures **4, 5, 6, 7**, and decide whether **Figure 8** should be kept as supplement or replaced by a stronger functional figure.
