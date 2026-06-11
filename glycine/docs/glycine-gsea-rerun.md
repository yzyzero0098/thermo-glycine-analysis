# Glycine RNA-seq GSEA rerun in the asset folder

Updated: 2026-05-13

## What was done
I reran preranked GSEA into the **lab-meeting asset folder** and generated both **result tables** and **plot files**.

### Execution environment
- Runner: `/home/yzyzero/miniforge3/envs/rnaseqfunc/bin/Rscript`
- Script: `file:///home/yzyzero/feynman/experiments/glycine-gsea-rerun/run_fgsea_with_plots.R`
- GSEA engine: `fgsea::fgseaMultilevel`
- Gene sets: `msigdbr` **Gallus gallus** `C5 / GO:BP`
- Ranking statistic: Wald `stat`

### Input contrasts used
- `NC_vs_G25` — archived DESeq2 full table
- `NC_vs_G50` — archived DESeq2 full table
- `NC_vs_PC` — archived DESeq2 full table
- `G50_vs_PC` — saved PyDESeq2 direct-contrast full table

### Important caveat on contrast naming
The archived DE tables are named `NC_vs_G25`, `NC_vs_G50`, and `NC_vs_PC`, but the project continuity notes have historically discussed them biologically as treatment-side contrasts relative to NC. The file naming and sign convention are therefore not perfectly transparent from the filename alone. For interpretation below, I rely on the **project's existing biological framing plus the direct `G50_vs_PC` rerun** rather than pretending the archived filenames are self-explanatory.

---

## Where the rerun outputs were saved

### Tables
Saved under:
- `outputs/glycine-labmeeting-package-assets/tables/rnaseq/`

Key files:
- `NC_vs_G25_fgsea_go_bp.csv`
- `NC_vs_G50_fgsea_go_bp.csv`
- `NC_vs_PC_fgsea_go_bp.csv`
- `G50_vs_PC_fgsea_go_bp.csv`
- `*_top30_sig.csv`
- `fgsea_go_bp_summary.csv`
- `fgsea_interpretation_summary.csv`
- `fgsea_rerun_README.md`

### Figures
Saved under:
- `outputs/glycine-labmeeting-package-assets/figures/rnaseq/`

New GSEA plots:
- `GSEA_summary_dotplot.png`
- `NC_vs_G25_GSEA_top_pathways.png`
- `NC_vs_G25_GSEA_enrichment_top_positive.png`
- `NC_vs_G25_GSEA_enrichment_top_negative.png`
- `NC_vs_G50_GSEA_top_pathways.png`
- `NC_vs_G50_GSEA_enrichment_top_positive.png`
- `NC_vs_G50_GSEA_enrichment_top_negative.png`
- `NC_vs_PC_GSEA_top_pathways.png`
- `NC_vs_PC_GSEA_enrichment_top_positive.png`
- `NC_vs_PC_GSEA_enrichment_top_negative.png`
- `G50_vs_PC_GSEA_top_pathways.png`
- `G50_vs_PC_GSEA_enrichment_top_positive.png`
- `G50_vs_PC_GSEA_enrichment_top_negative.png`

---

## Main rerun summary
From `fgsea_go_bp_summary.csv`:

| Contrast | Ranked genes | FDR < 0.05 pathways | Top positive pathway | Top negative pathway |
|---|---:|---:|---|---|
| NC_vs_G25 | 17,363 | 174 | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION | GOBP_GOLGI_VESICLE_TRANSPORT |
| NC_vs_G50 | 17,363 | 717 | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION | GOBP_GOLGI_VESICLE_TRANSPORT |
| NC_vs_PC | 17,363 | 275 | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION | GOBP_GOLGI_VESICLE_TRANSPORT |
| G50_vs_PC | 13,673 | 465 | GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION | GOBP_SMALL_MOLECULE_CATABOLIC_PROCESS |

---

## Interpretation

## 1. G50 still shows the strongest host pathway shift versus NC
Among the three archived NC-referenced contrasts, `NC_vs_G50` has the largest number of significant GO:BP sets:
- `NC_vs_G25`: **174**
- `NC_vs_G50`: **717**
- `NC_vs_PC`: **275**

This supports the same broad conclusion as the prior RNA-seq review:
> **G50 is the strongest transcriptomic remodeling state among the glycine groups.**

That does **not** mean every pathway in G50 simply becomes PC-like. It means the pathway-level host response is largest in G50.

---

## 2. The recurring biological themes are stable across the rerun
Across `NC_vs_G25`, `NC_vs_G50`, and `NC_vs_PC`, the rerun repeatedly highlights the same two broad pathway families.

### Recurring positive-side theme
- external encapsulating structure organization
- collagen fibril organization
- basement membrane organization
- cell-substrate adhesion
- epithelial / structural remodeling-related terms

### Recurring negative-side theme
- Golgi vesicle transport
- ER-to-Golgi vesicle-mediated transport
- vesicle organization
- endosomal transport
- macroautophagy / proteostasis-related terms
- in stronger contrasts, respiratory / OXPHOS-related terms

So the rerun is not generating a completely new story. It reinforces the existing pathway-centered interpretation:
> **glycine/PC-associated host remodeling involves extracellular/adhesion structure shifts together with trafficking, organelle, and metabolic reorganization.**

---

## 3. G25 is weak at DEG threshold but not null at pathway level
`NC_vs_G25` still has a much smaller pathway count than `NC_vs_G50`, but it is clearly not empty.

The top significant `NC_vs_G25` terms still include:
- extracellular encapsulating structure organization
- collagen fibril organization
- Golgi vesicle transport
- macroautophagy
- vesicle organization

So the safer interpretation remains:
> **G25 is weak by hard DEG count, but it is not transcriptomically null once the ranked-list pathway level is examined.**

---

## 4. G50 is not identical to PC
This rerun adds a direct `G50_vs_PC` fgsea result, which was missing from the previously copied result set.

`G50_vs_PC` still shows **465 significant GO:BP sets at FDR < 0.05**, which is too many to support a “G50 and PC are basically the same host state” claim.

The strongest `G50_vs_PC` negative-side pathways include:
- small molecule catabolic process
- organic acid catabolic process
- oxidative phosphorylation
- electron transport chain
- ATP biosynthetic process
- proton transmembrane transport

This is important because it suggests:
> **G50 moves toward PC on some adaptive axes, but the host transcriptome remains distinctly different from PC, especially in energy / catabolic / respiratory pathway structure.**

That is more defensible than saying “G50 fully mimics PC.”

---

## 5. The rerun is directionally consistent with the current GSVA module story
The current targeted integration outputs already emphasized modules such as:
- ECM / adhesion
- Golgi / vesicle transport
- autophagy / proteostasis
- respiration / OXPHOS

The rerun GSEA again surfaces exactly those axes. That is useful because it means the GSVA module design was **not arbitrary**; it matches the stronger contrast-level pathway signal in the rerun.

In practical terms:
- **GSEA** justifies the pathway/module choice,
- **GSVA** can still be used to score those modules sample-by-sample,
- and **WGCNA** remains the next unsupervised check if we want to see whether these same axes emerge without predefining the gene sets.

---

## 6. Why the rerun counts differ slightly from the old copied summary
The older copied summary in the workspace reported different significant-pathway counts for some contrasts.

This rerun produced:
- `NC_vs_G25`: 174 instead of the older 182
- `NC_vs_G50`: 717 instead of the older 813
- `NC_vs_PC`: 275 instead of the older 271

The most likely reason is that the rerun used the **current local R environment and current `msigdbr` package version**, so the exact GO set definitions and/or package behavior differ from the earlier saved run.

What survived the rerun despite those count differences is the more important point:
- `NC_vs_G50` remains strongest,
- ECM / adhesion terms remain prominent,
- Golgi / vesicle / trafficking terms remain prominent,
- and direct `G50_vs_PC` still shows substantial residual host-state difference.

---

## Bottom line
If I had to summarize the rerun in one sentence:

> **The rerun GSEA confirms that G50 produces the strongest host transcriptomic pathway remodeling, that G25 retains non-null pathway signal despite weak DEG counts, and that G50 is only partially PC-like because substantial energy, catabolic, and organelle-trafficking pathway differences remain in the direct G50_vs_PC comparison.**

---

## Recommended use in slides
For a lab meeting deck, I would use:
1. `GSEA_summary_dotplot.png`
2. `NC_vs_G50_GSEA_top_pathways.png`
3. `G50_vs_PC_GSEA_top_pathways.png`
4. if needed, one enrichment plot each for:
   - `NC_vs_G50_GSEA_enrichment_top_positive.png`
   - `G50_vs_PC_GSEA_enrichment_top_negative.png`

That gives both:
- the **largest host-remodeling contrast** (`NC_vs_G50`), and
- the **non-equivalence contrast** (`G50_vs_PC`).

---

## Open questions / next steps
- If needed, fix the archived contrast-sign labeling more explicitly so positive/negative side interpretation is unambiguous in the saved tables.
- Condense the repeated GO terms into 4–6 narrative modules for presentation.
- Reuse those condensed modules in GSVA and then compare them against future WGCNA module eigengenes.

---

## Sources
- Script: `file:///home/yzyzero/feynman/experiments/glycine-gsea-rerun/run_fgsea_with_plots.R`
- Summary table: `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/tables/rnaseq/fgsea_go_bp_summary.csv`
- Interpretation table: `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/tables/rnaseq/fgsea_interpretation_summary.csv`
- Full rerun table: `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/tables/rnaseq/G50_vs_PC_fgsea_go_bp.csv`
- Existing GSVA/module context: `file:///home/yzyzero/feynman/experiments/glycine-targeted-integration/module_pc_convergence_summary.csv`
- Existing GSVA/module context: `file:///home/yzyzero/feynman/experiments/glycine-targeted-integration/gsva_scores_matched17.csv`
