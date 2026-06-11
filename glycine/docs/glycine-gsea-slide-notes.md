# GSEA slide notes for the glycine lab-meeting deck

Updated: 2026-05-13

## Recommended main figure set
If you want the cleanest 2-slide GSEA story, use these figures in this order:

### Slide 1 — Overall pathway overview
- `figures/rnaseq/GSEA_summary_dotplot.png`

### Slide 2 — Key biological contrast and non-equivalence check
- `figures/rnaseq/NC_vs_G50_GSEA_top_pathways.png`
- `figures/rnaseq/G50_vs_PC_GSEA_top_pathways.png`

### Optional backup / appendix
- `figures/rnaseq/NC_vs_G50_GSEA_enrichment_top_positive.png`
- `figures/rnaseq/NC_vs_G50_GSEA_enrichment_top_negative.png`
- `figures/rnaseq/G50_vs_PC_GSEA_enrichment_top_negative.png`

---

## Why these are the best choices

### 1. `GSEA_summary_dotplot.png`
This is the best overview figure because it shows:
- all major contrasts together,
- relative pathway significance,
- pathway direction by NES,
- and which biology repeats across contrasts.

### 2. `NC_vs_G50_GSEA_top_pathways.png`
This is the best main contrast figure because `NC_vs_G50` is the strongest archived host transcriptome contrast at the pathway level.

### 3. `G50_vs_PC_GSEA_top_pathways.png`
This is the most important follow-up figure because it prevents overclaiming. It shows that `G50` is **not identical to `PC`**, even if some adaptive directions overlap.

---

## Slide-ready English captions and one-line interpretations

## Figure: `GSEA_summary_dotplot.png`
### Caption
**Cross-contrast preranked GSEA summary of jejunal host transcriptome remodeling under glycine supplementation and PC reference feeding.** Dot size indicates pathway significance (−log10 FDR), and color indicates the normalized enrichment score (NES).

### One-line interpretation
**The pathway-level host response was strongest in the G50 contrast, while recurring enrichment themes across contrasts included extracellular/adhesion remodeling and trafficking-organelle-metabolic reorganization.**

---

## Figure: `NC_vs_G50_GSEA_top_pathways.png`
### Caption
**Top significantly enriched GO Biological Process terms in the NC_vs_G50 preranked GSEA analysis.** Positive and negative NES pathways are shown separately to summarize the dominant biological axes in the strongest glycine-associated host contrast.

### One-line interpretation
**Compared with NC, the G50-associated host transcriptome showed the largest pathway-level shift, prominently involving extracellular structural remodeling together with vesicle-trafficking, proteostasis, and respiratory pathway changes.**

---

## Figure: `G50_vs_PC_GSEA_top_pathways.png`
### Caption
**Top significantly enriched GO Biological Process terms in the direct G50_vs_PC preranked GSEA analysis.** This comparison tests whether the high-dose glycine response is fully equivalent to the PC reference state at the pathway level.

### One-line interpretation
**Although G50 partially approached the PC reference pattern, substantial pathway differences remained, especially in catabolic, respiratory, and organelle-trafficking biology.**

---

## Figure: `NC_vs_G50_GSEA_enrichment_top_positive.png`
### Caption
**Representative enrichment plot for the top positive pathway in the NC_vs_G50 preranked GSEA analysis.** The running enrichment score indicates coordinated concentration of pathway genes toward one end of the ranked transcriptome.

### One-line interpretation
**This plot illustrates a coordinated shift in structural remodeling-related genes rather than a change driven by only a few extreme genes.**

---

## Figure: `NC_vs_G50_GSEA_enrichment_top_negative.png`
### Caption
**Representative enrichment plot for the top negative pathway in the NC_vs_G50 preranked GSEA analysis.**

### One-line interpretation
**This pattern supports coordinated redistribution of trafficking- and organelle-related programs in the G50-associated host response.**

---

## Figure: `G50_vs_PC_GSEA_enrichment_top_negative.png`
### Caption
**Representative enrichment plot for the top negative pathway in the direct G50_vs_PC preranked GSEA analysis.**

### One-line interpretation
**The direct G50_vs_PC comparison shows that the high-dose glycine state remains transcriptionally distinct from the PC reference state, particularly in metabolic and respiratory pathway structure.**

---

## Suggested spoken narrative for 1 slide
> We reran preranked GSEA using the RNA-seq full result tables and a Gallus gallus GO Biological Process gene-set collection. The strongest pathway remodeling was observed in the G50 contrast relative to NC. Across contrasts, recurrent themes included extracellular matrix and adhesion remodeling on one side and Golgi-vesicle, trafficking, proteostasis, and respiratory pathway reorganization on the other. Importantly, the direct G50_vs_PC comparison still retained a large number of significant pathways, indicating that G50 is only partially PC-like rather than fully equivalent to the PC reference state.

---

## Very short slide text version
- **G50 showed the strongest host pathway remodeling vs. NC.**
- **Recurring pathway themes: ECM/adhesion remodeling and Golgi-vesicle/metabolic reorganization.**
- **Direct G50_vs_PC GSEA indicates partial convergence, not full equivalence.**

---

## Recommended final wording to avoid overclaiming
Use:
- **partial convergence toward the PC reference state**
- **distinct host pathway remodeling in G50**
- **not fully equivalent to PC at the transcriptomic pathway level**

Avoid:
- **complete normalization**
- **fully mimicked PC**
- **host and PC became the same state**

---

## Sources
- `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/docs/glycine-gsea-rerun.md`
- `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/tables/rnaseq/fgsea_go_bp_summary.csv`
- `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/figures/rnaseq/GSEA_summary_dotplot.png`
- `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/figures/rnaseq/NC_vs_G50_GSEA_top_pathways.png`
- `file:///home/yzyzero/feynman/outputs/glycine-labmeeting-package-assets/figures/rnaseq/G50_vs_PC_GSEA_top_pathways.png`
