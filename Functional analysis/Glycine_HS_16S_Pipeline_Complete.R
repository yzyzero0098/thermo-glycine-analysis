###############################################################################
# Glycine × Heat Stress 16S rRNA Microbiome Analysis Pipeline (Complete)
# - Alpha/Beta diversity, Taxonomy, DA, Centroid distance, Mimicking analysis
# - Dose-response trend, Functional prediction preparation
###############################################################################

# =============================================================================
# 0) Configuration
# =============================================================================
setwd("C:/Users/Desktop/Desktop/지영/glycine/qiime2/function")

feature_fp  <- "./feature-table.tsv"
taxonomy_fp <- "./taxonomy.tsv"
meta_fp     <- "./sample-metadata.txt"
tree_fp     <- "./tree.nwk"

out_dir <- paste0("results_16S_", format(Sys.Date(), "%Y%m%d"))
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

# =============================================================================
# 1) Package loading
# =============================================================================
cran_pkgs <- c("dplyr","tidyr","ggplot2","stringr","tibble","vegan",
               "ggpubr","pheatmap","patchwork","ape","conflicted")
bioc_pkgs <- c("phyloseq","ANCOMBC","ALDEx2","microbiomeMarker","Maaslin2")

# install.packages(cran_pkgs)
# BiocManager::install(bioc_pkgs)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(stringr)
  library(tibble); library(vegan); library(ggpubr); library(pheatmap)
  library(patchwork); library(ape); library(phyloseq)
})

if (requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
}

# =============================================================================
# 2) Plot theme & colors
# =============================================================================
base_family <- ifelse(Sys.info()[["sysname"]] == "Windows", "Arial", "Helvetica")

Treat_colors <- c(
  "NC"  = "#1976D2",
  "PC"  = "#388E3C",
  "G25" = "#F57C00",
  "G50" = "#D32F2F"
)

theme_set(
  theme_bw(base_family = base_family) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

# =============================================================================
# 3) Load data
# =============================================================================
meta <- read.delim(meta_fp, sep = "\t", header = TRUE, check.names = FALSE,
                   stringsAsFactors = FALSE)
meta$Treat <- factor(meta$Treat, levels = c("NC", "PC", "G25", "G50"))

feature <- read.delim(feature_fp, sep = "\t", skip = 1, row.names = 1,
                      check.names = FALSE, stringsAsFactors = FALSE)
if ("taxonomy" %in% colnames(feature)) feature <- feature[, colnames(feature) != "taxonomy"]
colnames(feature) <- gsub("-", ".", colnames(feature))
feature_mat <- as.matrix(feature)
storage.mode(feature_mat) <- "numeric"
feature_mat[is.na(feature_mat)] <- 0

tax <- read.delim(taxonomy_fp, sep = "\t", header = TRUE, check.names = FALSE)
colnames(tax)[1:2] <- c("FeatureID", "Taxon")
tax_split <- str_split_fixed(tax$Taxon, ";", 7) %>% apply(2, trimws) %>% as.matrix()
colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_split) <- tax$FeatureID
tax_split[tax_split == ""] <- NA
# Remove prefix
for (i in 1:ncol(tax_split)) {
  tax_split[, i] <- gsub("^[a-z]__", "", tax_split[, i])
}
tax_split[tax_split == ""] <- NA

tree <- ape::read.tree(tree_fp)

# Match samples & taxa
common_samples <- intersect(colnames(feature_mat), meta$SampleID)
common_taxa    <- intersect(rownames(feature_mat), rownames(tax_split))
feature_mat <- feature_mat[common_taxa, common_samples]
tax_mat     <- tax_split[rownames(feature_mat), ]
meta_use    <- meta %>% filter(SampleID %in% common_samples) %>%
  column_to_rownames("SampleID")
meta_use <- meta_use[common_samples, , drop = FALSE]

# Remove zero taxa
keep <- rowSums(feature_mat) > 0
feature_mat <- feature_mat[keep, ]
tax_mat     <- tax_mat[rownames(feature_mat), ]

# Build phyloseq
ps <- phyloseq(
  otu_table(feature_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_use),
  phy_tree(tree)
)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

cat(sprintf("[INFO] phyloseq: %d ASVs × %d samples\n", ntaxa(ps), nsamples(ps)))

# =============================================================================
# 4) Alpha Diversity
# =============================================================================
alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha_df$SampleID <- rownames(alpha_df)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")

# Faith's PD
library(picante)
otu_t <- t(as(otu_table(ps), "matrix"))
alpha_df$Faith_PD <- pd(otu_t, phy_tree(ps), include.root = TRUE)$PD

alpha_long <- pivot_longer(alpha_df, cols = c("Observed","Shannon","Simpson","Faith_PD"),
                           names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long, aes(x = Treat, y = Value, fill = Treat)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = Treat_colors) +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  labs(title = "Alpha Diversity", x = "Treatment", y = "Diversity Index")

ggsave(file.path(fig_dir, "Fig1_AlphaDiversity.png"), p_alpha, width = 12, height = 6, dpi = 300)
write.csv(alpha_df, file.path(tab_dir, "alpha_diversity.csv"), row.names = FALSE)

# =============================================================================
# 5) Beta Diversity — PCoA with filled ellipses
# =============================================================================
plot_pcoa <- function(ps, dist_method, title) {
  if (dist_method == "bray") {
    dist <- phyloseq::distance(ps, method = "bray")
  } else if (dist_method == "unifrac") {
    dist <- phyloseq::distance(ps, method = "unifrac", weighted = FALSE)
  } else if (dist_method == "wunifrac") {
    dist <- phyloseq::distance(ps, method = "unifrac", weighted = TRUE)
  }
  
  ord <- ordinate(ps, method = "PCoA", distance = dist)
  ord_df <- plot_ordination(ps, ord, justDF = TRUE)
  ord_df$Treat <- factor(sample_data(ps)$Treat, levels = c("NC","PC","G25","G50"))
  eig <- ord$values$Relative_eig
  
  perm <- adonis2(dist ~ Treat, data = data.frame(sample_data(ps)))
  pval <- perm$`Pr(>F)`[1]
  fstat <- perm$F[1]
  
  p_str <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.4f", pval))
  
  p <- ggplot(ord_df, aes(Axis.1, Axis.2, color = Treat, fill = Treat)) +
    stat_ellipse(geom = "polygon", type = "norm", alpha = 0.15,
                 color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", type = "norm", linewidth = 0.8,
                 linetype = "dashed", show.legend = FALSE) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = Treat_colors) +
    scale_fill_manual(values = Treat_colors) +
    labs(
      title = title,
      subtitle = sprintf("PERMANOVA: F = %.2f, %s", fstat, p_str),
      x = sprintf("PCoA1 (%.1f%%)", eig[1] * 100),
      y = sprintf("PCoA2 (%.1f%%)", eig[2] * 100)
    ) +
    theme(panel.grid = element_blank())
  
  return(p)
}

p_bray <- plot_pcoa(ps_rel, "bray", "Bray-Curtis")
p_uu   <- plot_pcoa(ps_rel, "unifrac", "Unweighted UniFrac")
p_wu   <- plot_pcoa(ps_rel, "wunifrac", "Weighted UniFrac")

p_beta <- p_bray + p_uu + p_wu + plot_layout(ncol = 3, guides = "collect")

ggsave(file.path(fig_dir, "Fig2_BetaDiversity_PCoA.png"), p_beta,
       width = 18, height = 6, dpi = 300)

# =============================================================================
# 6) Taxonomy Composition — Phylum & Genus
# =============================================================================
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
df_genus <- psmelt(ps_genus)
df_genus$Genus <- ifelse(is.na(df_genus$Genus) | df_genus$Genus == "",
                         "Unclassified", df_genus$Genus)

top10 <- df_genus %>% group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>% slice_head(n = 15) %>% pull(Genus)

df_genus$Genus_plot <- ifelse(df_genus$Genus %in% top10, df_genus$Genus, "Others")

p_tax <- ggplot(df_genus, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Treat, scales = "free_x", space = "free_x") +
  labs(title = "Genus-level Composition", x = "Sample", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

ggsave(file.path(fig_dir, "Fig3_Taxonomy_Genus.png"), p_tax, width = 14, height = 6, dpi = 300)

# =============================================================================
# 7) Differential Abundance — ANCOM-BC2
# =============================================================================
if (requireNamespace("ANCOMBC", quietly = TRUE)) {
  library(ANCOMBC)
  
  ancom_out <- ancombc2(
    data = ps,
    fix_formula = "Treat",
    rand_formula = NULL,
    group = "Treat",
    p_adj_method = "BH",
    prv_cut = 0.10,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE,
    global = TRUE,
    pairwise = TRUE,
    dunnet = TRUE,
    trend = FALSE
  )
  
  write.csv(as.data.frame(ancom_out$res), file.path(tab_dir, "DA_ANCOMBC2_primary.csv"))
  if (!is.null(ancom_out$res_global))
    write.csv(as.data.frame(ancom_out$res_global), file.path(tab_dir, "DA_ANCOMBC2_global.csv"))
  if (!is.null(ancom_out$res_pair))
    write.csv(as.data.frame(ancom_out$res_pair), file.path(tab_dir, "DA_ANCOMBC2_pairwise.csv"))
}

# =============================================================================
# 8) Differential Abundance — ALDEx2
# =============================================================================
if (requireNamespace("ALDEx2", quietly = TRUE)) {
  library(ALDEx2)
  aldex_clr <- aldex.clr(feature_mat, meta_use$Treat, mc.samples = 128, denom = "all")
  aldex_kw  <- aldex.kw(aldex_clr)
  aldex_eff <- aldex.effect(aldex_clr)
  
  aldex_tbl <- cbind(FeatureID = rownames(aldex_kw), aldex_kw,
                     aldex_eff[rownames(aldex_kw), ])
  write.csv(aldex_tbl, file.path(tab_dir, "DA_ALDEx2.csv"), row.names = FALSE)
}

# =============================================================================
# 9) Differential Abundance — MaAsLin2
# =============================================================================
if (requireNamespace("Maaslin2", quietly = TRUE)) {
  library(Maaslin2)
  Maaslin2(
    input_data = as.data.frame(t(feature_mat)),
    input_metadata = meta_use,
    output = file.path(out_dir, "maaslin2_out"),
    fixed_effects = "Treat",
    reference = "Treat,NC",
    normalization = "TSS",
    transform = "LOG",
    analysis_method = "LM"
  )
}

# =============================================================================
# 10) Centroid Distance Analysis
# =============================================================================
dist_bc <- phyloseq::distance(ps_rel, method = "bray")
dist_wu <- phyloseq::distance(ps_rel, method = "unifrac", weighted = TRUE)

calc_intergroup_dist <- function(dist_mat, meta, g1, g2) {
  dm <- as.matrix(dist_mat)
  s1 <- rownames(meta)[meta$Treat == g1]
  s2 <- rownames(meta)[meta$Treat == g2]
  dists <- dm[s1, s2, drop = FALSE]
  c(mean = mean(dists), sd = sd(dists))
}

centroid_res <- data.frame()
for (dm_name in c("BrayCurtis", "WeightedUniFrac")) {
  dm <- if (dm_name == "BrayCurtis") dist_bc else dist_wu
  for (g in c("NC", "G25", "G50")) {
    d <- calc_intergroup_dist(dm, meta_use, g, "PC")
    centroid_res <- rbind(centroid_res,
                          data.frame(Distance = dm_name, Group = g,
                                     MeanDist = d["mean"], SD = d["sd"]))
  }
}
write.csv(centroid_res, file.path(tab_dir, "centroid_distance_to_PC.csv"), row.names = FALSE)

# Centroid distance barplot
p_centroid <- ggplot(centroid_res, aes(x = Group, y = MeanDist, fill = Group)) +
  geom_col(alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_errorbar(aes(ymin = MeanDist - SD, ymax = MeanDist + SD), width = 0.2) +
  facet_wrap(~ Distance, scales = "free_y") +
  scale_fill_manual(values = Treat_colors) +
  labs(title = "Mean Distance to PC (Energy Control)",
       y = "Mean Inter-group Distance") +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig5_Centroid_Distance.png"), p_centroid,
       width = 10, height = 5, dpi = 300)

# =============================================================================
# 11) Glycine Mimicking Taxa
# =============================================================================
genus_rel <- psmelt(tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE))
genus_rel$Genus[is.na(genus_rel$Genus) | genus_rel$Genus == ""] <- "Unclassified"

genus_wide <- genus_rel %>%
  group_by(Genus, Sample, Treat) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)

pseudo <- 1e-6
nc_s  <- rownames(meta_use)[meta_use$Treat == "NC"]
pc_s  <- rownames(meta_use)[meta_use$Treat == "PC"]
g25_s <- rownames(meta_use)[meta_use$Treat == "G25"]
g50_s <- rownames(meta_use)[meta_use$Treat == "G50"]

# Calculate LFC at genus level using aggregated relative abundance
genus_agg <- as.data.frame(tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE) %>%
                             otu_table() %>% as("matrix"))

genus_tax <- as.data.frame(tax_table(tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)))
rownames(genus_agg) <- ifelse(is.na(genus_tax$Genus) | genus_tax$Genus == "",
                              "Unclassified", genus_tax$Genus)

lfc_pc  <- log2((rowMeans(genus_agg[, pc_s])  + pseudo) / (rowMeans(genus_agg[, nc_s]) + pseudo))
lfc_g25 <- log2((rowMeans(genus_agg[, g25_s]) + pseudo) / (rowMeans(genus_agg[, nc_s]) + pseudo))
lfc_g50 <- log2((rowMeans(genus_agg[, g50_s]) + pseudo) / (rowMeans(genus_agg[, nc_s]) + pseudo))

mimick <- data.frame(
  Genus = names(lfc_pc),
  LFC_PC_vs_NC  = lfc_pc,
  LFC_G25_vs_NC = lfc_g25,
  LFC_G50_vs_NC = lfc_g50,
  Delta_G25 = abs(lfc_g25 - lfc_pc),
  Delta_G50 = abs(lfc_g50 - lfc_pc)
)
mimick$Min_delta <- pmin(mimick$Delta_G25, mimick$Delta_G50)
mimick <- mimick[order(mimick$Min_delta), ]

write.csv(mimick, file.path(tab_dir, "glycine_mimicking_taxa.csv"), row.names = FALSE)

# =============================================================================
# 12) Dose-Response Trend (Spearman on glycine axis)
# =============================================================================
dose_map <- c("NC" = 0, "G25" = 0.341, "G50" = 0.683)
dose_meta <- meta_use[meta_use$Treat %in% c("NC","G25","G50"), , drop = FALSE]
dose_meta$Dose <- dose_map[as.character(dose_meta$Treat)]

genus_dose <- genus_agg[, rownames(dose_meta)]
trend_res <- data.frame()

for (g in rownames(genus_dose)) {
  ct <- cor.test(dose_meta$Dose, as.numeric(genus_dose[g, ]), method = "spearman",
                 exact = FALSE)
  trend_res <- rbind(trend_res,
                     data.frame(Genus = g, Spearman_rho = ct$estimate,
                                Spearman_p = ct$p.value))
}
trend_res$BH_p <- p.adjust(trend_res$Spearman_p, method = "BH")
trend_res <- trend_res[order(trend_res$Spearman_p), ]
write.csv(trend_res, file.path(tab_dir, "dose_response_trend.csv"), row.names = FALSE)

cat("\n[DONE] Complete pipeline finished.\n")
