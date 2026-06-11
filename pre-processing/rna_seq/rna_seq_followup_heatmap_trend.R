options(stringsAsFactors = FALSE)

required_pkgs <- c("ggplot2", "dplyr", "readr", "tibble", "tidyr", "pheatmap")
install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, install_if_missing))

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(pheatmap)

base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(base_dir, "visualization_results")
font_family <- "sans"
palette_gly <- c(NC = "#75ACE4", PC = "#87BB8A", G25 = "#F9B066", G50 = "#E48282")

logcpm_path <- file.path(out_dir, "logTMM_CPM.csv")
sample_info_path <- file.path(out_dir, "sample_info.csv")
g50_path <- file.path(out_dir, "G50_vs_NC_edgeR_results.csv")
ann_path <- file.path(out_dir, "gene_annotation_from_gtf.csv")

stopifnot(file.exists(logcpm_path), file.exists(sample_info_path), file.exists(g50_path))

logcpm <- read.csv(logcpm_path, check.names = FALSE) %>% column_to_rownames("Geneid") %>% as.matrix()
sample_info <- read.csv(sample_info_path)
g50 <- read.csv(g50_path, check.names = FALSE)
ann <- if (file.exists(ann_path)) read.csv(ann_path, check.names = FALSE) else NULL

if (!is.null(ann) && !"Gene_name" %in% colnames(g50)) {
  g50 <- g50 %>% left_join(ann, by = "Geneid")
}
if (!"Gene_name" %in% colnames(g50)) g50$Gene_name <- NA

theme_dw <- function() {
  theme_minimal(base_family = font_family, base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.35),
      axis.title = element_text(color = "#222222"),
      axis.text = element_text(color = "#222222"),
      plot.title = element_text(face = "bold", size = 14, color = "#111111"),
      plot.subtitle = element_text(size = 10, color = "#555555"),
      legend.title = element_blank(),
      legend.position = "right"
    )
}

sample_info$Treat <- factor(sample_info$Treat, levels = c("NC", "G25", "G50", "PC"))
sample_info <- sample_info %>% arrange(Treat, sample)
logcpm <- logcpm[, sample_info$sample, drop = FALSE]

anno_col <- data.frame(Treat = sample_info$Treat)
rownames(anno_col) <- sample_info$sample
anno_colors <- list(Treat = palette_gly)

# Top DEG heatmap
heat_df <- g50 %>%
  mutate(label = ifelse(!is.na(Gene_name) & Gene_name != "", Gene_name, Geneid)) %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

if (nrow(heat_df) >= 10) {
  top_n <- min(60, nrow(heat_df))
  top_genes <- unique(heat_df$Geneid[1:top_n])
  heat_mat <- logcpm[top_genes, , drop = FALSE]
  row_labels <- heat_df$label[match(rownames(heat_mat), heat_df$Geneid)]
  rownames(heat_mat) <- make.unique(row_labels)
  pheatmap(
    heat_mat,
    filename = file.path(out_dir, "TopDEG_heatmap_G50_vs_NC.png"),
    width = 8.5,
    height = 12,
    scale = "row",
    clustering_method = "complete",
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    color = colorRampPalette(c("#4C78A8", "#F7F7F7", "#D65F5F"))(100),
    border_color = NA,
    main = "Top DEGs heatmap: G50 vs NC",
    fontsize = 10,
    fontsize_row = 7,
    fontsize_col = 8
  )
}

# Dose trend figure
trend_genes <- g50 %>%
  mutate(label = ifelse(!is.na(Gene_name) & Gene_name != "", Gene_name, Geneid)) %>%
  filter(FDR < 0.05, abs(logFC) >= 1) %>%
  arrange(FDR, desc(abs(logFC))) %>%
  distinct(Geneid, .keep_all = TRUE) %>%
  slice_head(n = 8) %>%
  select(Geneid, label)

if (nrow(trend_genes) > 0) {
  trend_mat <- logcpm[trend_genes$Geneid, , drop = FALSE]
  trend_df <- as.data.frame(t(trend_mat)) %>%
    rownames_to_column("sample") %>%
    left_join(sample_info, by = "sample") %>%
    pivot_longer(cols = all_of(trend_genes$Geneid), names_to = "Geneid", values_to = "logCPM") %>%
    left_join(trend_genes, by = "Geneid")

  trend_summary <- trend_df %>%
    group_by(Treat, label) %>%
    summarise(mean = mean(logCPM, na.rm = TRUE), se = sd(logCPM, na.rm = TRUE)/sqrt(sum(!is.na(logCPM))), .groups = "drop")

  p_trend <- ggplot(trend_summary, aes(Treat, mean, group = 1, color = Treat)) +
    geom_line(linewidth = 0.7, color = "#777777") +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12, linewidth = 0.4) +
    facet_wrap(~ label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = palette_gly) +
    labs(title = "Dose trend figure", subtitle = "Group mean logCPM ± SE for top G50-responsive genes", x = NULL, y = "mean logCPM") +
    theme_dw() +
    theme(strip.text = element_text(face = "bold", size = 9), legend.position = "none")

  ggsave(file.path(out_dir, "Dose_trend_top_genes.png"), p_trend, width = 12, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "Dose_trend_top_genes.pdf"), p_trend, width = 12, height = 7)
}

message("Follow-up plots completed: ", out_dir)
