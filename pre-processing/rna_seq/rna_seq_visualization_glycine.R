options(stringsAsFactors = FALSE)

required_pkgs <- c(
  "BiocManager", "edgeR", "limma", "ggplot2", "gridExtra",
  "dplyr", "readr", "tibble", "ggrepel", "tidyr", "pheatmap", "stringr"
)

install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("edgeR", "limma")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

invisible(lapply(required_pkgs, install_if_missing))

library(edgeR)
library(limma)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(readr)
library(tibble)
library(ggrepel)
library(tidyr)
library(pheatmap)
library(stringr)

base_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
counts_dir <- file.path(base_dir, "7.counts")
out_dir <- file.path(base_dir, "visualization_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf_candidates <- c(
  file.path(base_dir, "..", "reference", "gallus", "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.gtf"),
  file.path(base_dir, "..", "reference", "chicken", "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.107.gtf")
)
gtf_path <- gtf_candidates[file.exists(gtf_candidates)][1]
if (is.na(gtf_path)) gtf_path <- NULL

count_files <- list.files(counts_dir, pattern = "_Sorted_count\\.txt$", full.names = TRUE)
count_files <- count_files[!grepl("\\.summary$", count_files)]
if (length(count_files) == 0) stop("No count files found in 7.counts")

extract_sample_name <- function(path) sub("_Sorted_count\\.txt$", "", basename(path))
extract_treat <- function(sample_name) {
  parts <- strsplit(sample_name, "-")[[1]]
  if (length(parts) >= 2 && tail(parts, 1) == 'R') {
    tail(parts, 2)[1]
  } else {
    tail(parts, 1)
  }
}
extract_block <- function(sample_name) {
  parts <- strsplit(sample_name, "-")[[1]]
  paste(head(parts, 2), collapse = "-")
}

sample_names <- vapply(count_files, extract_sample_name, character(1))
group_levels <- c("NC", "G25", "G50", "PC")

sample_info <- tibble(
  sample = sample_names,
  file = count_files
)

read_featurecounts <- function(path) {
  df <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "#")
  count_col <- ncol(df)
  tibble(Geneid = df[[1]], Count = df[[count_col]])
}

count_list <- lapply(sample_info$file, read_featurecounts)
merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count)
for (i in 2:length(count_list)) {
  merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by = "Geneid")
}

expr <- merged %>% column_to_rownames("Geneid") %>% as.matrix()
storage.mode(expr) <- "integer"

sample_info <- tibble(sample = colnames(expr)) %>%
  mutate(
    file = file.path(counts_dir, paste0(sample, "_Sorted_count.txt")),
    Block = vapply(sample, extract_block, character(1)),
    Treat = factor(vapply(sample, extract_treat, character(1)), levels = group_levels)
  ) %>%
  arrange(Treat, sample)

expr <- expr[, sample_info$sample, drop = FALSE]
sample_info$Treat <- relevel(sample_info$Treat, ref = "NC")
message('expr dim before DGEList: ', paste(dim(expr), collapse=' x '))
message('sample_info rows: ', nrow(sample_info))
message('first expr cols: ', paste(head(colnames(expr)), collapse=', '))
message('first sample_info samples: ', paste(head(sample_info$sample), collapse=', '))
if (ncol(expr) != nrow(sample_info)) stop('Mismatch before DGEList: ncol(expr)=', ncol(expr), ', nrow(sample_info)=', nrow(sample_info))
y <- DGEList(counts = expr, group = sample_info$Treat)
keep <- rowSums(cpm(y) > 1) >= 3
if (sum(keep) < 100) {
  keep <- rowSums(y$counts >= 10) >= 2
}
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
message('Retained genes after filtering: ', nrow(y))

message('Treat counts:')
print(table(sample_info$Treat, useNA = 'ifany'))
design <- model.matrix(~ Treat, data = sample_info)
colnames(design) <- make.names(colnames(design))
message('design dim: ', paste(dim(design), collapse=' x '))
message('ncol(y): ', ncol(y))
if (nrow(design) != ncol(y)) stop('Design/count mismatch before dispersion: nrow(design)=', nrow(design), ', ncol(y)=', ncol(y))
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)

logcpm <- cpm(y, log = TRUE, prior.count = 2)
raw_cpm <- cpm(y, log = FALSE)

parse_gtf_annotation <- function(gtf_file) {
  lines <- readLines(gtf_file, warn = FALSE)
  lines <- lines[!startsWith(lines, "#")]
  gene_lines <- lines[grepl("\tgene\t", lines, fixed = FALSE)]
  if (length(gene_lines) == 0) return(NULL)
  split_tab <- strsplit(gene_lines, "\t")
  attrs <- vapply(split_tab, function(x) x[9], character(1))
  gene_id <- str_match(attrs, 'gene_id "([^"]+)"')[,2]
  gene_name <- str_match(attrs, 'gene_name "([^"]+)"')[,2]
  gene_bio <- str_match(attrs, 'gene_biotype "([^"]+)"')[,2]
  tibble(Geneid = gene_id, Gene_name = gene_name, Biotype = gene_bio) %>% distinct(Geneid, .keep_all = TRUE)
}

annotation_tbl <- NULL
if (!is.null(gtf_path)) {
  message("Parsing annotation from: ", gtf_path)
  annotation_tbl <- parse_gtf_annotation(gtf_path)
  if (!is.null(annotation_tbl)) {
    write.csv(annotation_tbl, file.path(out_dir, "gene_annotation_from_gtf.csv"), row.names = FALSE, quote = FALSE)
  }
}

annotate_results <- function(df) {
  out <- df %>% rownames_to_column("Geneid")
  if (!is.null(annotation_tbl)) {
    out <- out %>% left_join(annotation_tbl, by = "Geneid") %>% relocate(Gene_name, Biotype, .after = Geneid)
  }
  out
}

contrast_defs <- list(G25_vs_NC = 2, G50_vs_NC = 3, PC_vs_NC = 4)
all_results <- list()
for (nm in names(contrast_defs)) {
  lrt <- glmLRT(fit, coef = contrast_defs[[nm]])
  tab <- topTags(lrt, n = Inf, sort.by = "PValue")$table
  tab2 <- annotate_results(tab)
  all_results[[nm]] <- tab2
  write.csv(tab2, file.path(out_dir, paste0(nm, "_edgeR_results.csv")), row.names = FALSE, quote = FALSE)
  sig <- tab2 %>% filter(FDR < 0.05, abs(logFC) >= 1)
  write.csv(sig, file.path(out_dir, paste0(nm, "_DEG_FDR0.05_logFC1.csv")), row.names = FALSE, quote = FALSE)
}

palette_gly <- c(NC = "#75ACE4", PC = "#87BB8A", G25 = "#F9B066", G50 = "#E48282")
font_family <- "sans"

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

write.csv(raw_cpm %>% as.data.frame() %>% rownames_to_column("Geneid"), file.path(out_dir, "TMM_CPM.csv"), row.names = FALSE, quote = FALSE)
write.csv(logcpm %>% as.data.frame() %>% rownames_to_column("Geneid"), file.path(out_dir, "logTMM_CPM.csv"), row.names = FALSE, quote = FALSE)
write.csv(sample_info, file.path(out_dir, "sample_info.csv"), row.names = FALSE, quote = FALSE)

# MDS
mds <- plotMDS(y, top = 1000, plot = FALSE)
mds_df <- tibble(sample = colnames(y), Dim1 = mds$x, Dim2 = mds$y) %>% left_join(sample_info, by = "sample")
p_mds <- ggplot(mds_df, aes(Dim1, Dim2, color = Treat)) +
  geom_point(size = 3.8, alpha = 0.95) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3.1, family = font_family, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = palette_gly) +
  labs(title = "RNA-seq MDS plot", subtitle = "edgeR log-fold-change distance on TMM-normalized counts", x = "Dimension 1", y = "Dimension 2") +
  theme_dw()
ggsave(file.path(out_dir, "RNAseq_MDS_plot.png"), p_mds, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(out_dir, "RNAseq_MDS_plot.pdf"), p_mds, width = 7.2, height = 5.6)

# PCA
pca <- prcomp(t(logcpm), scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>% rownames_to_column("sample") %>% left_join(sample_info, by = "sample")
exp_var <- summary(pca)$importance[2, 1:2] * 100
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Treat)) +
  geom_point(size = 3.8, alpha = 0.95) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3.1, family = font_family, max.overlaps = 50, show.legend = FALSE) +
  scale_color_manual(values = palette_gly) +
  labs(title = "RNA-seq PCA plot", subtitle = "Principal components from logCPM matrix", x = paste0("PC1 (", round(exp_var[1], 1), "%)"), y = paste0("PC2 (", round(exp_var[2], 1), "%)")) +
  theme_dw()
ggsave(file.path(out_dir, "RNAseq_PCA_plot.png"), p_pca, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(out_dir, "RNAseq_PCA_plot.pdf"), p_pca, width = 7.2, height = 5.6)

make_volcano <- function(df, title_text, out_prefix) {
  label_col <- if ("Gene_name" %in% colnames(df)) "Gene_name" else "Geneid"
  df2 <- df %>% mutate(sig = case_when(FDR < 0.05 & logFC >= 1 ~ "Up", FDR < 0.05 & logFC <= -1 ~ "Down", TRUE ~ "NS"), neglog10FDR = -log10(pmax(FDR, 1e-300)))
  top_labels <- df2 %>% filter(sig != "NS") %>% arrange(FDR) %>% slice_head(n = 15)
  p <- ggplot(df2, aes(logFC, neglog10FDR)) +
    geom_point(aes(color = sig), size = 1.5, alpha = 0.85) +
    scale_color_manual(values = c(Down = "#75ACE4", NS = "#BDBDBD", Up = "#E48282")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
    ggrepel::geom_text_repel(data = top_labels, aes(label = .data[[label_col]]), family = font_family, size = 3, box.padding = 0.3, max.overlaps = 50, show.legend = FALSE) +
    labs(title = title_text, subtitle = "Thresholds: FDR < 0.05 and |log2FC| >= 1", x = expression(log[2]*" fold change"), y = expression(-log[10]*" FDR")) +
    theme_dw()
  ggsave(file.path(out_dir, paste0(out_prefix, "_volcano.png")), p, width = 6.3, height = 5.4, dpi = 300)
  ggsave(file.path(out_dir, paste0(out_prefix, "_volcano.pdf")), p, width = 6.3, height = 5.4)
}
make_volcano(all_results$G25_vs_NC, "G25 vs NC", "G25_vs_NC")
make_volcano(all_results$G50_vs_NC, "G50 vs NC", "G50_vs_NC")
make_volcano(all_results$PC_vs_NC, "PC vs NC", "PC_vs_NC")

# DEG summary
summary_tbl <- bind_rows(lapply(names(all_results), function(nm) {
  df <- all_results[[nm]]
  tibble(Contrast = nm, Up = sum(df$FDR < 0.05 & df$logFC >= 1, na.rm = TRUE), Down = sum(df$FDR < 0.05 & df$logFC <= -1, na.rm = TRUE))
}))
write.csv(summary_tbl, file.path(out_dir, "DEG_summary.csv"), row.names = FALSE, quote = FALSE)
summary_long <- summary_tbl %>% pivot_longer(cols = c("Up", "Down"), names_to = "Direction", values_to = "Count")
p_bar <- ggplot(summary_long, aes(Contrast, Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = c(Up = "#E48282", Down = "#75ACE4")) +
  labs(title = "Differentially expressed genes", subtitle = "FDR < 0.05 and |log2FC| >= 1", x = NULL, y = "Gene count") +
  theme_dw()
ggsave(file.path(out_dir, "DEG_summary_barplot.png"), p_bar, width = 6.8, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "DEG_summary_barplot.pdf"), p_bar, width = 6.8, height = 4.8)

# Correlation heatmap
anno_col <- data.frame(Treat = sample_info$Treat)
rownames(anno_col) <- sample_info$sample
anno_colors <- list(Treat = palette_gly)
sample_cor <- cor(logcpm, method = "pearson")
pheatmap(sample_cor, filename = file.path(out_dir, "RNAseq_sample_correlation_heatmap.png"), width = 8, height = 7, annotation_col = anno_col, annotation_row = anno_col, annotation_colors = anno_colors, color = colorRampPalette(c("#1F4E79", "#F7F7F7", "#E48282"))(100), border_color = NA, main = "Sample correlation heatmap", fontsize = 10, fontsize_row = 8, fontsize_col = 8)

# Top DEG heatmap using G50_vs_NC
heat_df <- all_results$G50_vs_NC %>% mutate(label = ifelse(!is.na(Gene_name) & Gene_name != "", Gene_name, Geneid)) %>% filter(FDR < 0.05) %>% arrange(FDR)
if (nrow(heat_df) >= 10) {
  top_n <- min(50, nrow(heat_df))
  top_genes <- unique(heat_df$Geneid[1:top_n])
  heat_mat <- logcpm[top_genes, , drop = FALSE]
  row_labels <- heat_df$label[match(rownames(heat_mat), heat_df$Geneid)]
  rownames(heat_mat) <- make.unique(row_labels)
  pheatmap(heat_mat, filename = file.path(out_dir, "TopDEG_heatmap_G50_vs_NC.png"), width = 8, height = 11, scale = "row", clustering_method = "complete", annotation_col = anno_col, annotation_colors = anno_colors, color = colorRampPalette(c("#4C78A8", "#F7F7F7", "#D65F5F"))(100), border_color = NA, main = "Top DEGs heatmap: G50 vs NC", fontsize = 10, fontsize_row = 7, fontsize_col = 8)
}

# Dose trend figure based on candidate genes from ordered contrasts
candidate_pool <- bind_rows(
  all_results$G25_vs_NC %>% mutate(Contrast = "G25_vs_NC"),
  all_results$G50_vs_NC %>% mutate(Contrast = "G50_vs_NC"),
  all_results$PC_vs_NC %>% mutate(Contrast = "PC_vs_NC")
) %>%
  mutate(label = ifelse(!is.na(Gene_name) & Gene_name != "", Gene_name, Geneid))

trend_genes <- candidate_pool %>%
  group_by(Geneid, label) %>%
  summarise(minFDR = min(FDR, na.rm = TRUE), maxAbsFC = max(abs(logFC), na.rm = TRUE), .groups = "drop") %>%
  arrange(minFDR, desc(maxAbsFC)) %>%
  slice_head(n = 8)

if (nrow(trend_genes) > 0) {
  trend_mat <- logcpm[trend_genes$Geneid, , drop = FALSE]
  trend_df <- as.data.frame(t(trend_mat)) %>% rownames_to_column("sample") %>% left_join(sample_info, by = "sample") %>%
    pivot_longer(cols = all_of(trend_genes$Geneid), names_to = "Geneid", values_to = "logCPM") %>%
    left_join(trend_genes, by = "Geneid")

  trend_summary <- trend_df %>% group_by(Treat, label) %>% summarise(mean = mean(logCPM, na.rm = TRUE), se = sd(logCPM, na.rm = TRUE)/sqrt(sum(!is.na(logCPM))), .groups = "drop")
  trend_summary$Treat <- factor(trend_summary$Treat, levels = c("NC", "G25", "G50", "PC"))

  p_trend <- ggplot(trend_summary, aes(Treat, mean, group = 1, color = Treat)) +
    geom_line(linewidth = 0.7, color = "#777777") +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12, linewidth = 0.4) +
    facet_wrap(~ label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = palette_gly) +
    labs(title = "Dose trend figure", subtitle = "Group mean logCPM ± SE for top candidate genes", x = NULL, y = "mean logCPM") +
    theme_dw() +
    theme(strip.text = element_text(face = "bold", size = 9), legend.position = "none")

  ggsave(file.path(out_dir, "Dose_trend_top_genes.png"), p_trend, width = 12, height = 7, dpi = 300)
  ggsave(file.path(out_dir, "Dose_trend_top_genes.pdf"), p_trend, width = 12, height = 7)
}

message("RNA-seq visualization pipeline completed. Outputs written to: ", out_dir)
