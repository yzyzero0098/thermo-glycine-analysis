options(stringsAsFactors = FALSE)

required_pkgs <- c("BiocManager", "edgeR", "limma", "ggplot2", "dplyr", "readr", "tibble", "tidyr", "ggrepel", "pheatmap", "stringr")
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
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(stringr)

base_dir <- normalizePath('.', winslash='/', mustWork=TRUE)
counts_dir <- file.path(base_dir, '7.counts')
out_dir <- file.path(base_dir, 'visualization_results_upgrade')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

palette_gly <- c(NC = '#75ACE4', PC = '#87BB8A', G25 = '#F9B066', G50 = '#E48282')
font_family <- 'sans'

extract_sample_name <- function(path) sub('_Sorted_count\\.txt$', '', basename(path))
extract_treat <- function(sample_name) {
  parts <- strsplit(sample_name, '-')[[1]]
  if (length(parts) >= 2 && tail(parts, 1) == 'R') tail(parts, 2)[1] else tail(parts, 1)
}
extract_block <- function(sample_name) paste(head(strsplit(sample_name, '-')[[1]], 2), collapse='-')

read_featurecounts <- function(path) {
  df <- read.delim(path, header=TRUE, sep='\t', check.names=FALSE, comment.char='#')
  tibble(Geneid = df[[1]], Count = df[[ncol(df)]])
}

theme_dw <- function() {
  theme_minimal(base_family = font_family, base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = '#EAEAEA', linewidth = 0.35),
      axis.title = element_text(color = '#111111', face='plain'),
      axis.text = element_text(color = '#111111'),
      axis.line = element_line(color = '#444444', linewidth = 0.4),
      axis.ticks = element_line(color = '#444444', linewidth = 0.35),
      plot.title = element_text(face = 'bold', size = 14, color = '#111111'),
      plot.subtitle = element_text(size = 10, color = '#555555'),
      legend.title = element_blank(),
      legend.position = 'right'
    )
}

count_files <- list.files(counts_dir, pattern = '_Sorted_count\\.txt$', full.names = TRUE)
count_files <- count_files[!grepl('\\.summary$', count_files)]
stopifnot(length(count_files) > 0)

sample_info <- tibble(sample = vapply(count_files, extract_sample_name, character(1)), file = count_files) %>%
  mutate(Block = vapply(sample, extract_block, character(1)), Treat = factor(vapply(sample, extract_treat, character(1)), levels = c('NC','G25','G50','PC'))) %>%
  arrange(Treat, sample)

count_list <- lapply(sample_info$file, read_featurecounts)
merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count)
for (i in 2:length(count_list)) merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by='Geneid')
expr <- merged %>% column_to_rownames('Geneid') %>% as.matrix()
storage.mode(expr) <- 'integer'
expr <- expr[, sample_info$sample, drop=FALSE]

# basic filtering
sample_info$Treat <- relevel(sample_info$Treat, ref='NC')
y <- DGEList(counts=expr, group=sample_info$Treat)
keep <- rowSums(cpm(y) > 1) >= 3
if (sum(keep) < 100) keep <- rowSums(y$counts >= 10) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method='TMM')
logcpm <- cpm(y, log=TRUE, prior.count=2)

# QC metrics
qc_tbl <- tibble(
  sample = colnames(y),
  lib_size = y$samples$lib.size,
  norm_factor = y$samples$norm.factors,
  detected_genes = colSums(y$counts > 0)
) %>% left_join(sample_info, by='sample')
write.csv(qc_tbl, file.path(out_dir, 'sample_qc_metrics.csv'), row.names=FALSE, quote=FALSE)

# MDS/PCA helper
make_group_polygon_plot <- function(df, xvar, yvar, title_text, subtitle_text, out_prefix) {
  p <- ggplot(df, aes(.data[[xvar]], .data[[yvar]], color=Treat, fill=Treat)) +
    stat_ellipse(geom='polygon', alpha=0.14, linewidth=0.3, level=0.80, show.legend=FALSE) +
    geom_point(size=3.2, alpha=0.95) +
    ggrepel::geom_text_repel(aes(label=sample), size=3.1, family=font_family, max.overlaps=50, show.legend=FALSE) +
    scale_color_manual(values=palette_gly) +
    scale_fill_manual(values=palette_gly) +
    labs(title=title_text, subtitle=subtitle_text, x=xvar, y=yvar) +
    theme_dw()
  ggsave(file.path(out_dir, paste0(out_prefix, '.png')), p, width=7.4, height=5.8, dpi=300)
  ggsave(file.path(out_dir, paste0(out_prefix, '.pdf')), p, width=7.4, height=5.8, dpi=300)
}

mds <- plotMDS(y, top=1000, plot=FALSE)
mds_df <- tibble(sample=colnames(y), Dim1=mds$x, Dim2=mds$y) %>% left_join(sample_info, by='sample')
make_group_polygon_plot(mds_df, 'Dim1', 'Dim2', 'RNA-seq MDS plot', 'edgeR log-fold-change distance with group polygons', 'RNAseq_MDS_polygon')

pca <- prcomp(t(logcpm), scale.=TRUE)
pca_df <- as.data.frame(pca$x[,1:2]) %>% rownames_to_column('sample') %>% left_join(sample_info, by='sample')
colnames(pca_df)[2:3] <- c('PC1','PC2')
make_group_polygon_plot(pca_df, 'PC1', 'PC2', 'RNA-seq PCA plot', 'PC scores with group polygons', 'RNAseq_PCA_polygon')

# QC figures
p_lib <- ggplot(qc_tbl, aes(sample, lib_size/1e6, fill=Treat)) + geom_col(width=0.7) + scale_fill_manual(values=palette_gly) + labs(title='Library size', x=NULL, y='Millions of reads') + theme_dw() + theme(axis.text.x=element_text(angle=60, hjust=1), legend.position='none')
ggsave(file.path(out_dir, 'QC_library_size.png'), p_lib, width=9, height=4.8, dpi=300)

p_detect <- ggplot(qc_tbl, aes(sample, detected_genes, fill=Treat)) + geom_col(width=0.7) + scale_fill_manual(values=palette_gly) + labs(title='Detected genes per sample', x=NULL, y='Detected genes') + theme_dw() + theme(axis.text.x=element_text(angle=60, hjust=1), legend.position='none')
ggsave(file.path(out_dir, 'QC_detected_genes.png'), p_detect, width=9, height=4.8, dpi=300)

box_df <- as.data.frame(logcpm) %>% rownames_to_column('Geneid') %>% pivot_longer(-Geneid, names_to='sample', values_to='logCPM') %>% left_join(sample_info, by='sample')
p_box <- ggplot(box_df, aes(sample, logCPM, fill=Treat)) + geom_boxplot(outlier.size=0.25, linewidth=0.25) + scale_fill_manual(values=palette_gly) + labs(title='logCPM distribution', x=NULL, y='logCPM') + theme_dw() + theme(axis.text.x=element_text(angle=60, hjust=1), legend.position='none')
ggsave(file.path(out_dir, 'QC_logCPM_boxplot.png'), p_box, width=10, height=5.2, dpi=300)

# correlation heatmap
anno_col <- data.frame(Treat = sample_info$Treat)
rownames(anno_col) <- sample_info$sample
anno_colors <- list(Treat = palette_gly)
pheatmap(cor(logcpm), filename=file.path(out_dir, 'QC_sample_correlation_heatmap.png'), width=8, height=7, annotation_col=anno_col, annotation_row=anno_col, annotation_colors=anno_colors, color=colorRampPalette(c('#1F4E79','#F7F7F7','#E48282'))(100), border_color=NA, main='Sample correlation heatmap', fontsize=10, fontsize_row=8, fontsize_col=8)

# models
run_edger_lrt <- function(y_obj, design, coef_name_or_idx, label) {
  y2 <- estimateGLMCommonDisp(y_obj, design)
  y2 <- estimateGLMTrendedDisp(y2, design)
  y2 <- estimateGLMTagwiseDisp(y2, design)
  fit <- glmFit(y2, design)
  coef_idx <- if (is.character(coef_name_or_idx)) which(colnames(design) == coef_name_or_idx) else coef_name_or_idx
  lrt <- glmLRT(fit, coef=coef_idx)
  tab <- topTags(lrt, n=Inf, sort.by='PValue')$table %>% rownames_to_column('Geneid')
  write.csv(tab, file.path(out_dir, paste0(label, '_LRT.csv')), row.names=FALSE, quote=FALSE)
  c(up = sum(tab$FDR < 0.05 & tab$logFC >= 1, na.rm=TRUE), down = sum(tab$FDR < 0.05 & tab$logFC <= -1, na.rm=TRUE))
}

run_edger_qlf <- function(y_obj, design, coef_name_or_idx, label) {
  y2 <- estimateDisp(y_obj, design)
  fit <- glmQLFit(y2, design, robust=TRUE)
  coef_idx <- if (is.character(coef_name_or_idx)) which(colnames(design) == coef_name_or_idx) else coef_name_or_idx
  qlf <- glmQLFTest(fit, coef=coef_idx)
  tab <- topTags(qlf, n=Inf, sort.by='PValue')$table %>% rownames_to_column('Geneid')
  write.csv(tab, file.path(out_dir, paste0(label, '_QLF.csv')), row.names=FALSE, quote=FALSE)
  c(up = sum(tab$FDR < 0.05 & tab$logFC >= 1, na.rm=TRUE), down = sum(tab$FDR < 0.05 & tab$logFC <= -1, na.rm=TRUE))
}

# model 1: Treat only
sample_info$Treat <- relevel(sample_info$Treat, ref='NC')
design_treat <- model.matrix(~ Treat, data=sample_info)
colnames(design_treat) <- make.names(colnames(design_treat))

# model 2: Block + Treat (drop incomplete blocks impossible? keep additive block model)
design_block <- model.matrix(~ Block + Treat, data=sample_info)
colnames(design_block) <- make.names(colnames(design_block))

# optional outlier exclusion candidate
candidate_outlier <- 'P1-7-JM-G50-R'
keep_samples <- colnames(y) != candidate_outlier
y_no <- y[, keep_samples, keep.lib.sizes=FALSE]
sample_info_no <- sample_info[keep_samples, , drop=FALSE]
design_treat_no <- model.matrix(~ Treat, data=sample_info_no)
colnames(design_treat_no) <- make.names(colnames(design_treat_no))

summary_rows <- list(
  cbind(Model='Treat_only', Method='LRT', Contrast='G25_vs_NC', t(run_edger_lrt(y, design_treat, 'TreatG25', 'TreatOnly_G25_vs_NC'))),
  cbind(Model='Treat_only', Method='LRT', Contrast='G50_vs_NC', t(run_edger_lrt(y, design_treat, 'TreatG50', 'TreatOnly_G50_vs_NC'))),
  cbind(Model='Treat_only', Method='LRT', Contrast='PC_vs_NC', t(run_edger_lrt(y, design_treat, 'TreatPC', 'TreatOnly_PC_vs_NC'))),
  cbind(Model='Treat_only', Method='QLF', Contrast='G25_vs_NC', t(run_edger_qlf(y, design_treat, 'TreatG25', 'TreatOnly_G25_vs_NC'))),
  cbind(Model='Treat_only', Method='QLF', Contrast='G50_vs_NC', t(run_edger_qlf(y, design_treat, 'TreatG50', 'TreatOnly_G50_vs_NC'))),
  cbind(Model='Treat_only', Method='QLF', Contrast='PC_vs_NC', t(run_edger_qlf(y, design_treat, 'TreatPC', 'TreatOnly_PC_vs_NC'))),
  cbind(Model='Block_plus_Treat', Method='LRT', Contrast='G50_vs_NC', t(run_edger_lrt(y, design_block, 'TreatG50', 'BlockTreat_G50_vs_NC'))),
  cbind(Model='Block_plus_Treat', Method='QLF', Contrast='G50_vs_NC', t(run_edger_qlf(y, design_block, 'TreatG50', 'BlockTreat_G50_vs_NC'))),
  cbind(Model='Treat_only_minus_G50_7', Method='LRT', Contrast='G50_vs_NC', t(run_edger_lrt(y_no, design_treat_no, 'TreatG50', 'TreatOnly_minusG50_7_G50_vs_NC'))),
  cbind(Model='Treat_only_minus_G50_7', Method='QLF', Contrast='G50_vs_NC', t(run_edger_qlf(y_no, design_treat_no, 'TreatG50', 'TreatOnly_minusG50_7_G50_vs_NC')))
)
summary_tbl <- bind_rows(lapply(summary_rows, as.data.frame))
summary_tbl$up <- as.integer(summary_tbl$up)
summary_tbl$down <- as.integer(summary_tbl$down)
summary_tbl$total <- summary_tbl$up + summary_tbl$down
write.csv(summary_tbl, file.path(out_dir, 'DEG_model_comparison_summary.csv'), row.names=FALSE, quote=FALSE)

# integrated heatmap across all groups based on G50 vs NC from Treat-only LRT
res_g50 <- read.csv(file.path(out_dir, 'TreatOnly_G50_vs_NC_LRT.csv'))
sel <- res_g50 %>% filter(FDR < 0.05) %>% arrange(FDR) %>% slice_head(n=60)
if (nrow(sel) >= 10) {
  heat_mat <- logcpm[sel$Geneid, , drop=FALSE]
  rownames(heat_mat) <- make.unique(sel$Geneid)
  pheatmap(heat_mat, filename=file.path(out_dir, 'Integrated_heatmap_all_groups_from_G50_selection.png'), width=8.5, height=12, scale='row', annotation_col=anno_col, annotation_colors=anno_colors, color=colorRampPalette(c('#4C78A8','#F7F7F7','#D65F5F'))(100), border_color=NA, main='Integrated heatmap across NC/G25/G50/PC')
}

message('Upgrade QC/block/outlier workflow completed: ', out_dir)
