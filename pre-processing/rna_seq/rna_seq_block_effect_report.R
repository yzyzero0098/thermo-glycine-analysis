options(stringsAsFactors = FALSE)

required_pkgs <- c('BiocManager','edgeR','ggplot2','dplyr','tibble','tidyr')
install_if_missing <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg == 'edgeR') {
      if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')
      BiocManager::install(pkg, ask=FALSE, update=FALSE)
    } else {
      install.packages(pkg, repos='https://cloud.r-project.org')
    }
  }
}
invisible(lapply(required_pkgs, install_if_missing))

library(edgeR)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

base_dir <- normalizePath('.', winslash='/', mustWork=TRUE)
counts_dir <- file.path(base_dir, '7.counts')
out_dir <- file.path(base_dir, 'block_effect_report')
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

palette_comp <- c('Treat only' = '#B0B0B0', 'Block + Treat' = '#2D9CDB')
font_family <- 'sans'

extract_sample_name <- function(path) sub('_Sorted_count\\.txt$', '', basename(path))
extract_treat <- function(sample_name) { parts <- strsplit(sample_name, '-')[[1]]; if (length(parts) >= 2 && tail(parts,1) == 'R') tail(parts,2)[1] else tail(parts,1) }
extract_block <- function(sample_name) paste(head(strsplit(sample_name, '-')[[1]], 2), collapse='-')
read_featurecounts <- function(path) { df <- read.delim(path, header=TRUE, sep='\t', check.names=FALSE, comment.char='#'); tibble(Geneid = df[[1]], Count = df[[ncol(df)]]) }

theme_dw <- function() {
  theme_minimal(base_family = font_family, base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = '#EAEAEA', linewidth = 0.35),
      axis.title = element_text(color = '#111111'),
      axis.text = element_text(color = '#111111'),
      axis.line = element_line(color = '#444444', linewidth = 0.4),
      axis.ticks = element_line(color = '#444444', linewidth = 0.35),
      plot.title = element_text(face='bold', size=14, color='#111111'),
      plot.subtitle = element_text(size=10, color='#555555'),
      legend.title = element_blank(),
      legend.position = 'right'
    )
}

count_files <- list.files(counts_dir, pattern='_Sorted_count\\.txt$', full.names=TRUE)
count_files <- count_files[!grepl('\\.summary$', count_files)]
sample_info <- tibble(sample=vapply(count_files, extract_sample_name, character(1)), file=count_files) %>%
  mutate(Treat=factor(vapply(sample, extract_treat, character(1)), levels=c('NC','G25','G50','PC')),
         Block=vapply(sample, extract_block, character(1)))
count_list <- lapply(sample_info$file, read_featurecounts)
merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count)
for (i in 2:length(count_list)) merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by='Geneid')
expr <- merged %>% column_to_rownames('Geneid') %>% as.matrix(); storage.mode(expr) <- 'integer'; expr <- expr[, sample_info$sample, drop=FALSE]

run_model <- function(use_block=FALSE, method='LRT') {
  info <- sample_info
  info$Treat <- relevel(info$Treat, ref='NC')
  y <- DGEList(counts=expr, group=info$Treat)
  keep <- rowSums(cpm(y) > 1) >= 3
  if (sum(keep) < 100) keep <- rowSums(y$counts >= 10) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  design <- if (use_block) model.matrix(~ Block + Treat, data=info) else model.matrix(~ Treat, data=info)
  colnames(design) <- make.names(colnames(design))

  if (method == 'LRT') {
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    test_fun <- function(coef_name) glmLRT(fit, coef=coef_name)
  } else {
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    test_fun <- function(coef_name) glmQLFTest(fit, coef=coef_name)
  }

  res <- bind_rows(lapply(c('G25','G50','PC'), function(tr) {
    coef_name <- paste0('Treat', tr)
    tt <- topTags(test_fun(coef_name), n=Inf, sort.by='PValue')$table
    tibble(
      Contrast = paste0(tr, '_vs_NC'),
      Up = sum(tt$FDR < 0.05 & tt$logFC >= 1, na.rm=TRUE),
      Down = sum(tt$FDR < 0.05 & tt$logFC <= -1, na.rm=TRUE)
    ) %>% mutate(Total = Up + Down)
  }))
  res
}

lrt_treat <- run_model(use_block=FALSE, method='LRT') %>% mutate(Model='Treat only', Method='LRT')
lrt_block <- run_model(use_block=TRUE, method='LRT') %>% mutate(Model='Block + Treat', Method='LRT')
qlf_treat <- run_model(use_block=FALSE, method='QLF') %>% mutate(Model='Treat only', Method='QLF')
qlf_block <- run_model(use_block=TRUE, method='QLF') %>% mutate(Model='Block + Treat', Method='QLF')

summary_tbl <- bind_rows(lrt_treat, lrt_block, qlf_treat, qlf_block)
write.csv(summary_tbl, file.path(out_dir, 'DEG_counts_treat_vs_block.csv'), row.names=FALSE, quote=FALSE)

# report table focused on LRT first
report_tbl <- bind_rows(lrt_treat, lrt_block)
write.csv(report_tbl, file.path(out_dir, 'DEG_counts_LRT_report.csv'), row.names=FALSE, quote=FALSE)

plot_df <- report_tbl %>% select(Contrast, Model, Total)
p <- ggplot(plot_df, aes(Contrast, Total, fill=Model)) +
  geom_col(position=position_dodge(width=0.72), width=0.64) +
  scale_fill_manual(values=palette_comp) +
  labs(title='DEG counts before and after block adjustment', subtitle='Cutoff: FDR < 0.05 and |log2FC| >= 1 (LRT)', x=NULL, y='Number of DEGs') +
  theme_dw()
ggsave(file.path(out_dir, 'DEG_counts_before_after_block_LRT.png'), p, width=7.5, height=5.2, dpi=300)
ggsave(file.path(out_dir, 'DEG_counts_before_after_block_LRT.pdf'), p, width=7.5, height=5.2, dpi=300)

plot_df2 <- summary_tbl %>% mutate(Label = paste(Model, Method, sep=' / '))
p2 <- ggplot(plot_df2, aes(Contrast, Total, fill=Label)) +
  geom_col(position=position_dodge(width=0.8), width=0.68) +
  labs(title='DEG count comparison across model settings', subtitle='LRT and QLF, with and without block adjustment', x=NULL, y='Number of DEGs') +
  theme_dw()
ggsave(file.path(out_dir, 'DEG_counts_model_settings.png'), p2, width=9.2, height=5.5, dpi=300)

# report text
report_lines <- c(
  'Batch effect interpretation',
  '- In this dataset, Block corresponds to the repeated sample set identifier (e.g., P1-1, P1-2, P1-3).',
  '- If these identifiers reflect shared non-treatment variation (batch/run/paired set/cage-level background), then ~ Block + Treat adjusts for that source of variation before estimating treatment effects.',
  '',
  'Model formulas',
  '- Treat only: y ~ Treat',
  '- Block-adjusted: y ~ Block + Treat',
  '',
  'DEG cutoff',
  '- FDR < 0.05 and |log2FC| >= 1',
  '',
  'Interpretation note',
  '- Treat only maximizes sensitivity but may absorb block-level heterogeneity into treatment effects.',
  '- Block + Treat is more conservative and can be preferable when the repeated identifier captures real batch-like structure.'
)
writeLines(report_lines, file.path(out_dir, 'report_notes.txt'))

message('Block effect report completed: ', out_dir)
