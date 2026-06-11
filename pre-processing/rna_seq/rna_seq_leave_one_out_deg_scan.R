options(stringsAsFactors = FALSE)

required_pkgs <- c('BiocManager','edgeR','dplyr','tibble','readr')
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
library(dplyr)
library(tibble)
library(readr)

base_dir <- normalizePath('.', winslash='/', mustWork=TRUE)
counts_dir <- file.path(base_dir, '7.counts')
out_dir <- file.path(base_dir, 'leave_one_out_scan')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

extract_sample_name <- function(path) sub('_Sorted_count\\.txt$', '', basename(path))
extract_treat <- function(sample_name) {
  parts <- strsplit(sample_name, '-')[[1]]
  if (length(parts) >= 2 && tail(parts, 1) == 'R') tail(parts, 2)[1] else tail(parts, 1)
}

read_featurecounts <- function(path) {
  df <- read.delim(path, header=TRUE, sep='\t', check.names=FALSE, comment.char='#')
  tibble(Geneid = df[[1]], Count = df[[ncol(df)]])
}

count_files <- list.files(counts_dir, pattern = '_Sorted_count\\.txt$', full.names = TRUE)
count_files <- count_files[!grepl('\\.summary$', count_files)]
stopifnot(length(count_files) > 0)

sample_info <- tibble(sample = vapply(count_files, extract_sample_name, character(1)), file = count_files) %>%
  mutate(Treat = factor(vapply(sample, extract_treat, character(1)), levels = c('NC','G25','G50','PC')))

count_list <- lapply(sample_info$file, read_featurecounts)
merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count)
for (i in 2:length(count_list)) merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by='Geneid')
expr <- merged %>% column_to_rownames('Geneid') %>% as.matrix()
storage.mode(expr) <- 'integer'
expr <- expr[, sample_info$sample, drop=FALSE]

run_scan <- function(drop_sample = NULL) {
  keep_samples <- if (is.null(drop_sample)) rep(TRUE, ncol(expr)) else colnames(expr) != drop_sample
  expr2 <- expr[, keep_samples, drop=FALSE]
  info2 <- sample_info[match(colnames(expr2), sample_info$sample), , drop=FALSE]
  info2$Treat <- relevel(info2$Treat, ref='NC')

  y <- DGEList(counts=expr2, group=info2$Treat)
  keep <- rowSums(cpm(y) > 1) >= 3
  if (sum(keep) < 100) keep <- rowSums(y$counts >= 10) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method='TMM')

  design <- model.matrix(~ Treat, data=info2)
  colnames(design) <- make.names(colnames(design))

  y_lrt <- estimateGLMCommonDisp(y, design)
  y_lrt <- estimateGLMTrendedDisp(y_lrt, design)
  y_lrt <- estimateGLMTagwiseDisp(y_lrt, design)
  fit_lrt <- glmFit(y_lrt, design)
  lrt <- glmLRT(fit_lrt, coef='TreatG50')
  tab_lrt <- topTags(lrt, n=Inf, sort.by='PValue')$table

  y_qlf <- estimateDisp(y, design)
  fit_qlf <- glmQLFit(y_qlf, design, robust=TRUE)
  qlf <- glmQLFTest(fit_qlf, coef='TreatG50')
  tab_qlf <- topTags(qlf, n=Inf, sort.by='PValue')$table

  tibble(
    dropped = ifelse(is.null(drop_sample), 'NONE', drop_sample),
    lrt_up = sum(tab_lrt$FDR < 0.05 & tab_lrt$logFC >= 1, na.rm=TRUE),
    lrt_down = sum(tab_lrt$FDR < 0.05 & tab_lrt$logFC <= -1, na.rm=TRUE),
    lrt_total = lrt_up + lrt_down,
    qlf_up = sum(tab_qlf$FDR < 0.05 & tab_qlf$logFC >= 1, na.rm=TRUE),
    qlf_down = sum(tab_qlf$FDR < 0.05 & tab_qlf$logFC <= -1, na.rm=TRUE),
    qlf_total = qlf_up + qlf_down,
    n_samples = ncol(expr2)
  )
}

all_runs <- list()
all_runs[[1]] <- run_scan(NULL)
message('Completed: NONE')
for (i in seq_along(colnames(expr))) {
  s <- colnames(expr)[i]
  message('Running leave-one-out for: ', s)
  all_runs[[i + 1]] <- run_scan(s)
  write.csv(bind_rows(all_runs), file.path(out_dir, 'leave_one_out_G50_vs_NC_summary_partial.csv'), row.names=FALSE, quote=FALSE)
}
results <- bind_rows(all_runs) %>% arrange(desc(lrt_total), desc(qlf_total))
write.csv(results, file.path(out_dir, 'leave_one_out_G50_vs_NC_summary.csv'), row.names=FALSE, quote=FALSE)
print(results)
