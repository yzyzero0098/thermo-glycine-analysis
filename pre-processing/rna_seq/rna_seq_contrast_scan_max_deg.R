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
out_dir <- file.path(base_dir, 'contrast_scan_max_deg')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

count_files <- list.files(counts_dir, pattern = '_Sorted_count\\.txt$', full.names = TRUE)
count_files <- count_files[!grepl('\\.summary$', count_files)]
stopifnot(length(count_files) > 0)

sample_info <- tibble(sample = vapply(count_files, extract_sample_name, character(1)), file = count_files) %>%
  mutate(Treat = factor(vapply(sample, extract_treat, character(1)), levels = c('NC','G25','G50','PC')),
         Block = vapply(sample, extract_block, character(1)))

count_list <- lapply(sample_info$file, read_featurecounts)
merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count)
for (i in 2:length(count_list)) merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by='Geneid')
expr <- merged %>% column_to_rownames('Geneid') %>% as.matrix()
storage.mode(expr) <- 'integer'
expr <- expr[, sample_info$sample, drop=FALSE]

run_one <- function(info, counts, contrast_name, method='LRT', use_block=FALSE) {
  info$Treat <- relevel(info$Treat, ref='NC')
  y <- DGEList(counts=counts, group=info$Treat)
  keep <- rowSums(cpm(y) > 1) >= 3
  if (sum(keep) < 100) keep <- rowSums(y$counts >= 10) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method='TMM')
  design <- if (use_block) model.matrix(~ Block + Treat, data=info) else model.matrix(~ Treat, data=info)
  colnames(design) <- make.names(colnames(design))
  coef_name <- paste0('Treat', contrast_name)
  if (!(coef_name %in% colnames(design))) return(NULL)

  if (method == 'LRT') {
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTrendedDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    test <- glmLRT(fit, coef=coef_name)
  } else {
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    test <- glmQLFTest(fit, coef=coef_name)
  }
  tab <- topTags(test, n=Inf, sort.by='PValue')$table
  tibble(
    up = sum(tab$FDR < 0.05 & tab$logFC >= 1, na.rm=TRUE),
    down = sum(tab$FDR < 0.05 & tab$logFC <= -1, na.rm=TRUE)
  ) %>% mutate(total = up + down)
}

scan_contrast <- function(contrast_name, method='LRT', use_block=FALSE, restrict_groups=NULL) {
  info0 <- sample_info
  counts0 <- expr
  if (!is.null(restrict_groups)) {
    keep_idx <- info0$Treat %in% restrict_groups
    info0 <- info0[keep_idx, , drop=FALSE]
    counts0 <- counts0[, info0$sample, drop=FALSE]
    info0$Treat <- droplevels(info0$Treat)
  }

  runs <- list()
  base <- run_one(info0, counts0, contrast_name, method=method, use_block=use_block)
  if (is.null(base)) return(NULL)
  runs[[1]] <- cbind(dropped='NONE', base)

  for (s in info0$sample) {
    info2 <- info0[info0$sample != s, , drop=FALSE]
    counts2 <- counts0[, info2$sample, drop=FALSE]
    out <- tryCatch(run_one(info2, counts2, contrast_name, method=method, use_block=use_block), error=function(e) NULL)
    if (!is.null(out)) runs[[length(runs)+1]] <- cbind(dropped=s, out)
  }
  bind_rows(runs) %>% arrange(desc(total), desc(up))
}

# 1) Main scans for G25 vs NC / PC vs NC / G50 vs PC
jobs <- list(
  list(name='G25_vs_NC_LRT', contrast='G25', method='LRT', use_block=FALSE, groups=c('NC','G25','G50','PC')),
  list(name='G25_vs_NC_QLF', contrast='G25', method='QLF', use_block=FALSE, groups=c('NC','G25','G50','PC')),
  list(name='G25_vs_NC_Block_LRT', contrast='G25', method='LRT', use_block=TRUE, groups=c('NC','G25','G50','PC')),
  list(name='G25_vs_NC_Block_QLF', contrast='G25', method='QLF', use_block=TRUE, groups=c('NC','G25','G50','PC')),
  list(name='PC_vs_NC_LRT', contrast='PC', method='LRT', use_block=FALSE, groups=c('NC','G25','G50','PC')),
  list(name='PC_vs_NC_QLF', contrast='PC', method='QLF', use_block=FALSE, groups=c('NC','G25','G50','PC')),
  list(name='G50_vs_PC_onlygroups_LRT', contrast='TreatG50', method='LRT', use_block=FALSE, groups=c('G50','PC'))
)

all_summaries <- list()
for (job in jobs) {
  message('Running ', job$name)
  if (job$name == 'G50_vs_PC_onlygroups_LRT') {
    info0 <- sample_info[sample_info$Treat %in% c('G50','PC'), , drop=FALSE]
    counts0 <- expr[, info0$sample, drop=FALSE]
    info0$Treat <- factor(info0$Treat, levels=c('PC','G50'))
    run_gp <- function(info, counts) {
      y <- DGEList(counts=counts, group=info$Treat)
      keep <- rowSums(cpm(y) > 1) >= 3
      if (sum(keep) < 100) keep <- rowSums(y$counts >= 10) >= 2
      y <- y[keep, , keep.lib.sizes=FALSE]
      y <- calcNormFactors(y)
      design <- model.matrix(~ Treat, data=info)
      colnames(design) <- make.names(colnames(design))
      y <- estimateGLMCommonDisp(y, design)
      y <- estimateGLMTrendedDisp(y, design)
      y <- estimateGLMTagwiseDisp(y, design)
      fit <- glmFit(y, design)
      lrt <- glmLRT(fit, coef='TreatG50')
      tab <- topTags(lrt, n=Inf, sort.by='PValue')$table
      tibble(up=sum(tab$FDR < 0.05 & tab$logFC >= 1, na.rm=TRUE), down=sum(tab$FDR < 0.05 & tab$logFC <= -1, na.rm=TRUE)) %>% mutate(total=up+down)
    }
    runs <- list(cbind(dropped='NONE', run_gp(info0, counts0)))
    for (s in info0$sample) {
      info2 <- info0[info0$sample != s, , drop=FALSE]
      counts2 <- counts0[, info2$sample, drop=FALSE]
      out <- tryCatch(run_gp(info2, counts2), error=function(e) NULL)
      if (!is.null(out)) runs[[length(runs)+1]] <- cbind(dropped=s, out)
    }
    res <- bind_rows(runs) %>% arrange(desc(total), desc(up))
  } else {
    res <- scan_contrast(job$contrast, method=job$method, use_block=job$use_block, restrict_groups=job$groups)
  }
  if (!is.null(res)) {
    write.csv(res, file.path(out_dir, paste0(job$name, '.csv')), row.names=FALSE, quote=FALSE)
    best <- res %>% slice(1) %>% mutate(scan=job$name)
    all_summaries[[length(all_summaries)+1]] <- best
  }
}

# 2) Pairwise-only subset for G25 vs NC (drop G50/PC entirely)
run_pair_only <- function(method='LRT', use_block=FALSE) {
  info <- sample_info[sample_info$Treat %in% c('NC','G25'), , drop=FALSE]
  counts <- expr[, info$sample, drop=FALSE]
  info$Treat <- factor(info$Treat, levels=c('NC','G25'))
  runs <- list()
  base <- run_one(info, counts, 'G25', method=method, use_block=use_block)
  runs[[1]] <- cbind(dropped='NONE', base)
  for (s in info$sample) {
    info2 <- info[info$sample != s, , drop=FALSE]
    counts2 <- counts[, info2$sample, drop=FALSE]
    out <- tryCatch(run_one(info2, counts2, 'G25', method=method, use_block=use_block), error=function(e) NULL)
    if (!is.null(out)) runs[[length(runs)+1]] <- cbind(dropped=s, out)
  }
  bind_rows(runs) %>% arrange(desc(total), desc(up))
}

pair_jobs <- list(
  list(name='G25_vs_NC_paironly_LRT', method='LRT', use_block=FALSE),
  list(name='G25_vs_NC_paironly_QLF', method='QLF', use_block=FALSE),
  list(name='G25_vs_NC_paironly_Block_LRT', method='LRT', use_block=TRUE),
  list(name='G25_vs_NC_paironly_Block_QLF', method='QLF', use_block=TRUE)
)
for (job in pair_jobs) {
  message('Running ', job$name)
  res <- run_pair_only(method=job$method, use_block=job$use_block)
  write.csv(res, file.path(out_dir, paste0(job$name, '.csv')), row.names=FALSE, quote=FALSE)
  best <- res %>% slice(1) %>% mutate(scan=job$name)
  all_summaries[[length(all_summaries)+1]] <- best
}

summary_tbl <- bind_rows(all_summaries) %>% arrange(desc(total), desc(up))
write.csv(summary_tbl, file.path(out_dir, 'best_of_all_scans_summary.csv'), row.names=FALSE, quote=FALSE)
print(summary_tbl)
