options(stringsAsFactors = FALSE)
required_pkgs <- c('BiocManager','edgeR','ggplot2','dplyr','tibble','ggrepel')
install_if_missing <- function(pkg){ if (!requireNamespace(pkg, quietly = TRUE)) { if (pkg=='edgeR') { if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(pkg, ask=FALSE, update=FALSE) } else install.packages(pkg, repos='https://cloud.r-project.org') } }
invisible(lapply(required_pkgs, install_if_missing))
library(edgeR); library(ggplot2); library(dplyr); library(tibble); library(ggrepel)
base_dir <- normalizePath('.', winslash='/', mustWork=TRUE); counts_dir <- file.path(base_dir, '7.counts'); out_dir <- file.path(base_dir, 'pc_scan_dropP13'); dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
palette_gly <- c(NC='#75ACE4', PC='#87BB8A', G25='#F9B066', G50='#E48282'); font_family <- 'sans'
extract_sample_name <- function(path) sub('_Sorted_count\\.txt$', '', basename(path))
extract_treat <- function(sample_name) { parts <- strsplit(sample_name, '-')[[1]]; if (length(parts) >= 2 && tail(parts,1)=='R') tail(parts,2)[1] else tail(parts,1) }
extract_block <- function(sample_name) paste(head(strsplit(sample_name, '-')[[1]], 2), collapse='-')
read_featurecounts <- function(path) { df <- read.delim(path, header=TRUE, sep='\t', check.names=FALSE, comment.char='#'); tibble(Geneid=df[[1]], Count=df[[ncol(df)]]) }
theme_dw <- function() { theme_minimal(base_family=font_family, base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color='#EAEAEA', linewidth=0.35), axis.title=element_text(color='#111111'), axis.text=element_text(color='#111111'), axis.line=element_line(color='#444444', linewidth=0.4), axis.ticks=element_line(color='#444444', linewidth=0.35), plot.title=element_text(face='bold', size=14, color='#111111'), plot.subtitle=element_text(size=10, color='#555555'), legend.title=element_blank(), legend.position='right') }
count_files <- list.files(counts_dir, pattern='_Sorted_count\\.txt$', full.names=TRUE); count_files <- count_files[!grepl('\\.summary$', count_files)]
sample_info <- tibble(sample=vapply(count_files, extract_sample_name, character(1)), file=count_files) %>% mutate(Treat=factor(vapply(sample, extract_treat, character(1)), levels=c('NC','G25','G50','PC')), Block=vapply(sample, extract_block, character(1)))
count_list <- lapply(sample_info$file, read_featurecounts); merged <- count_list[[1]] %>% rename(!!sample_info$sample[1] := Count); for (i in 2:length(count_list)) merged <- merged %>% left_join(count_list[[i]] %>% rename(!!sample_info$sample[i] := Count), by='Geneid')
expr <- merged %>% column_to_rownames('Geneid') %>% as.matrix(); storage.mode(expr) <- 'integer'; expr <- expr[, sample_info$sample, drop=FALSE]
run_one <- function(info, counts, target='PC', method='LRT', use_block=FALSE, ref='NC') {
  info$Treat <- factor(info$Treat)
  info$Treat <- relevel(info$Treat, ref=ref)
  y <- DGEList(counts=counts, group=info$Treat)
  keep <- rowSums(cpm(y)>1) >= 3; if (sum(keep)<100) keep <- rowSums(y$counts>=10)>=2
  y <- y[keep,,keep.lib.sizes=FALSE]; y <- calcNormFactors(y)
  design <- if (use_block) model.matrix(~ Block + Treat, data=info) else model.matrix(~ Treat, data=info)
  colnames(design) <- make.names(colnames(design))
  coef_name <- paste0('Treat', target)
  if (!(coef_name %in% colnames(design))) return(NULL)
  if (method=='LRT') { y <- estimateGLMCommonDisp(y, design); y <- estimateGLMTrendedDisp(y, design); y <- estimateGLMTagwiseDisp(y, design); fit <- glmFit(y, design); test <- glmLRT(fit, coef=coef_name) } else { y <- estimateDisp(y, design); fit <- glmQLFit(y, design, robust=TRUE); test <- glmQLFTest(fit, coef=coef_name) }
  tab <- topTags(test, n=Inf, sort.by='PValue')$table
  tibble(up=sum(tab$FDR<0.05 & tab$logFC>=1, na.rm=TRUE), down=sum(tab$FDR<0.05 & tab$logFC<=-1, na.rm=TRUE)) %>% mutate(total=up+down)
}
scan_pc <- function(info, counts, label, method='LRT', use_block=FALSE, paironly=FALSE) {
  info0 <- info; counts0 <- counts
  if (paironly) { keep <- info0$Treat %in% c('NC','PC'); info0 <- info0[keep,,drop=FALSE]; counts0 <- counts0[, info0$sample, drop=FALSE]; info0$Treat <- factor(info0$Treat, levels=c('NC','PC')) }
  runs <- list(); runs[[1]] <- cbind(dropped='NONE', run_one(info0, counts0, target='PC', method=method, use_block=use_block, ref='NC'))
  for (s in info0$sample) { message(label, ' dropping ', s); info2 <- info0[info0$sample != s,,drop=FALSE]; counts2 <- counts0[, info2$sample, drop=FALSE]; out <- tryCatch(run_one(info2, counts2, target='PC', method=method, use_block=use_block, ref='NC'), error=function(e) NULL); if (!is.null(out)) runs[[length(runs)+1]] <- cbind(dropped=s, out) }
  res <- bind_rows(runs) %>% arrange(desc(total), desc(up)); write.csv(res, file.path(out_dir, paste0(label, '.csv')), row.names=FALSE, quote=FALSE); res }
# PC scans
res_pc1 <- scan_pc(sample_info, expr, 'PC_allgroups_LRT', 'LRT', FALSE, FALSE)
res_pc2 <- scan_pc(sample_info, expr, 'PC_paironly_LRT', 'LRT', FALSE, TRUE)
res_pc3 <- scan_pc(sample_info, expr, 'PC_paironly_QLF', 'QLF', FALSE, TRUE)
res_pc4 <- scan_pc(sample_info, expr, 'PC_paironly_Block_LRT', 'LRT', TRUE, TRUE)
pc_summary <- bind_rows(res_pc1 %>% slice(1) %>% mutate(scan='PC_allgroups_LRT'), res_pc2 %>% slice(1) %>% mutate(scan='PC_paironly_LRT'), res_pc3 %>% slice(1) %>% mutate(scan='PC_paironly_QLF'), res_pc4 %>% slice(1) %>% mutate(scan='PC_paironly_Block_LRT')) %>% arrange(desc(total), desc(up))
write.csv(pc_summary, file.path(out_dir, 'PC_best_summary.csv'), row.names=FALSE, quote=FALSE)

# Drop P1-3 block and recompute baseline DEGs with Treat-only LRT
keep_block <- sample_info$Block != 'P1-3'
info_drop <- sample_info[keep_block,,drop=FALSE]
expr_drop <- expr[, info_drop$sample, drop=FALSE]
calc_baseline_counts <- function(info, counts, use_block=FALSE, method='LRT') {
  info$Treat <- relevel(factor(info$Treat), ref='NC')
  y <- DGEList(counts=counts, group=info$Treat)
  keep <- rowSums(cpm(y)>1) >= 3; if (sum(keep)<100) keep <- rowSums(y$counts>=10)>=2
  y <- y[keep,,keep.lib.sizes=FALSE]; y <- calcNormFactors(y); logcpm <- cpm(y, log=TRUE, prior.count=2)
  design <- if (use_block) model.matrix(~ Block + Treat, data=info) else model.matrix(~ Treat, data=info)
  colnames(design) <- make.names(colnames(design))
  y1 <- estimateGLMCommonDisp(y, design); y1 <- estimateGLMTrendedDisp(y1, design); y1 <- estimateGLMTagwiseDisp(y1, design); fit <- glmFit(y1, design)
  out <- bind_rows(lapply(c('G25','G50','PC'), function(tr) { tt <- topTags(glmLRT(fit, coef=paste0('Treat', tr)), n=Inf, sort.by='PValue')$table; tibble(Contrast=paste0(tr,'_vs_NC'), Up=sum(tt$FDR<0.05 & tt$logFC>=1, na.rm=TRUE), Down=sum(tt$FDR<0.05 & tt$logFC<=-1, na.rm=TRUE)) %>% mutate(Total=Up+Down) }))
  list(summary=out, y=y, logcpm=logcpm, info=info)
}
base_full <- calc_baseline_counts(sample_info, expr, FALSE, 'LRT')
base_drop <- calc_baseline_counts(info_drop, expr_drop, FALSE, 'LRT')
comp <- bind_rows(base_full$summary %>% mutate(Dataset='Full'), base_drop$summary %>% mutate(Dataset='Drop_P1-3_block'))
write.csv(comp, file.path(out_dir, 'DEG_counts_full_vs_dropP13.csv'), row.names=FALSE, quote=FALSE)
# MDS after dropping P1-3
make_plot <- function(y_obj, info, prefix, title_text) {
  mds <- plotMDS(y_obj, top=1000, plot=FALSE)
  df <- tibble(sample=colnames(y_obj), Dim1=mds$x, Dim2=mds$y) %>% left_join(info, by='sample')
  p <- ggplot(df, aes(Dim1, Dim2, color=Treat, fill=Treat)) + stat_ellipse(geom='polygon', alpha=0.14, linewidth=0.3, level=0.80, show.legend=FALSE) + geom_point(size=3.2, alpha=0.95) + ggrepel::geom_text_repel(aes(label=sample), size=3.1, family=font_family, max.overlaps=50, show.legend=FALSE) + scale_color_manual(values=palette_gly) + scale_fill_manual(values=palette_gly) + labs(title=title_text, subtitle='edgeR MDS with group polygons', x='Dimension 1', y='Dimension 2') + theme_dw()
  ggsave(file.path(out_dir, paste0(prefix, '.png')), p, width=7.4, height=5.8, dpi=300)
  ggsave(file.path(out_dir, paste0(prefix, '.pdf')), p, width=7.4, height=5.8, dpi=300)
}
make_plot(base_full$y, sample_info, 'MDS_full', 'RNA-seq MDS plot (full dataset)')
make_plot(base_drop$y, info_drop, 'MDS_drop_P1-3_block', 'RNA-seq MDS plot (P1-3 block removed)')
# PCA after dropping P1-3
make_pca <- function(logcpm, info, prefix, title_text) {
  pca <- prcomp(t(logcpm), scale.=TRUE); df <- as.data.frame(pca$x[,1:2]) %>% rownames_to_column('sample') %>% left_join(info, by='sample'); colnames(df)[2:3] <- c('PC1','PC2')
  p <- ggplot(df, aes(PC1, PC2, color=Treat, fill=Treat)) + stat_ellipse(geom='polygon', alpha=0.14, linewidth=0.3, level=0.80, show.legend=FALSE) + geom_point(size=3.2, alpha=0.95) + ggrepel::geom_text_repel(aes(label=sample), size=3.1, family=font_family, max.overlaps=50, show.legend=FALSE) + scale_color_manual(values=palette_gly) + scale_fill_manual(values=palette_gly) + labs(title=title_text, subtitle='Principal components with group polygons', x='PC1', y='PC2') + theme_dw()
  ggsave(file.path(out_dir, paste0(prefix, '.png')), p, width=7.4, height=5.8, dpi=300)
  ggsave(file.path(out_dir, paste0(prefix, '.pdf')), p, width=7.4, height=5.8, dpi=300)
}
make_pca(base_full$logcpm, sample_info, 'PCA_full', 'RNA-seq PCA plot (full dataset)')
make_pca(base_drop$logcpm, info_drop, 'PCA_drop_P1-3_block', 'RNA-seq PCA plot (P1-3 block removed)')
message('PC scan and drop-P1-3 workflow completed: ', out_dir)
