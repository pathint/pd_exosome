# ============================================================
# 04_MultiOmics.R —— FJMU miRNA + mRNA + lncRNA + 蛋白质组四层整合
# ============================================================

library(tidyverse)
library(limma)
library(edgeR)
library(mixOmics)
library(WGCNA)
library(igraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(reticulate)

use_python("/home/data/software/python/3.7.7/bin/python3", required = TRUE)

py_config()

py_run_string("import mofapy2; print(mofapy2.__version__)")
py_run_string("
import scipy
import numpy as np
if not hasattr(scipy, 'random'):
    scipy.random = np.random
")

library(MOFA2)

set.seed(42)
proc <- readRDS("data/processed/processed_data.rds")

s1_meta <- proc$s1_meta
s1_meta$Group <- factor(s1_meta$Group, levels = c("HC", "PD"))

# ================================================================
# 5.0 分层样本对齐（蛋白组少1个样本）
# ================================================================

common_3omics <- Reduce(intersect, list(
  rownames(s1_meta),
  colnames(proc$s1_mirna_raw),
  colnames(proc$s1_mrna_raw),
  colnames(proc$s1_lncrna_raw)
))

common_4omics <- intersect(common_3omics, colnames(proc$s1_protein))

cat("\n========== Sample Alignment ==========\n")
cat("3-omics (miRNA+mRNA+lncRNA):", length(common_3omics),
    "(PD:", sum(s1_meta[common_3omics, "Group"] == "PD"),
    "HC:", sum(s1_meta[common_3omics, "Group"] == "HC"), ")\n")
cat("4-omics (+protein):          ", length(common_4omics),
    "(PD:", sum(s1_meta[common_4omics, "Group"] == "PD"),
    "HC:", sum(s1_meta[common_4omics, "Group"] == "HC"), ")\n")

if (length(common_3omics) != length(common_4omics)) {
  missing_prot <- setdiff(common_3omics, common_4omics)
  cat("  Samples missing from protein:", paste(missing_prot, collapse = ", "), "\n")
}

# --- 3-omics 对齐（miRNA / mRNA / lncRNA DE 用） ---
s1_meta_3o   <- s1_meta[common_3omics, ]
s1_mirna_3o  <- proc$s1_mirna_raw[, common_3omics]
s1_mrna_3o   <- proc$s1_mrna_raw[, common_3omics]
s1_lncrna_3o <- proc$s1_lncrna_raw[, common_3omics]

stopifnot(nrow(s1_meta_3o) == ncol(s1_mirna_3o))
stopifnot(nrow(s1_meta_3o) == ncol(s1_mrna_3o))
stopifnot(nrow(s1_meta_3o) == ncol(s1_lncrna_3o))

design_3o <- model.matrix(~ Group, data = s1_meta_3o)
cat("Design 3-omics:", nrow(design_3o), "x", ncol(design_3o), "\n")

# --- 4-omics 对齐（蛋白 DE / MOFA / DIABLO 用） ---
s1_meta_4o    <- s1_meta[common_4omics, ]
s1_mirna_4o   <- proc$s1_mirna_raw[, common_4omics]
s1_mrna_4o    <- proc$s1_mrna_raw[, common_4omics]
s1_lncrna_4o  <- proc$s1_lncrna_raw[, common_4omics]
s1_protein_4o <- proc$s1_protein[, common_4omics]

stopifnot(nrow(s1_meta_4o) == ncol(s1_mirna_4o))
stopifnot(nrow(s1_meta_4o) == ncol(s1_protein_4o))

design_4o <- model.matrix(~ Group, data = s1_meta_4o)
cat("Design 4-omics:", nrow(design_4o), "x", ncol(design_4o), "\n")

# ================================================================
# 5.1 差异分析（统一 limma 框架，阈值适配小样本）
# ================================================================

cat("\n========== FJMU Differential Analysis ==========\n")

# --- 通用 limma 差异分析函数 ---
run_limma_de <- function(counts, design, feature_name = "feature",
                          fdr = 0.1, lfc = 0.5, is_count = TRUE) {
  if (is_count) {
    dge <- DGEList(counts = round(counts))
    keep <- filterByExpr(dge, design = design, min.count = 1)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge) == 0) {
      cat("  WARNING: no features passed filter, relaxing filter\n")
      keep <- rowSums(counts >= 1) >= 2
      dge <- DGEList(counts = round(counts[keep, ]))
      dge <- calcNormFactors(dge, method = "TMM")
    } else {
      dge <- calcNormFactors(dge, method = "TMM")
    }
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
  } else {
    fit <- lmFit(counts, design)
  }

  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf) %>%
    rownames_to_column(feature_name) %>%
    mutate(
      significance = case_when(
        adj.P.Val < fdr & logFC >  lfc ~ "Up",
        adj.P.Val < fdr & logFC < -lfc ~ "Down",
        TRUE ~ "NS"
      )
    )

  cat(sprintf("  %s: %d tested, %d Up, %d Down (FDR<%.2f, |LFC|>%.1f)\n",
              feature_name, nrow(tt),
              sum(tt$significance == "Up"),
              sum(tt$significance == "Down"),
              fdr, lfc))
  return(tt)
}

# --- miRNA（3-omics，13 samples）---
cat("\n--- miRNA ---\n")
res_mirna <- run_limma_de(s1_mirna_3o, design_3o, "miRNA",
                           fdr = 0.1, lfc = 0.5, is_count = TRUE)

# --- mRNA（3-omics，13 samples）---
cat("\n--- mRNA ---\n")
res_mrna <- run_limma_de(s1_mrna_3o, design_3o, "gene",
                          fdr = 0.1, lfc = 0.5, is_count = TRUE)
if (sum(res_mrna$significance != "NS") == 0) {
  cat("  No mRNA DEGs at FDR<0.1 & |LFC|>0.5, trying FDR<0.2 & |LFC|>0\n")
  res_mrna <- run_limma_de(s1_mrna_3o, design_3o, "gene",
                            fdr = 0.2, lfc = 0, is_count = TRUE)
}

# --- lncRNA（3-omics，13 samples）---
cat("\n--- lncRNA ---\n")
res_lnc <- run_limma_de(s1_lncrna_3o, design_3o, "lncRNA",
                         fdr = 0.1, lfc = 0.5, is_count = TRUE)

# --- 蛋白质组（4-omics，12 samples）---
cat("\n--- Protein ---\n")
prot_mat <- as.matrix(s1_protein_4o)

detect_rate <- rowSums(!is.na(prot_mat)) / ncol(prot_mat)
keep_prot <- detect_rate >= 0.7
cat(sprintf("  Protein: %d/%d passed detection (>70%%)\n",
            sum(keep_prot), nrow(prot_mat)))
prot_mat <- prot_mat[keep_prot, ]

if (any(is.na(prot_mat))) {
  min_val <- min(prot_mat, na.rm = TRUE)
  prot_mat[is.na(prot_mat)] <- min_val * 0.5
  cat("  Filled NA with half-minimum:", sprintf("%.2f", min_val * 0.5), "\n")
}

cat(sprintf("  Protein matrix: %d x %d, Design: %d rows\n",
            nrow(prot_mat), ncol(prot_mat), nrow(design_4o)))
stopifnot(ncol(prot_mat) == nrow(design_4o))

res_prot <- run_limma_de(prot_mat, design_4o, "protein",
                          fdr = 0.05, lfc = 0.5, is_count = FALSE)

# --- 保存 ---
write_csv(res_mirna, "results/04_MultiOmics/FJMU_DE_miRNA.csv")
write_csv(res_mrna,  "results/04_MultiOmics/FJMU_DE_mRNA.csv")
write_csv(res_lnc,   "results/04_MultiOmics/FJMU_DE_lncRNA.csv")
write_csv(res_prot,  "results/04_MultiOmics/FJMU_DE_protein.csv")

cat("\n========== DE Summary ==========\n")
cat(sprintf("  miRNA:  %d DE (FDR<0.1, |LFC|>0.5)\n", sum(res_mirna$significance != "NS")))
cat(sprintf("  mRNA:   %d DE\n", sum(res_mrna$significance != "NS")))
cat(sprintf("  lncRNA: %d DE\n", sum(res_lnc$significance != "NS")))
cat(sprintf("  Protein:%d DE (FDR<0.05, |LFC|>0.5)\n", sum(res_prot$significance != "NS")))

# ================================================================
# 5.2 Volcano Plots
# ================================================================

plot_de_volcano <- function(df, feature_col, title, fdr_line = 0.1) {
  df <- df %>% mutate(adj.P.Val = ifelse(is.na(adj.P.Val), 1, adj.P.Val))
  top <- df %>% filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(10)
  df$label <- ifelse(df[[feature_col]] %in% top[[feature_col]], df[[feature_col]], "")

  ggplot(df, aes(x = logFC, y = -log10(P.Value), color = significance)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(aes(label = label), size = 2, max.overlaps = 15, segment.alpha = 0.3) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    scale_color_manual(values = c("Up" = "#C0392B", "Down" = "#2980B9", "NS" = "#BDC3C7")) +
    labs(title = title, x = "logFC", y = "-log10(p)") +
    theme_bw(base_size = 11) + theme(legend.position = "none")
}

ggsave("results/04_MultiOmics/volcano_miRNA.pdf",
       plot_de_volcano(res_mirna, "miRNA", "FJMU miRNA"), width = 7, height = 6)
ggsave("results/04_MultiOmics/volcano_mRNA.pdf",
       plot_de_volcano(res_mrna, "gene", "FJMU mRNA"), width = 7, height = 6)
ggsave("results/04_MultiOmics/volcano_lncRNA.pdf",
       plot_de_volcano(res_lnc, "lncRNA", "FJMU lncRNA"), width = 7, height = 6)
ggsave("results/04_MultiOmics/volcano_protein.pdf",
       plot_de_volcano(res_prot, "protein", "FJMU Protein"), width = 7, height = 6)

# ================================================================
# 5.3 MOFA+ 四组学因子分析
# ================================================================

cat("\n========== MOFA+ ==========\n")

# voom 标准化各组学（均用 4-omics 共同样本）
dge_m <- DGEList(counts = round(s1_mirna_4o))
dge_m <- calcNormFactors(dge_m)
v_mirna_mofa <- voom(dge_m, design_4o, plot = FALSE)

dge_r <- DGEList(counts = round(s1_mrna_4o))
keep_r <- filterByExpr(dge_r, design = design_4o)
dge_r <- dge_r[keep_r, ]
dge_r <- calcNormFactors(dge_r)
v_mrna_mofa <- voom(dge_r, design_4o, plot = FALSE)

dge_l <- DGEList(counts = round(s1_lncrna_4o))
keep_l <- filterByExpr(dge_l, design = design_4o, min.count = 1)
if (sum(keep_l) < 10) {
  keep_l <- rowSums(round(s1_lncrna_4o) >= 1) >= 2
}
dge_l <- dge_l[keep_l, ]
dge_l <- calcNormFactors(dge_l)
v_lnc_mofa <- voom(dge_l, design_4o, plot = FALSE)

prot_scaled <- t(scale(t(prot_mat)))

mofa_data <- list(
  miRNA   = v_mirna_mofa$E,
  mRNA    = v_mrna_mofa$E,
  lncRNA  = v_lnc_mofa$E,
  Protein = prot_scaled
)

for (v in names(mofa_data)) {
  rv <- apply(mofa_data[[v]], 1, var, na.rm = TRUE)
  top_n <- min(3000, nrow(mofa_data[[v]]))
  mofa_data[[v]] <- mofa_data[[v]][order(rv, decreasing = TRUE)[1:top_n], ]
  cat(sprintf("  MOFA %s: %d features\n", v, nrow(mofa_data[[v]])))
}

mofa_obj <- create_mofa(mofa_data)
model_opts <- get_default_model_options(mofa_obj)
#model_opts$num_factors <- min(10, length(common_4omics) - 2)
model_opts$num_factors <- 3 

data_opts <- get_default_data_options(mofa_obj)

train_opts <- get_default_training_options(mofa_obj)
train_opts$maxiter <- 1000
train_opts$seed <- 43
train_opts$save_interrupted <- TRUE


mofa_obj <- prepare_mofa(
  object       = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts,
  #save_data    = FALSE
)


mofa_model <- run_mofa(mofa_obj,
                        outfile = "results/04_MultiOmics/mofa_model.hdf5",
                        save_data = TRUE,
						use_basilisk = FALSE)


reticulate::py_last_error()


# 方差解释
p_var <- plot_variance_explained(mofa_model, max_r2 = 10, plot_total = TRUE)
pdf("results/04_MultiOmics/MOFA_variance.pdf", width = 10, height = 6)
print(p_var)
dev.off()

# 因子分布
p_fac <- plot_factor(mofa_model, factors = 1:min(3, model_opts$num_factors),
                     color_by = "group", add_violin = TRUE, dot_size = 3)
pdf("results/04_MultiOmics/MOFA_factors.pdf", width = 10, height = 8)
print(p_fac)
dev.off()

# 各组学 top features
for (v in names(mofa_data)) {
  p_top <- plot_top_weights(mofa_model, view = v, factor = 1, nfeatures = 15)
  pdf(paste0("results/04_MultiOmics/MOFA_top_", v, ".pdf"), width = 8, height = 6)
  print(p_top)
  dev.off()
}

# ================================================================
# 5.4 DIABLO 四组学有监督整合
# ================================================================

cat("\n========== DIABLO ==========\n")

# 重新对齐 voom 输出（确保行名顺序与 design_4o 一致）
X_list <- list(
  miRNA   = t(v_mirna_mofa$E),
  mRNA    = t(v_mrna_mofa$E),
  lncRNA  = t(v_lnc_mofa$E),
  Protein = t(prot_scaled)
)

Y <- s1_meta_4o$Group
cat("DIABLO Y:", paste(table(Y), collapse = " "), "\n")

# 检查各 view 行数 = 样本数
for (nm in names(X_list)) {
  cat(sprintf("  %s: %d samples x %d features\n", nm, nrow(X_list[[nm]]), ncol(X_list[[nm]])))
  stopifnot(nrow(X_list[[nm]]) == length(Y))
}

# 设计矩阵
n_views <- length(X_list)
design_mat <- matrix(0.1, n_views, n_views,
                     dimnames = list(names(X_list), names(X_list)))
diag(design_mat) <- 0

# 初始模型
cat("  Building initial block.splsda...\n")
init_model <- block.splsda(X = X_list, Y = Y, ncomp = 2, design = design_mat)

# keepX 调优
test_kx <- list(
  miRNA   = c(3, 5, 10, min(20, ncol(X_list$miRNA))),
  mRNA    = c(5, 10, 20, min(50, ncol(X_list$mRNA))),
  lncRNA  = c(3, 5, 10, min(15, ncol(X_list$lncRNA))),
  Protein = c(3, 5, 10, min(15, ncol(X_list$Protein)))
)

opt_keepX <- list(
  miRNA   = min(10, ncol(X_list$miRNA)),
  mRNA    = min(20, ncol(X_list$mRNA)),
  lncRNA  = min(10, ncol(X_list$lncRNA)),
  Protein = min(10, ncol(X_list$Protein))
)

tryCatch({
  cat("  Tuning keepX...\n")
  tune_res <- tune.block.splsda(
    X = X_list, Y = Y, ncomp = 2,
    design = design_mat, test.keepX = test_kx,
    validation = "loo", nrepeat = 20,
    dist = "centroids.dist", progressBar = FALSE
  )
  opt_keepX <- tune_res$choice.keepX
  cat("  Optimal keepX:", paste(sapply(names(opt_keepX), function(n)
    paste0(n, "=", opt_keepX[[n]])), collapse = ", "), "\n")
}, error = function(e) {
  cat("  tune.block.splsda failed:", e$message, "\n")
  cat("  Using default keepX\n")
})

# 最终模型
cat("  Building final DIABLO model...\n")
diablo <- block.splsda(X = X_list, Y = Y, ncomp = 2,
                        keepX = opt_keepX, design = design_mat)

# 可视化: 样本空间
pdf("results/04_MultiOmics/DIABLO_sample.pdf", width = 10, height = 8)
plotIndiv(diablo, ind.names = TRUE, legend = TRUE,
          title = "DIABLO Sample Plot (comp 1-2)")
dev.off()

# 可视化: 相关圈图
pdf("results/04_MultiOmics/DIABLO_correlation.pdf", width = 8, height = 8)
plotVar(diablo, comp = 1:2, legend = TRUE,
        title = "DIABLO Correlation Circle")
dev.off()


# --- 选中特征（安全版） ---
selected <- selectVar(diablo, comp = 1)

diablo_feats <- data.frame(
  view = character(0), feature = character(0),
  value = numeric(0), stringsAsFactors = FALSE
)

for (v in names(selected)) {
  feat_info <- selected[[v]]
  # 安全检查：必须是 list 且包含 name 和 value
  if (is.list(feat_info) &&
      !is.null(feat_info$name) &&
      length(feat_info$name) > 0) {
    diablo_feats <- rbind(diablo_feats, data.frame(
      view    = v,
      feature = feat_info$name,
      value   = feat_info$value,
      stringsAsFactors = FALSE
    ))
  } else {
    cat("  No features selected for view:", v, "\n")
  }
}

write_csv(diablo_feats, "results/04_MultiOmics/DIABLO_features.csv")
cat("  DIABLO selected features:\n")
if (nrow(diablo_feats) > 0) print(diablo_feats) else cat("  (none)\n")

# --- 热图 ---
tryCatch({
  pdf("results/04_MultiOmics/DIABLO_heatmap.pdf", width = 14, height = 10)
  cimDiablo(diablo, margin = c(8, 20))
  dev.off()
}, error = function(e) {
  cat("  cimDiablo skipped:", e$message, "\n")
})

# ================================================================
# 5.5 WGCNA 共表达网络（mRNA 层面）
# ================================================================

cat("\n========== WGCNA ==========\n")

allowWGCNAThreads(4)

# 使用 4-omics 对齐的 mRNA（因 WGCNA 样本量差异1个影响不大）
# 或者用 3-omics 的 mRNA 获得多1个样本
# 此处统一用 4-omics，保持与下游网络分析一致

mrna_expr_wgcna <- t(v_mrna_mofa$E)
cat("WGCNA input:", nrow(mrna_expr_wgcna), "samples x", ncol(mrna_expr_wgcna), "genes\n")

# 软阈值选择
powers <- 1:20
sft <- pickSoftThreshold(mrna_expr_wgcna, powerVector = powers, verbose = 5)

# 绘制软阈值图
pdf("results/04_MultiOmics/WGCNA_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold", ylab = "Scale Free Topology Model Fit",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "blue", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

soft_pow <- sft$powerEstimate
if (is.na(soft_pow) || soft_pow < 3) soft_pow <- 6
cat("Soft threshold power:", soft_pow, "\n")

# 构建网络
net <- blockwiseModules(
  mrna_expr_wgcna,
  power = soft_pow,
  maxBlockSize = 20000,
  TOMType = "unsigned",
  minModuleSize = min(20, ncol(mrna_expr_wgcna) / 5),
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

mod_colors <- labels2colors(net$colors)
n_modules <- length(unique(mod_colors))
cat("Number of modules:", n_modules, "\n")

# 模块与 PD 关联
group_num <- as.numeric(s1_meta_4o$Group == "PD")

mod_trait <- map_dfr(unique(mod_colors), function(mc) {
  genes_in_mod <- which(mod_colors == mc)
  if (length(genes_in_mod) < 2) {
    me <- mrna_expr_wgcna[, genes_in_mod]
  } else {
    me <- rowMeans(mrna_expr_wgcna[, genes_in_mod, drop = FALSE])
  }
  ct <- cor.test(me, group_num, method = "spearman")
  data.frame(module = mc, n_genes = length(genes_in_mod),
             rho = ct$estimate, pvalue = ct$p.value,
             stringsAsFactors = FALSE)
})

mod_trait$fdr <- p.adjust(mod_trait$pvalue, method = "BH")
mod_trait <- mod_trait[order(mod_trait$pvalue), ]

cat("\nModules most associated with PD:\n")
print(head(mod_trait, 10))
write_csv(mod_trait, "results/04_MultiOmics/WGCNA_module_PD.csv")

# 模块树状图
pdf("results/04_MultiOmics/WGCNA_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  mod_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
dev.off()

# 模块与性状关系图
if (n_modules > 1) {
  mod_trait_mat <- matrix(mod_trait$rho, nrow = 1,
                          dimnames = list("PD vs HC", mod_trait$module))

  pdf("results/04_MultiOmics/WGCNA_module_trait.pdf", width = max(6, n_modules * 0.8), height = 3)
  Heatmap(mod_trait_mat,
          name = "Spearman rho",
          col = colorRamp2(c(-0.5, 0, 0.5), c("#3498DB", "white", "#E74C3C")),
          cluster_columns = TRUE, cluster_rows = FALSE,
          column_names_rot = 45,
          row_names_gp = gpar(fontsize = 11),
          column_names_gp = gpar(fontsize = 9),
          cell_fun = function(j, i, x, y, w, h, fill) {
            pval <- mod_trait$pvalue[j]
            stars <- ifelse(pval < 0.001, "***",
                     ifelse(pval < 0.01, "**",
                     ifelse(pval < 0.05, "*", "")))
            grid.text(stars, x, y, gp = gpar(fontsize = 10))
          },
          column_title = "WGCNA Module-PD Association")
  dev.off()
}

# ================================================================
# 5.6 miRNA-mRNA 靶标分析与相关性
# ================================================================

cat("\n========== miRNA-mRNA Correlation ==========\n")

sig_mir <- res_mirna %>% filter(significance != "NS") %>% pull(miRNA)
sig_mrna_genes <- res_mrna %>% filter(significance != "NS") %>% pull(gene)

# 在 3-omics 样本中计算相关性（多1个样本）
mir_expr_corr <- v_mirna_mofa$E[sig_mir[1:min(20, length(sig_mir))],
                                  , drop = FALSE]
mr_expr_corr  <- v_mrna_mofa$E[sig_mrna_genes[1:min(50, length(sig_mrna_genes))],
                                 , drop = FALSE]

# 共有样本中
corr_samples <- intersect(colnames(mir_expr_corr), colnames(mr_expr_corr))

if (length(sig_mir) > 0 & length(sig_mrna_genes) > 0 & length(corr_samples) >= 6) {
  pairs <- expand.grid(
    miRNA = sig_mir[1:min(20, length(sig_mir))],
    mRNA  = sig_mrna_genes[1:min(50, length(sig_mrna_genes))],
    stringsAsFactors = FALSE
  )

  corr_res <- pmap_dfr(pairs, function(miRNA, mRNA) {
    if (miRNA %in% rownames(mir_expr_corr) & mRNA %in% rownames(mr_expr_corr)) {
      tryCatch({
        ct <- cor.test(mir_expr_corr[miRNA, corr_samples],
                        mr_expr_corr[mRNA, corr_samples],
                        method = "spearman")
        data.frame(miRNA = miRNA, mRNA = mRNA,
                   rho = ct$estimate, pvalue = ct$p.value,
                   stringsAsFactors = FALSE)
      }, error = function(e) NULL)
    }
  })

  if (nrow(corr_res) > 0) {
    corr_res$fdr <- p.adjust(corr_res$pvalue, method = "BH")
    sig_neg <- corr_res %>% filter(fdr < 0.2, rho < 0)

    cat("Total pairs tested:", nrow(corr_res), "\n")
    cat("Significant negative correlations (FDR<0.2):", nrow(sig_neg), "\n")

    write_csv(corr_res, "results/04_MultiOmics/miRNA_mRNA_correlation.csv")

    # 相关性热图（top negative pairs）
    if (nrow(sig_neg) > 1) {
      top_neg <- sig_neg %>% arrange(rho) %>% head(30)

      mir_unique <- unique(top_neg$miRNA)
      mr_unique  <- unique(top_neg$mRNA)
      corr_mat <- matrix(NA, nrow = length(mir_unique), ncol = length(mr_unique),
                         dimnames = list(mir_unique, mr_unique))
      for (i in seq_len(nrow(top_neg))) {
        corr_mat[top_neg$miRNA[i], top_neg$mRNA[i]] <- top_neg$rho[i]
      }

      pdf("results/04_MultiOmics/miRNA_mRNA_corr_heatmap.pdf",
          width = max(8, length(mr_unique) * 0.4),
          height = max(5, length(mir_unique) * 0.4))
      Heatmap(corr_mat,
              name = "Spearman rho",
              col = colorRamp2(c(-0.8, 0, 0.8), c("#3498DB", "white", "#E74C3C")),
              na_col = "grey90",
              cluster_rows = TRUE, cluster_columns = TRUE,
              row_names_gp = gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 7),
              column_names_rot = 45,
              row_title = "miRNA", column_title = "Target mRNA",
              heatmap_legend_param = list(title = "rho"))
      dev.off()
    }
  }
} else {
  cat("Insufficient DE features or samples for correlation analysis\n")
}


# ================================================================
# 5.7 miRNA-mRNA 负调控验证（纯 base R 版本）
# ================================================================

cat("\n========== Negative Regulation Validation ==========\n")

if (length(sig_mir) > 0 & length(sig_mrna_genes) > 0) {

  # 提取 miRNA 方向
  mirna_dir <- res_mirna[res_mirna$significance != "NS", c("miRNA", "logFC")]
  colnames(mirna_dir) <- c("miRNA", "mirna_logFC")

  # 提取 mRNA 方向
  mrna_dir <- res_mrna[res_mrna$significance != "NS", c("gene", "logFC")]
  colnames(mrna_dir) <- c("mRNA", "mrna_logFC")

  # 全组合
  direction_pairs <- expand.grid(
    miRNA = mirna_dir$miRNA,
    mRNA  = mrna_dir$mRNA,
    stringsAsFactors = FALSE
  )

  # 合并（merge 代替 left_join）
  direction_pairs <- merge(direction_pairs, mirna_dir, by = "miRNA", all.x = TRUE)
  direction_pairs <- merge(direction_pairs, mrna_dir,  by = "mRNA",  all.x = TRUE)
  direction_pairs$opposite_direction <- sign(direction_pairs$mirna_logFC) !=
                                         sign(direction_pairs$mrna_logFC)

  cat("miRNA-mRNA direction pairs:", nrow(direction_pairs), "\n")
  cat("Opposite direction (expected negative regulation):",
      sum(direction_pairs$opposite_direction, na.rm = TRUE), "/",
      nrow(direction_pairs), "\n")
}

# ================================================================
# 5.8 四组学一致性通路概览
# ================================================================

cat("\n========== Cross-Omics DE Overlap ==========\n")

# 统计各层 DE 数量
de_summary <- data.frame(
  Layer = c("miRNA", "mRNA", "lncRNA", "Protein"),
  DE_count = c(
    sum(res_mirna$significance != "NS"),
    sum(res_mrna$significance != "NS"),
    sum(res_lnc$significance != "NS"),
    sum(res_prot$significance != "NS")
  ),
  Up = c(
    sum(res_mirna$significance == "Up"),
    sum(res_mrna$significance == "Up"),
    sum(res_lnc$significance == "Up"),
    sum(res_prot$significance == "Up")
  ),
  Down = c(
    sum(res_mirna$significance == "Down"),
    sum(res_mrna$significance == "Down"),
    sum(res_lnc$significance == "Down"),
    sum(res_prot$significance == "Down")
  ),
  stringsAsFactors = FALSE
)
print(de_summary)

# DE 数量柱状图
de_long <- de_summary %>%
  pivot_longer(cols = c(Up, Down), names_to = "direction", values_to = "count")

p_de_bar <- ggplot(de_long, aes(x = Layer, y = count, fill = direction)) +
  geom_col(position = "stack", width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB")) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 4) +
  labs(title = "DE Features by Omics Layer (FJMU)",
       x = "", y = "Number of DE Features") +
  theme_bw(base_size = 13)
ggsave("results/04_MultiOmics/DE_count_by_layer.pdf", p_de_bar, width = 7, height = 5)

# ================================================================
# 5.9 保存中间结果
# ================================================================

# 保存各组学 voom 标准化后的表达矩阵（供下游模块使用）
multiomics_processed <- list(
  # 3-omics 对齐
  common_3omics = common_3omics,
  s1_meta_3o = s1_meta_3o,
  design_3o = design_3o,
  # 4-omics 对齐
  common_4omics = common_4omics,
  s1_meta_4o = s1_meta_4o,
  design_4o = design_4o,
  # 标准化表达矩阵
  v_mirna = v_mirna_mofa,
  v_mrna = v_mrna_mofa,
  v_lnc = v_lnc_mofa,
  prot_mat = prot_mat,
  prot_scaled = prot_scaled,
  # DE 结果
  res_mirna = res_mirna,
  res_mrna = res_mrna,
  res_lnc = res_lnc,
  res_prot = res_prot,
  # WGCNA
  wgcna_colors = mod_colors,
  wgcna_net = net,
  # DIABLO
  diablo_model = diablo
)

saveRDS(multiomics_processed, "data/processed/multiomics_processed.rds")

cat("\n========== Multi-Omics Summary ==========\n")
cat(sprintf("  Samples: %d (3-omics) / %d (4-omics)\n",
            length(common_3omics), length(common_4omics)))
cat(sprintf("  DE miRNA: %d | mRNA: %d | lncRNA: %d | Protein: %d\n",
            sum(res_mirna$significance != "NS"),
            sum(res_mrna$significance != "NS"),
            sum(res_lnc$significance != "NS"),
            sum(res_prot$significance != "NS")))
cat(sprintf("  WGCNA modules: %d\n", n_modules))
cat(sprintf("  DIABLO features: %d\n", nrow(diablo_feats)))
cat("\nMulti-omics module complete.\n")
