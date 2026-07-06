# ============================================================
# 02_QC.R —— 质控、标准化、批次效应诊断与校正
# ============================================================

library(tidyverse)
library(edgeR)
library(limma)
library(sva)
library(pheatmap)
library(cowplot)
library(ggrepel)

set.seed(42)
proc <- readRDS("data/processed/proc_raw.rds")

# ================================================================
# 2.1 临床变量预处理
# ================================================================

clean_meta <- function(meta_df) {

  # 安全数值转换函数：非数值型统一变 NA
  safe_numeric <- function(x) {
    x <- as.character(x)
    x[x == "" | x == "-" | x == "NA" | is.na(x)] <- NA_character_
    suppressWarnings(as.numeric(x))
  }

  meta_df <- meta_df %>%
    mutate(
      Group     = factor(Group, levels = c("HC", "PD")),
      Dataset   = factor(Dataset),
      Gender    = factor(Gender, levels = c("M", "F")),
      Age       = safe_numeric(Age),
      Yahr      = safe_numeric(Yahr),
      UPDRS.III = safe_numeric(UPDRS.III),
      LEDD      = safe_numeric(LEDD),
      early.H.M  = safe_numeric(early.H.M),
      delayed.H.M = safe_numeric(delayed.H.M),
      MIBG.WR    = safe_numeric(MIBG.WR),
      DaTscan.R   = safe_numeric(DaTscan.R),
      DaTscan.L   = safe_numeric(DaTscan.L),
      DaTscan.Ave = safe_numeric(DaTscan.Ave),
      DaTscan.AI  = safe_numeric(DaTscan.AI...)
    )
  return(meta_df)
}

proc$g1_meta <- clean_meta(proc$g1_meta)
proc$g2_meta <- clean_meta(proc$g2_meta)
proc$s1_meta <- clean_meta(proc$s1_meta)

# --- 诊断：各变量缺失率 ---
cat("\n========== Clinical Data Availability ==========\n")
for (ds_name in c("g1_meta", "g2_meta", "s1_meta")) {
  m <- proc[[ds_name]]
  cat("\n---", ds_name, "(n=", nrow(m), ":",
      sum(m$Group == "PD"), "PD,", sum(m$Group == "HC"), "HC) ---\n")

  check_vars <- c("Age", "Gender", "Yahr", "UPDRS.III", "LEDD",
                   "early.H.M", "delayed.H.M", "MIBG.WR",
                   "DaTscan.R", "DaTscan.L", "DaTscan.Ave", "DaTscan.AI")
  for (v in check_vars) {
    if (v %in% colnames(m)) {
      n_total <- sum(!is.na(m[[v]]))
      n_pd <- sum(!is.na(m[[v]][m$Group == "PD"]))
      n_hc <- sum(!is.na(m[[v]][m$Group == "HC"]))
      cat(sprintf("  %-14s: total=%2d  PD=%2d  HC=%2d\n", v, n_total, n_pd, n_hc))
    }
  }
}

# ================================================================
# 2.2 各数据集独立质控
# ================================================================

run_qc <- function(counts, meta, dataset_name, outdir = "results/01_QC/") {
  cat("\n========== QC:", dataset_name, "==========\n")
  cat("Raw:", nrow(counts), "features x", ncol(counts), "samples\n")

  # 低表达过滤
  dge <- DGEList(counts = counts, samples = meta)
  keep <- filterByExpr(dge, group = meta$Group, min.count = 1, min.prop = 0.8)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("After filter:", nrow(dge$counts), "features\n")

  # TMM 归一化
  dge <- calcNormFactors(dge, method = "TMM")
  cpm_mat <- cpm(dge, log = TRUE, prior.count = 1)

  # --- PCA: Group ---
  pca_res <- prcomp(t(cpm_mat), scale. = TRUE)
  var_exp <- summary(pca_res)$importance[2, 1:4] * 100
  pca_df <- data.frame(
    PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], PC3 = pca_res$x[, 3],
    Group = meta$Group, sample_id = rownames(meta)
  )

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = 2) +
    labs(title = paste0(dataset_name, " — PCA by Group"),
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2])) +
    scale_color_manual(values = c("PD" = "#E74C3C", "HC" = "#3498DB")) +
    theme_bw(base_size = 13)
  ggsave(file.path(outdir, paste0(dataset_name, "_PCA_group.pdf")), p_pca, width = 7, height = 6)

  # --- PCA: Gender（仅当有有效 Gender 数据时）---
  has_gender <- sum(!is.na(meta$Gender)) > 2
  if (has_gender) {
    pca_df$Gender <- meta$Gender
    p_sex <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Gender, shape = Group)) +
      geom_point(size = 2.5, alpha = 0.8) +
      labs(title = paste0(dataset_name, " — PCA by Gender")) +
      theme_bw(base_size = 12)
    ggsave(file.path(outdir, paste0(dataset_name, "_PCA_gender.pdf")), p_sex, width = 7, height = 6)
  } else {
    cat("  Skipping PCA by Gender (insufficient data)\n")
  }

  # --- 表达分布箱线图 ---
  cpm_long <- as.data.frame(cpm_mat) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample_id", values_to = "log2CPM") %>%
    left_join(meta %>% select(ID, Group), by = c("sample_id" = "ID"))

  p_box <- ggplot(cpm_long, aes(x = sample_id, y = log2CPM, fill = Group)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.7) +
    scale_fill_manual(values = c("PD" = "#E74C3C", "HC" = "#3498DB")) +
    labs(title = paste0(dataset_name, " — Expression Distribution")) +
    theme_bw(base_size = 11) +
    #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # 纵向旋转
		          axis.ticks.x = element_line())  # 可选：保留刻度线
  ggsave(file.path(outdir, paste0(dataset_name, "_expr_dist.pdf")), p_box, width = 12, height = 5)

  # --- 样本相关性热图 ---
  cor_mat <- cor(cpm_mat, method = "spearman")
  anno_col <- data.frame(Group = meta$Group, row.names = rownames(meta))

  pdf(file.path(outdir, paste0(dataset_name, "_cor_heatmap.pdf")), width = 10, height = 9)
  pheatmap(cor_mat, annotation_col = anno_col,
           annotation_colors = list(Group = c(PD = "#E74C3C", HC = "#3498DB")),
           color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
           show_rownames = FALSE, show_colnames = FALSE,
           main = paste0(dataset_name, " — Sample Correlation (Spearman)"),
           clustering_method = "ward.D2")
  dev.off()

  # --- 检测饱和度 ---
  detect_rate <- colSums(cpm(dge) > 1) / nrow(dge$counts) * 100
  depth <- colSums(dge$counts)
  qc_df <- data.frame(depth = depth, detect_pct = detect_rate, Group = meta$Group)

  p_sat <- ggplot(qc_df, aes(x = depth, y = detect_pct, color = Group)) +
    geom_point(size = 2) +
    labs(title = paste0(dataset_name, " — Detection Saturation"),
         x = "Library Size", y = "% Detected (CPM > 1)") +
    scale_color_manual(values = c("PD" = "#E74C3C", "HC" = "#3498DB")) +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, paste0(dataset_name, "_saturation.pdf")), p_sat, width = 7, height = 5)

  list(dge = dge, cpm = cpm_mat, meta = meta, pca = pca_df)
}

# suppressWarnings 避免 filterByExpr 的非关键警告干扰输出
qc_g1 <- suppressWarnings(run_qc(proc$g1_counts, proc$g1_meta, "GSE269775"))
qc_g2 <- suppressWarnings(run_qc(proc$g2_counts, proc$g2_meta, "GSE269776"))
qc_s1 <- suppressWarnings(run_qc(proc$s1_mirna_counts, proc$s1_meta, "FJMU"))

# ================================================================
# 2.3 跨数据集批次效应诊断
# ================================================================

cat("\n========== Batch Effect Diagnostics ==========\n")

common_mirnas <- Reduce(intersect, list(
  rownames(qc_g1$cpm), rownames(qc_g2$cpm), rownames(qc_s1$cpm)
))
g1g2_common <- intersect(rownames(qc_g1$cpm), rownames(qc_g2$cpm))

cat("G1 miRNAs after filter:", nrow(qc_g1$cpm), "\n")
cat("G2 miRNAs after filter:", nrow(qc_g2$cpm), "\n")
cat("S1 miRNAs after filter:", nrow(qc_s1$cpm), "\n")
cat("G1 ∩ G2:", length(g1g2_common), "\n")
cat("G1 ∩ G2 ∩ S1:", length(common_mirnas), "\n")

# --- 合并矩阵（三数据集，用于全景 PCA）---
combined_cpm <- cbind(
  qc_g1$cpm[common_mirnas, ],
  qc_g2$cpm[common_mirnas, ],
  qc_s1$cpm[common_mirnas, ]
)
combined_meta <- bind_rows(
  qc_g1$meta %>% select(ID, Group, Dataset),
  qc_g2$meta %>% select(ID, Group, Dataset),
  qc_s1$meta %>% select(ID, Group, Dataset)
)
rownames(combined_meta) <- combined_meta$ID

# 确保行名顺序一致
stopifnot(all(colnames(combined_cpm) == rownames(combined_meta)))

# --- 校正前 PCA ---
pca_raw <- prcomp(t(combined_cpm), scale. = TRUE)
pca_raw_df <- data.frame(
  PC1 = pca_raw$x[, 1], PC2 = pca_raw$x[, 2],
  Group = combined_meta$Group, Dataset = combined_meta$Dataset
)
var_raw <- summary(pca_raw)$importance[2, 1:2] * 100

p_before <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, color = Dataset, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Before Batch Correction (All 3 Datasets)",
       x = sprintf("PC1 (%.1f%%)", var_raw[1]),
       y = sprintf("PC2 (%.1f%%)", var_raw[2])) +
  scale_color_manual(values = c("GSE269775" = "#E67E22", "GSE269776" = "#27AE60", "FJMU" = "#8E44AD")) +
  theme_bw(base_size = 13)

# ================================================================
# 2.4 批次效应校正（ComBat）
# ================================================================

cat("\n========== ComBat Batch Correction ==========\n")

# --- 2.4.1 G1 + G2 校正 ---
g1g2_cpm <- cbind(qc_g1$cpm[g1g2_common, ], qc_g2$cpm[g1g2_common, ])
g1g2_meta <- bind_rows(
  qc_g1$meta %>% select(ID, Group, Dataset, Gender, Age),
  qc_g2$meta %>% select(ID, Group, Dataset, Gender, Age)
)
rownames(g1g2_meta) <- g1g2_meta$ID

stopifnot(all(colnames(g1g2_cpm) == rownames(g1g2_meta)))

model_g1g2 <- model.matrix(~ Group, data = g1g2_meta)
cat("ComBat model matrix dim:", dim(model_g1g2), "\n")
cat("Expression matrix dim:", dim(g1g2_cpm), "\n")

# 确认维度一致
stopifnot(nrow(model_g1g2) == ncol(g1g2_cpm))

g1g2_combat <- ComBat(
  dat = g1g2_cpm,
  batch = g1g2_meta$Dataset,
  mod = model_g1g2,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 校正后 PCA
pca_g1g2 <- prcomp(t(g1g2_combat), scale. = TRUE)
pca_g1g2_df <- data.frame(
  PC1 = pca_g1g2$x[, 1], PC2 = pca_g1g2$x[, 2],
  Group = g1g2_meta$Group, Dataset = g1g2_meta$Dataset
)
var_g1g2 <- summary(pca_g1g2)$importance[2, 1:2] * 100

p_after <- ggplot(pca_g1g2_df, aes(x = PC1, y = PC2, color = Dataset, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(linetype = Group), level = 0.95) +
  labs(title = "After ComBat (G1 + G2)",
       x = sprintf("PC1 (%.1f%%)", var_g1g2[1]),
       y = sprintf("PC2 (%.1f%%)", var_g1g2[2])) +
  scale_color_manual(values = c("GSE269775" = "#E67E22", "GSE269776" = "#27AE60")) +
  theme_bw(base_size = 13)

# 对比图
p_batch_g1g2 <- plot_grid(p_before %+% {
    data = pca_raw_df %>% filter(Dataset != "FJMU")
  }, p_after, ncol = 2, labels = c("A", "B"))
ggsave("results/01_QC/batch_correction_G1G2.pdf", p_batch_g1g2, width = 14, height = 6)

# --- 2.4.2 三数据集校正 ---
# FJMU 样本量极小(13)，三数据集直接 ComBat 可能导致过度校正
# 但仍提供全景视角

model_full <- model.matrix(~ Group, data = combined_meta)
stopifnot(nrow(model_full) == ncol(combined_cpm))

combat_all <- ComBat(
  dat = combined_cpm,
  batch = combined_meta$Dataset,
  mod = model_full,
  par.prior = TRUE,
  prior.plots = FALSE
)

# 校正后 PCA（三数据集）
pca_all <- prcomp(t(combat_all), scale. = TRUE)
pca_all_df <- data.frame(
  PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2],
  Group = combined_meta$Group, Dataset = combined_meta$Dataset
)
var_all <- summary(pca_all)$importance[2, 1:2] * 100

p_all_after <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = Dataset, shape = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(linetype = Group), level = 0.95) +
  labs(title = "After ComBat (All 3 Datasets)",
       x = sprintf("PC1 (%.1f%%)", var_all[1]),
       y = sprintf("PC2 (%.1f%%)", var_all[2])) +
  scale_color_manual(values = c("GSE269775" = "#E67E22", "GSE269776" = "#27AE60", "FJMU" = "#8E44AD")) +
  theme_bw(base_size = 13)

# 三数据集校正前后对比
p_all_compare <- plot_grid(p_before, p_all_after, ncol = 2, labels = c("Before", "After"))
ggsave("results/01_QC/batch_correction_all3.pdf", p_all_compare, width = 14, height = 6)

# ================================================================
# 2.5 方差分解分析
# ================================================================

cat("\n========== Variance Decomposition ==========\n")

# 只用 Group 和 Dataset（不含 Gender/Age，因缺失太多）
var_decomp <- apply(combined_cpm, 1, function(y) {
  fit <- lm(y ~ Group + Dataset, data = combined_meta)
  a <- anova(fit)
  ss <- sum(a[, "Sum Sq"])
  data.frame(
    Group  = a["Group", "Sum Sq"] / ss,
    Batch  = a["Dataset", "Sum Sq"] / ss,
    Residual = a["Residuals", "Sum Sq"] / ss
  )
})
var_decomp_df <- do.call(rbind, var_decomp)

var_summary <- data.frame(
  Source = c("Group (PD/HC)", "Batch (Dataset)", "Residual"),
  Proportion = c(mean(var_decomp_df$Group), mean(var_decomp_df$Batch),
                  mean(var_decomp_df$Residual))
)

p_var <- ggplot(var_summary, aes(x = reorder(Source, -Proportion), y = Proportion, fill = Source)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion * 100)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#E74C3C", "#F39C12", "#95A5A6")) +
  labs(title = "Proportion of Variance Explained", x = "", y = "Mean Proportion") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  ylim(0, max(var_summary$Proportion) * 1.2)
ggsave("results/01_QC/variance_decomposition.pdf", p_var, width = 6, height = 5)

# ================================================================
# 2.6 G1+G2 校正后各数据集分离验证
# ================================================================

cat("\n========== Post-Correction Dataset Separation Check ==========\n")

# 校正后检查 G1 和 G2 的 PD 是否更可比
g1_post <- g1g2_combat[, g1g2_meta$Dataset == "GSE269775"]
g2_post <- g1g2_combat[, g1g2_meta$Dataset == "GSE269776"]

# 校正后 PD vs HC 在各数据集中的 PC1 分布
g1g2_pca_check <- data.frame(
  PC1 = pca_g1g2$x[, 1], PC2 = pca_g1g2$x[, 2],
  Group = g1g2_meta$Group, Dataset = g1g2_meta$Dataset
)

p_sep <- ggplot(g1g2_pca_check, aes(x = PC1, fill = interaction(Group, Dataset))) +
  geom_density(alpha = 0.5) +
  labs(title = "Post-ComBat: PC1 Distribution by Group × Dataset",
       fill = "Group × Dataset") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave("results/01_QC/post_combat_distribution.pdf", p_sep, width = 10, height = 6)

# ================================================================
# 2.7 保存处理后数据
# ================================================================

processed_data <- list(
  # 各数据集独立QC结果
  g1_dge = qc_g1$dge, g1_cpm = qc_g1$cpm, g1_meta = qc_g1$meta,
  g2_dge = qc_g2$dge, g2_cpm = qc_g2$cpm, g2_meta = qc_g2$meta,
  s1_dge = qc_s1$dge, s1_cpm = qc_s1$cpm, s1_meta = qc_s1$meta,
  # 跨数据集
  common_mirnas = common_mirnas,
  g1g2_common = g1g2_common,
  g1g2_cpm_combat = g1g2_combat,
  g1g2_meta = g1g2_meta,
  combined_cpm_combat = combat_all,
  combined_meta = combined_meta,
  # 原始计数（用于DESeq2）
  g1_counts_raw = proc$g1_counts,
  g2_counts_raw = proc$g2_counts,
  s1_mirna_raw = proc$s1_mirna_counts,
  # FJMU多组学
  s1_mrna_raw    = proc$s1_mrna_counts,
  s1_lncrna_raw  = proc$s1_lncrna_counts,
  s1_lncrna_info = proc$s1_lncrna_info,
  s1_protein     = proc$s1_protein,
  # 完整meta
  meta_all = proc$meta_all
)

saveRDS(processed_data, "data/processed/processed_data.rds")

cat("\n========== QC Summary ==========\n")
cat(sprintf("%-30s: %d miRNAs × %d samples\n", "GSE269775", nrow(qc_g1$cpm), ncol(qc_g1$cpm)))
cat(sprintf("%-30s: %d miRNAs × %d samples\n", "GSE269776", nrow(qc_g2$cpm), ncol(qc_g2$cpm)))
cat(sprintf("%-30s: %d miRNAs × %d samples\n", "FJMU", nrow(qc_s1$cpm), ncol(qc_s1$cpm)))
cat(sprintf("%-30s: %d\n", "G1 ∩ G2 common miRNAs", length(g1g2_common)))
cat(sprintf("%-30s: %d\n", "All 3 common miRNAs", length(common_mirnas)))
cat("\nQC module complete.\n")

write.csv(g1g2_common, file = "g1g2_common.csv", row.names = FALSE)
write.csv(common_mirnas, file = "common_mirnas.csv", row.names = FALSE)
