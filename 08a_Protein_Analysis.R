# ============================================================
# 08_Protein_Analysis.R（完整更新版）
# 双公共数据集 + FJMU独立验证 + 各数据集独立DEA + J1独立WGCNA
# 训练集：JPST003013 (n=120) + JPST003026 (n=80) → ComBat 合并
# 测试集：FJMU (独立)
# 富集分析见独立脚本 08b_Protein_Enrichment.R
# ============================================================

library(tidyverse)
library(limma)
library(sva)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)
library(pROC)
library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(pheatmap)
library(igraph)
library(ggrepel)

set.seed(42)
allowWGCNAThreads()

dir.create("results/08_Protein",          recursive = TRUE, showWarnings = FALSE)
dir.create("results/08_Protein/WGCNA",    showWarnings = FALSE)
dir.create("results/08_Protein/WGCNA_J1", showWarnings = FALSE)
dir.create("results/08_Protein/ML",       showWarnings = FALSE)
dir.create("results/08_Protein/Enrichment", showWarnings = FALSE)
dir.create("results/08_Protein/Integration", showWarnings = FALSE)

# ================================================================
# 8.1 数据读取（CSV 格式）
# ================================================================

cat("\n========== 8.1 Data Loading ==========\n")

j1_expr <- as.matrix(read.csv("data/protein/JPST003013_protein.csv",
                               row.names = 1, check.names = FALSE))
j1_meta <- read.csv("data/protein/meta_JPST003013.csv", row.names = 1,
                      check.names = FALSE)

j2_expr <- as.matrix(read.csv("data/protein/JPST003026_protein.csv",
                               row.names = 1, check.names = FALSE))
j2_meta <- read.csv("data/protein/meta_JPST003026.csv", row.names = 1,
                      check.names = FALSE)

fjmu_expr <- as.matrix(read.csv("data/protein/FJMU_protein.csv",
                                  row.names = 1, check.names = FALSE))
fjmu_meta <- read.csv("data/protein/FJMU_protein_meta.csv", row.names = 1,
                        check.names = FALSE)

cat("JPST003013:", nrow(j1_expr), "proteins x", ncol(j1_expr), "samples\n")
cat("JPST003026:", nrow(j2_expr), "proteins x", ncol(j2_expr), "samples\n")
cat("FJMU:       ", nrow(fjmu_expr), "proteins x", ncol(fjmu_expr), "samples\n")

clean_group <- function(g) {
  g <- trimws(as.character(g))
  ifelse(g %in% c("PDND", "PD-MCI", "PDD"), "PD",
  ifelse(g %in% c("MSA-C", "MSA-P"), "MSA",
  ifelse(g == "HC", "HC", g)))
}

j1_meta$Group   <- factor(clean_group(j1_meta$Group),   levels = c("HC", "PD", "MSA"))
j2_meta$Group   <- factor(clean_group(j2_meta$Group),   levels = c("HC", "PD", "MSA"))
fjmu_meta$Group <- factor(clean_group(fjmu_meta$Group), levels = c("HC", "PD", "MSA"))

cat("\nJPST003013 groups:\n"); print(table(j1_meta$Group))
cat("\nJPST003026 groups:\n"); print(table(j2_meta$Group))
cat("\nFJMU groups:\n");       print(table(fjmu_meta$Group))

# ================================================================
# 辅助函数
# ================================================================

preprocess_expr <- function(mat, min_detect = 0.5) {
  detected <- rowMeans(!is.na(mat) & mat > 0)
  mat <- mat[detected >= min_detect, ]
  min_val <- min(mat[mat > 0], na.rm = TRUE) / 2
  mat[is.na(mat) | mat <= 0] <- min_val
  mat <- log2(mat)
  return(mat)
}

run_limma_de <- function(mat, meta, group1, group2, p_cut = 0.05, fc_cut = 0.5) {
  keep <- meta$Group %in% c(group1, group2)
  expr_sub <- mat[, keep]
  grp <- droplevels(meta$Group[keep])
  cat(sprintf("  %s vs %s: %d + %d samples, %d proteins\n",
              group1, group2, sum(grp == group1), sum(grp == group2), nrow(expr_sub)))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)
  rownames(design) <- colnames(expr_sub)
  fit <- lmFit(expr_sub, design)
  fit <- contrasts.fit(fit, makeContrasts(contrasts = paste0(group1, "-", group2),
                                            levels = design))
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 1, number = Inf) %>%
    rownames_to_column("Protein") %>%
    mutate(
      significance = ifelse(adj.P.Val < p_cut & abs(logFC) > fc_cut, "DE", "NS"),
      direction = case_when(
        adj.P.Val < p_cut & logFC > fc_cut  ~ "Up",
        adj.P.Val < p_cut & logFC < -fc_cut ~ "Down",
        TRUE ~ "NS"
      )
    ) %>% arrange(adj.P.Val)
  return(tt)
}

# ================================================================
# 8.2 公共数据集预处理 + 三数据集 ComBat 合并
# ================================================================

cat("\n========== 8.2 Preprocessing & Batch Correction ==========\n")

j1_log2   <- preprocess_expr(j1_expr)
j2_log2   <- preprocess_expr(j2_expr)
fjmu_log2 <- preprocess_expr(fjmu_expr)

cat("After filtering: J1 =", nrow(j1_log2), ", J2 =", nrow(j2_log2),
    ", FJMU =", nrow(fjmu_log2), "\n")

shared_all <- Reduce(intersect, list(rownames(j1_log2),
                                      rownames(j2_log2),
                                      rownames(fjmu_log2)))
cat("Shared proteins (all 3):", length(shared_all), "\n")

j1_shared   <- j1_log2[shared_all, ]
j2_shared   <- j2_log2[shared_all, ]
fjmu_shared <- fjmu_log2[shared_all, ]

#/ write.csv(j1_shared, file = "JPST003013_protein_common.csv")
#/ write.csv(j2_shared, file = "JPST003026_protein_common.csv")
#/ write.csv(fjmu_shared, file = "FJMU_protein_common.csv")

shared_j1j2 <- Reduce(intersect, list(rownames(j1_log2),
                                      rownames(j2_log2)
                                     ))
cat("Shared proteins (j1 and j2):", length(shared_j1j2), "\n")

#/ j1_g1g2  <- j1_log2[shared_j1j2, ]
#/ j2_g1g2  <- j2_log2[shared_j1j2, ]
#/ write.csv(j1_g1g2, file = "JPST003013_protein_g1g2.csv")
#/ write.csv(j2_g1g2, file = "JPST003026_protein_g1g2.csv")

combined_raw <- cbind(j1_shared, j2_shared, fjmu_shared)
batch <- c(rep("JPST003013", ncol(j1_shared)),
           rep("JPST003026", ncol(j2_shared)),
           rep("FJMU",       ncol(fjmu_shared)))

cat("Running ComBat (3 datasets)...\n")
combined_combat <- ComBat(dat = combined_raw, batch = batch,
                            mod = NULL, par.prior = TRUE)

combined_norm <- normalizeBetweenArrays(combined_combat, method = "quantile")

n_j1   <- ncol(j1_shared)
n_j2   <- ncol(j2_shared)
n_fjmu <- ncol(fjmu_shared)

train_mat <- combined_norm[, 1:(n_j1 + n_j2)]
test_mat  <- combined_norm[, (n_j1 + n_j2 + 1):ncol(combined_norm)]

combined_meta <- rbind(
  data.frame(j1_meta, Dataset = "JPST003013", stringsAsFactors = FALSE),
  data.frame(j2_meta, Dataset = "JPST003026", stringsAsFactors = FALSE)
)

cat("Training:", nrow(train_mat), "proteins x", ncol(train_mat), "samples\n")
cat("Testing: ", nrow(test_mat),  "proteins x", ncol(test_mat),  "samples\n")

# ================================================================
# 8.3 QC
# ================================================================

cat("\n========== 8.3 Quality Control ==========\n")

# PCA：训练集 + FJMU 合并看批次效应
qc_combined <- cbind(train_mat, test_mat)
qc_meta <- data.frame(
  Dataset = c(combined_meta$Dataset, rep("FJMU", ncol(test_mat))),
  Group   = c(as.character(combined_meta$Group), as.character(fjmu_meta$Group)),
  stringsAsFactors = FALSE
)

pca_all <- prcomp(t(qc_combined), scale. = TRUE)
ve <- summary(pca_all)$importance[2, 1:4] * 100

pca_df <- data.frame(PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2], qc_meta)

p_pca_ds <- ggplot(pca_df, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = c("JPST003013" = "#E67E22",
                                 "JPST003026" = "#27AE60",
                                 "FJMU"       = "#8E44AD")) +
  labs(title = "PCA: All Datasets (After ComBat)",
       x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_bw(base_size = 13)
ggsave("results/08_Protein/QC_PCA_all_datasets.pdf", p_pca_ds, width = 9, height = 7)

# 训练集 PCA
pca_train <- prcomp(t(train_mat), scale. = TRUE)
ve_tr <- summary(pca_train)$importance[2, 1:2] * 100
pca_tr_df <- data.frame(PC1 = pca_train$x[, 1], PC2 = pca_train$x[, 2],
                          Group = combined_meta$Group)

p_pca_tr <- ggplot(pca_tr_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linetype = 2) +
  scale_color_manual(values = c("HC" = "#3498DB", "PD" = "#E74C3C", "MSA" = "#27AE60")) +
  labs(title = "Training Set PCA",
       x = sprintf("PC1 (%.1f%%)", ve_tr[1]),
       y = sprintf("PC2 (%.1f%%)", ve_tr[2])) +
  theme_bw(base_size = 13)
ggsave("results/08_Protein/QC_PCA_training.pdf", p_pca_tr, width = 8, height = 6)

# ================================================================
# 8.4 各数据集独立 DEA（未去批次，未取共有蛋白）
# ================================================================

cat("\n========== 8.4 Per-Dataset DEA (Raw, Independent) ==========\n")

# --- JPST003013 独立 DEA ---
cat("\n--- JPST003013 DEA ---\n")
j1_norm_raw <- normalizeBetweenArrays(j1_log2, method = "quantile")
de_j1 <- run_limma_de(j1_norm_raw, j1_meta, "PD", "HC")
cat(sprintf("  J1 DE proteins: %d\n", sum(de_j1$significance == "DE")))
write_csv(de_j1, "results/08_Protein/DE_JPST003013_independent.csv")

# --- JPST003026 独立 DEA ---
cat("\n--- JPST003026 DEA ---\n")
j2_norm_raw <- normalizeBetweenArrays(j2_log2, method = "quantile")
de_j2 <- run_limma_de(j2_norm_raw, j2_meta, "PD", "HC")
cat(sprintf("  J2 DE proteins: %d\n", sum(de_j2$significance == "DE")))
write_csv(de_j2, "results/08_Protein/DE_JPST003026_independent.csv")

# ================================================================
# 8.5 合并 DEA（PD vs HC，训练集 ComBat 后）
# ================================================================

cat("\n========== 8.5 Differential Expression (Combined, PD vs HC) ==========\n")

keep_ph <- combined_meta$Group %in% c("PD", "HC")
expr_ph <- train_mat[, keep_ph]
grp_ph  <- droplevels(combined_meta$Group[keep_ph])

cat("PD vs HC:", sum(grp_ph == "PD"), "PD +", sum(grp_ph == "HC"), "HC\n")

design_ph <- model.matrix(~ 0 + grp_ph)
colnames(design_ph) <- levels(grp_ph)
rownames(design_ph) <- colnames(expr_ph)

fit_ph <- lmFit(expr_ph, design_ph)
fit_ph <- contrasts.fit(fit_ph, makeContrasts(PD - HC, levels = design_ph))
fit_ph <- eBayes(fit_ph)

de_protein <- topTable(fit_ph, coef = 1, number = Inf) %>%
  rownames_to_column("Protein") %>%
  mutate(
    significance = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.5, "DE", "NS"),
    direction = case_when(
      adj.P.Val < 0.05 & logFC > 0.5  ~ "Up",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>% arrange(adj.P.Val)

cat(sprintf("Combined DE proteins: %d Up, %d Down\n",
            sum(de_protein$direction == "Up"),
            sum(de_protein$direction == "Down")))
write_csv(de_protein, "results/08_Protein/DE_proteins_PD_vs_HC.csv")

# 火山图
top_label <- de_protein %>% filter(significance == "DE") %>% head(15)
p_vol <- ggplot(de_protein, aes(logFC, -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_label, aes(label = Protein),
                   size = 3, max.overlaps = 20, show.legend = FALSE) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#3498DB", "NS" = "#BDC3C7")) +
  labs(title = "DE Proteins: PD vs HC (Combined)", x = "log2(FC)", y = "-log10(FDR)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave("results/08_Protein/volcano_PD_vs_HC.pdf", p_vol, width = 8, height = 6)

# ================================================================
# 8.6 PD vs MSA（鉴别诊断）
# ================================================================

cat("\n========== 8.6 PD vs MSA ==========\n")

keep_pm <- combined_meta$Group %in% c("PD", "MSA")
if (sum(keep_pm) > 10) {
  expr_pm <- train_mat[, keep_pm]
  grp_pm  <- droplevels(combined_meta$Group[keep_pm])
  design_pm <- model.matrix(~ 0 + grp_pm)
  colnames(design_pm) <- levels(grp_pm)
  rownames(design_pm) <- colnames(expr_pm)
  fit_pm <- lmFit(expr_pm, design_pm)
  fit_pm <- contrasts.fit(fit_pm, makeContrasts(PD - MSA, levels = design_pm))
  fit_pm <- eBayes(fit_pm)
  de_pd_msa <- topTable(fit_pm, coef = 1, number = Inf) %>%
    rownames_to_column("Protein") %>%
    mutate(significance = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.5, "DE", "NS"))
  cat("PD vs MSA DE proteins:", sum(de_pd_msa$significance == "DE"), "\n")
  write_csv(de_pd_msa, "results/08_Protein/DE_proteins_PD_vs_MSA.csv")
}

# ================================================================
# 8.7 DE 三路交集分析
# ================================================================

cat("\n========== 8.7 DE Intersection Analysis ==========\n")

de_j1_up     <- de_j1 %>% filter(direction == "Up") %>% pull(Protein)
de_j1_down   <- de_j1 %>% filter(direction == "Down") %>% pull(Protein)
de_j2_up     <- de_j2 %>% filter(direction == "Up") %>% pull(Protein)
de_j2_down   <- de_j2 %>% filter(direction == "Down") %>% pull(Protein)
de_comb_up   <- de_protein %>% filter(direction == "Up") %>% pull(Protein)
de_comb_down <- de_protein %>% filter(direction == "Down") %>% pull(Protein)

shared_up   <- Reduce(intersect, list(de_j1_up, de_j2_up, de_comb_up))
shared_down <- Reduce(intersect, list(de_j1_down, de_j2_down, de_comb_down))
cat(sprintf("Shared Up   (J1 ∩ J2 ∩ Combined): %d\n", length(shared_up)))
cat(sprintf("Shared Down (J1 ∩ J2 ∩ Combined): %d\n", length(shared_down)))

j1_j2_up     <- intersect(de_j1_up, de_j2_up)
j1_comb_up   <- intersect(de_j1_up, de_comb_up)
j2_comb_up   <- intersect(de_j2_up, de_comb_up)
j1_j2_down   <- intersect(de_j1_down, de_j2_down)
j1_comb_down <- intersect(de_j1_down, de_comb_down)
j2_comb_down <- intersect(de_j2_down, de_comb_down)

cat(sprintf("  Up:   J1∩J2=%d, J1∩Comb=%d, J2∩Comb=%d\n",
            length(j1_j2_up), length(j1_comb_up), length(j2_comb_up)))
cat(sprintf("  Down: J1∩J2=%d, J1∩Comb=%d, J2∩Comb=%d\n",
            length(j1_j2_down), length(j1_comb_down), length(j2_comb_down)))

intersect_df <- data.frame(
  Comparison = c("J1∩J2∩Combined (Up)", "J1∩J2∩Combined (Down)",
                  "J1∩J2 (Up)", "J1∩Combined (Up)", "J2∩Combined (Up)",
                  "J1∩J2 (Down)", "J1∩Combined (Down)", "J2∩Combined (Down)"),
  Count = c(length(shared_up), length(shared_down),
            length(j1_j2_up), length(j1_comb_up), length(j2_comb_up),
            length(j1_j2_down), length(j1_comb_down), length(j2_comb_down)),
  stringsAsFactors = FALSE
)

#write_csv(intersect_df, "results/08_Protein/DE_intersection_summary.csv")
#write_csv(data.frame(Protein = shared_up,   Direction = "Up"),
#          "results/08_Protein/DE_shared_up_all3.csv")
#write_csv(data.frame(Protein = shared_down, Direction = "Down"),
#          "results/08_Protein/DE_shared_down_all3.csv")

# UpSet 图
pdf("results/08_Protein/DE_upset_up.pdf", width = 8, height = 5)
UpSetR::upset(UpSetR::fromList(list(JPST003013 = de_j1_up,
                                      JPST003026 = de_j2_up,
                                      Combined   = de_comb_up)),
              sets = c("JPST003013", "JPST003026", "Combined"),
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              text.scale = 1.1, main.bar.color = "#E74C3C")
dev.off()

pdf("results/08_Protein/DE_upset_down.pdf", width = 8, height = 5)
UpSetR::upset(UpSetR::fromList(list(JPST003013 = de_j1_down,
                                      JPST003026 = de_j2_down,
                                      Combined   = de_comb_down)),
              sets = c("JPST003013", "JPST003026", "Combined"),
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              text.scale = 1.1, main.bar.color = "#3498DB")
dev.off()

# logFC 一致性散点图
fc_compare <- de_j1 %>%
  dplyr::select(Protein, logFC_j1 = logFC) %>%
  inner_join(de_j2 %>% dplyr::select(Protein, logFC_j2 = logFC), by = "Protein")

cor_j1_j2 <- cor.test(fc_compare$logFC_j1, fc_compare$logFC_j2, method = "spearman")

p_fc <- ggplot(fc_compare, aes(logFC_j1, logFC_j2)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_smooth(method = "lm", color = "#E74C3C", se = FALSE, linewidth = 0.8) +
  annotate("text", x = min(fc_compare$logFC_j1), y = max(fc_compare$logFC_j2),
           label = sprintf("Spearman rho = %.3f\np = %.2e",
                           cor_j1_j2$estimate, cor_j1_j2$p.value),
           hjust = 0, size = 3.5) +
  labs(title = "logFC Consistency: JPST003013 vs JPST003026",
       x = "logFC (JPST003013)", y = "logFC (JPST003026)") +
  theme_bw(base_size = 13)
ggsave("results/08_Protein/DE_logFC_consistency_J1_vs_J2.pdf", p_fc, width = 7, height = 6)

# ================================================================
# 8.8 WGCNA（合并数据集，训练集）
# ================================================================

cat("\n========== 8.8 WGCNA (Combined Training Set) ==========\n")

datExpr <- as.data.frame(t(expr_ph))
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

traitDF <- data.frame(Disease = ifelse(grp_ph[rownames(datExpr)] == "PD", 1, 0))
rownames(traitDF) <- rownames(datExpr)

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- ifelse(!is.na(sft$powerEstimate), sft$powerEstimate, 6)

p_sft <- ggplot(
  data.frame(Power = sft$fitIndices$Power,
             R2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq),
  aes(Power, R2)) +
  geom_line(color = "#2C3E50", linewidth = 1) +
  geom_point(size = 2, color = ifelse(sft$fitIndices$Power == softPower,
                                       "#E74C3C", "#2C3E50")) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey50") +
  annotate("text", x = softPower + 1, y = 0.5,
           label = sprintf("Selected: %d", softPower), color = "#E74C3C") +
  labs(title = "WGCNA Soft Threshold (Combined)",
       x = "Power", y = "Scale-Free Topology Model Fit") +
  theme_bw(base_size = 12)
ggsave("results/08_Protein/WGCNA/soft_threshold.pdf", p_sft, width = 7, height = 5)

net <- blockwiseModules(
  datExpr, power = softPower, TOMType = "unsigned",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, verbose = 3
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs <- net$MEs

cat("Modules:", length(unique(moduleColors)), "\n")
cat("Module sizes:\n")
print(sort(table(moduleColors), decreasing = TRUE))

# 模块-trait（手动对齐）
me_names <- colnames(MEs)
moduleTraitCor <- matrix(NA, length(me_names), 1, dimnames = list(me_names, "Disease"))
moduleTraitPvalue <- moduleTraitCor

for (i in seq_along(me_names)) {
  me_v <- as.numeric(MEs[, i]); names(me_v) <- rownames(MEs)
  tr_v <- as.numeric(traitDF$Disease); names(tr_v) <- rownames(traitDF)
  common <- intersect(names(me_v), names(tr_v))
  ok <- is.finite(me_v[common]) & is.finite(tr_v[common])
  if (sum(ok) >= 10 && sd(me_v[common][ok]) > 0 && sd(tr_v[common][ok]) > 0) {
    t <- cor.test(me_v[common][ok], tr_v[common][ok])
    moduleTraitCor[i, 1] <- t$estimate; moduleTraitPvalue[i, 1] <- t$p.value
  }
}

valid <- apply(moduleTraitCor, 1, function(x) all(is.finite(x)))
moduleTraitCor    <- as.data.frame(moduleTraitCor[valid, , drop = FALSE])
moduleTraitPvalue <- as.data.frame(moduleTraitPvalue[valid, , drop = FALSE])

if (nrow(moduleTraitCor) > 0) {
  pdf("results/08_Protein/WGCNA/module_trait_heatmap.pdf", width = 6, height = 8)
  pheatmap(as.matrix(moduleTraitCor), display_numbers = TRUE,
           number_format = "%.3f", cluster_rows = FALSE, cluster_cols = FALSE,
           main = "Module-Trait Correlation (PD vs HC)")
  dev.off()
}

sig_modules <- rownames(moduleTraitCor)[
  abs(moduleTraitCor[, 1]) > 0.3 & moduleTraitPvalue[, 1] < 0.05
]
if (length(sig_modules) == 0) {
  sig_modules <- rownames(moduleTraitCor)[
    abs(moduleTraitCor[, 1]) > 0.2 & moduleTraitPvalue[, 1] < 0.1
  ]
}
cat("Significant modules:", length(sig_modules), "\n")
for (m in sig_modules) cat(sprintf("  %s: cor=%.3f, p=%.4f\n",
                                    m, moduleTraitCor[m, 1], moduleTraitPvalue[m, 1]))

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
hub_list <- list()
hub_summary <- data.frame()

for (m in sig_modules) {
  mod_color <- sub("^ME", "", m)
  mod_genes <- colnames(datExpr)[moduleColors == mod_color]
  kme <- geneModuleMembership[mod_genes, m, drop = TRUE]
  names(kme) <- mod_genes; kme <- sort(kme[is.finite(kme)], decreasing = TRUE)
  hub_list[[m]] <- head(kme, 20)
  hub_summary <- rbind(hub_summary, data.frame(
    Module = m, Protein = names(hub_list[[m]]),
    kME = as.numeric(hub_list[[m]]), stringsAsFactors = FALSE))
}
write_csv(hub_summary, "results/08_Protein/WGCNA/hub_genes_all_modules.csv")

wgcna_results <- list(net = net, MEs = MEs, moduleColors = moduleColors,
                       moduleTraitCor = moduleTraitCor,
                       moduleTraitPvalue = moduleTraitPvalue,
                       hub_list = hub_list, sig_modules = sig_modules)
saveRDS(wgcna_results, "results/08_Protein/WGCNA/wgcna_results.rds")

# ================================================================
# 8.9 JPST003013 独立 WGCNA（未去批次，未取共有蛋白）
# ================================================================

cat("\n========== 8.9 WGCNA on JPST003013 (Independent) ==========\n")

keep_j1_ph <- j1_meta$Group %in% c("PD", "HC")
j1_ph_expr <- j1_norm_raw[, keep_j1_ph]
j1_ph_grp  <- droplevels(j1_meta$Group[keep_j1_ph])

cat("J1 PD+HC:", ncol(j1_ph_expr),
    "(PD:", sum(j1_ph_grp == "PD"), "HC:", sum(j1_ph_grp == "HC"), ")\n")

datExpr_j1 <- as.data.frame(t(j1_ph_expr))
gsg_j1 <- goodSamplesGenes(datExpr_j1, verbose = 3)
if (!gsg_j1$allOK) datExpr_j1 <- datExpr_j1[gsg_j1$goodSamples, gsg_j1$goodGenes]

traitDF_j1 <- data.frame(Disease = ifelse(j1_ph_grp[rownames(datExpr_j1)] == "PD", 1, 0))
rownames(traitDF_j1) <- rownames(datExpr_j1)

sft_j1 <- pickSoftThreshold(datExpr_j1, powerVector = c(1:20), verbose = 5)
softPower_j1 <- ifelse(!is.na(sft_j1$powerEstimate), sft_j1$powerEstimate, 6)
cat("J1 soft-thresholding power:", softPower_j1, "\n")

p_sft_j1 <- ggplot(
  data.frame(Power = sft_j1$fitIndices$Power,
             R2 = -sign(sft_j1$fitIndices$slope) * sft_j1$fitIndices$SFT.R.sq),
  aes(Power, R2)) +
  geom_line(color = "#2C3E50", linewidth = 1) +
  geom_point(size = 2, color = ifelse(sft_j1$fitIndices$Power == softPower_j1,
                                       "#E74C3C", "#2C3E50")) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey50") +
  labs(title = "WGCNA Soft Threshold (JPST003013)",
       x = "Power", y = "Scale-Free Topology Model Fit") +
  theme_bw(base_size = 12)
ggsave("results/08_Protein/WGCNA_J1/soft_threshold.pdf", p_sft_j1, width = 7, height = 5)

net_j1 <- blockwiseModules(
  datExpr_j1, power = softPower_j1, TOMType = "unsigned",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, verbose = 3
)

moduleLabels_j1 <- net_j1$colors
moduleColors_j1 <- labels2colors(moduleLabels_j1)
MEs_j1 <- net_j1$MEs

cat("J1 modules:", length(unique(moduleColors_j1)), "\n")
print(sort(table(moduleColors_j1), decreasing = TRUE))

# 模块-trait
me_names_j1 <- colnames(MEs_j1)
mtCor_j1 <- matrix(NA, length(me_names_j1), 1, dimnames = list(me_names_j1, "Disease"))
mtPval_j1 <- mtCor_j1

for (i in seq_along(me_names_j1)) {
  me_v <- as.numeric(MEs_j1[, i]); names(me_v) <- rownames(MEs_j1)
  tr_v <- as.numeric(traitDF_j1$Disease); names(tr_v) <- rownames(traitDF_j1)
  common <- intersect(names(me_v), names(tr_v))
  ok <- is.finite(me_v[common]) & is.finite(tr_v[common])
  if (sum(ok) >= 10 && sd(me_v[common][ok]) > 0 && sd(tr_v[common][ok]) > 0) {
    t <- cor.test(me_v[common][ok], tr_v[common][ok])
    mtCor_j1[i, 1] <- t$estimate; mtPval_j1[i, 1] <- t$p.value
  }
}

valid_j1 <- apply(mtCor_j1, 1, function(x) all(is.finite(x)))
mtCor_j1  <- as.data.frame(mtCor_j1[valid_j1, , drop = FALSE])
mtPval_j1 <- as.data.frame(mtPval_j1[valid_j1, , drop = FALSE])

if (nrow(mtCor_j1) > 0) {
  pdf("results/08_Protein/WGCNA_J1/module_trait_heatmap.pdf", width = 6, height = 8)
  pheatmap(as.matrix(mtCor_j1), display_numbers = TRUE,
           number_format = "%.3f", cluster_rows = FALSE, cluster_cols = FALSE,
           main = "Module-Trait Correlation (J1: PD vs HC)")
  dev.off()
}

sig_mods_j1 <- rownames(mtCor_j1)[
  abs(mtCor_j1[, 1]) > 0.3 & mtPval_j1[, 1] < 0.05
]
if (length(sig_mods_j1) == 0) {
  sig_mods_j1 <- rownames(mtCor_j1)[
    abs(mtCor_j1[, 1]) > 0.2 & mtPval_j1[, 1] < 0.1
  ]
}
cat("\nJ1 significant modules:\n")
for (m in sig_mods_j1) cat(sprintf("  %s: cor=%.3f, p=%.4f\n", m, mtCor_j1[m, 1], mtPval_j1[m, 1]))

# Hub genes
geneMM_j1 <- as.data.frame(cor(datExpr_j1, MEs_j1, use = "p"))
hub_list_j1 <- list()
hub_sum_j1 <- data.frame()

for (m in sig_mods_j1) {
  mod_color <- sub("^ME", "", m)
  mod_genes <- colnames(datExpr_j1)[moduleColors_j1 == mod_color]
  kme <- geneMM_j1[mod_genes, m, drop = TRUE]
  names(kme) <- mod_genes; kme <- sort(kme[is.finite(kme)], decreasing = TRUE)
  hub_list_j1[[m]] <- head(kme, 20)
  hub_sum_j1 <- rbind(hub_sum_j1, data.frame(
    Module = m, Protein = names(hub_list_j1[[m]]),
    kME = as.numeric(hub_list_j1[[m]]), stringsAsFactors = FALSE))
}
write_csv(hub_sum_j1, "results/08_Protein/WGCNA_J1/hub_genes_all_modules.csv")

# Eigengene boxplot
for (m in sig_mods_j1) {
  eig_df <- data.frame(
    Sample = rownames(MEs_j1),
    ME = as.numeric(MEs_j1[, m]),
    Group = ifelse(traitDF_j1$Disease[rownames(MEs_j1)] == 1, "PD", "HC"))

  p_eig <- ggplot(eig_df, aes(x = Group, y = ME, color = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    scale_color_manual(values = c("HC" = "#3498DB", "PD" = "#E74C3C")) +
    labs(title = sprintf("%s Eigengene (J1)", m), y = "Eigengene") +
    theme_bw(base_size = 13) + theme(legend.position = "none")

  fname <- gsub(" ", "_", paste0("results/08_Protein/WGCNA_J1/", m, "_eigengene.pdf"))
  ggsave(fname, p_eig, width = 5, height = 5)
}

# J1 模块基因导出
wgcna_j1_module_genes <- data.frame(
  Protein = colnames(datExpr_j1), Module = moduleColors_j1,
  stringsAsFactors = FALSE)
write_csv(wgcna_j1_module_genes, "results/08_Protein/WGCNA_J1/all_module_genes.csv")

# J1 vs Combined hub 对比
j1_hub_proteins <- unique(hub_sum_j1$Protein)
combined_hub_proteins <- unique(hub_summary$Protein)
shared_hubs <- intersect(j1_hub_proteins, combined_hub_proteins)
cat(sprintf("\nShared hub proteins (J1 ∩ Combined): %d\n", length(shared_hubs)))
write_csv(data.frame(Protein = shared_hubs),
          "results/08_Protein/WGCNA_J1/shared_hubs_J1_vs_combined.csv")

wgcna_j1_results <- list(
  net = net_j1, MEs = MEs_j1, moduleColors = moduleColors_j1,
  moduleLabels = moduleLabels_j1,
  moduleTraitCor = mtCor_j1, moduleTraitPvalue = mtPval_j1,
  hub_list = hub_list_j1, sig_modules = sig_mods_j1,
  datExpr = datExpr_j1, traitDF = traitDF_j1)
saveRDS(wgcna_j1_results, "results/08_Protein/WGCNA_J1/wgcna_j1_results.rds")

# ================================================================
# 8.10 多方法特征选择
# ================================================================

cat("\n========== 8.10 Feature Selection ==========\n")

X_tr <- t(expr_ph)
y_tr <- ifelse(grp_ph == "PD", 1, 0)

# --- LASSO ---
cat("LASSO...\n")
cv_l <- cv.glmnet(as.matrix(X_tr), y_tr, family = "binomial", alpha = 1,
                    type.measure = "auc",
                    nfolds = min(10, min(sum(y_tr == 0), sum(y_tr == 1))))
coefs_l <- coef(cv_l, s = "lambda.min")
coefs_l_trimmed <- coefs_l[-1, , drop = FALSE]
lasso_f <- rownames(coefs_l_trimmed)[which(coefs_l_trimmed != 0)]

# --- RF Importance ---
cat("RF Importance...\n")
rf_tmp <- randomForest(X_tr, as.factor(y_tr), ntree = 1000, importance = TRUE)
rf_imp <- importance(rf_tmp)
rf_df <- data.frame(Protein = rownames(rf_imp), MDA = rf_imp[, "MeanDecreaseAccuracy"],
                     MDG = rf_imp[, "MeanDecreaseGini"], stringsAsFactors = FALSE) %>%
  arrange(desc(MDG))
gini <- sort(rf_df$MDG, decreasing = TRUE)
elbow <- which(diff(diff(gini)) > 0)[1]
if (is.na(elbow) || elbow < 10) elbow <- min(50, length(gini))
rf_f <- rf_df$Protein[1:elbow]

# --- Permutation Importance ---
cat("Permutation Importance...\n")
perm_max <- numeric(100)
for (i in 1:100) {
  rf_s <- randomForest(X_tr, as.factor(sample(y_tr)), ntree = 500, importance = TRUE)
  perm_max[i] <- max(importance(rf_s)[, "MeanDecreaseAccuracy"])
}
thresh <- quantile(perm_max, 0.95)
perm_f <- rf_df %>% filter(MDA > thresh) %>% pull(Protein)

# --- Elastic Net ---
cat("Elastic Net...\n")
cv_e <- cv.glmnet(as.matrix(X_tr), y_tr, family = "binomial", alpha = 0.5,
                    type.measure = "auc")
coefs_e <- coef(cv_e, s = "lambda.min")
coefs_e_trimmed <- coefs_e[-1, , drop = FALSE]
enet_f <- rownames(coefs_e_trimmed)[which(coefs_e_trimmed != 0)]

# --- WGCNA ---
wgcna_f <- unique(unlist(lapply(sig_modules, function(m) {
  colnames(datExpr)[moduleColors == sub("^ME", "", m)]
})))

# --- 共识 ---
feat_list <- list(LASSO = lasso_f, RF = rf_f, PermImp = perm_f, ENet = enet_f, WGCNA = wgcna_f)
fc_tab <- table(unlist(feat_list)); fc_v <- as.integer(fc_tab); names(fc_v) <- names(fc_tab)
consensus_prot <- names(fc_v[fc_v >= 2])
cat("Consensus features (>=2 methods):", length(consensus_prot), "\n")

if (length(consensus_prot) < 3) {
  cat("Supplementing from top DE proteins\n")
  top_de <- de_protein %>% filter(significance == "DE") %>%
    arrange(adj.P.Val) %>% head(10) %>% pull(Protein)
  consensus_prot <- unique(c(consensus_prot, top_de))
}

pdf("results/08_Protein/ML/feature_upset.pdf", width = 10, height = 6)
UpSetR::upset(UpSetR::fromList(feat_list), sets = names(feat_list),
              order.by = "freq", mb.ratio = c(0.6, 0.4), text.scale = 1.1)
dev.off()

feat_sum <- data.frame(Protein = consensus_prot, stringsAsFactors = FALSE)
feat_sum$n_methods <- fc_v[feat_sum$Protein]
feat_sum$in_LASSO  <- feat_sum$Protein %in% lasso_f
feat_sum$in_RF     <- feat_sum$Protein %in% rf_f
feat_sum$in_PermImp<- feat_sum$Protein %in% perm_f
feat_sum$in_ENet   <- feat_sum$Protein %in% enet_f
feat_sum$in_WGCNA  <- feat_sum$Protein %in% wgcna_f
feat_sum <- feat_sum[order(feat_sum$n_methods, decreasing = TRUE), ]
write_csv(feat_sum, "results/08_Protein/ML/consensus_proteins.csv")

# ================================================================
# 8.11 分类模型：5-Fold CV + FJMU 独立验证
# ================================================================

cat("\n========== 8.11 Classification Models ==========\n")

X_tr_sel <- train_mat[consensus_prot, ]
y_tr_all <- ifelse(combined_meta$Group == "PD", 1, 0)
X_te_sel <- test_mat[consensus_prot, ]
y_te     <- ifelse(fjmu_meta$Group == "PD", 1, 0)

cat("Train:", length(y_tr_all), "(PD:", sum(y_tr_all == 1), "HC:", sum(y_tr_all == 0), ")\n")
cat("Test: ", length(y_te),     "(PD:", sum(y_te == 1),     "HC:", sum(y_te == 0),     ")\n")

# --- 5-Fold CV ---
cv_folds <- createFolds(y_tr_all, k = 5, returnTrain = TRUE)
cv_aucs <- list(LR = c(), SVM = c(), RF = c(), XGBoost = c())

for (fi in seq_along(cv_folds)) {
  tr <- cv_folds[[fi]]; vl <- setdiff(1:length(y_tr_all), tr)
  Xcv_tr <- X_tr_sel[, tr, drop = FALSE]; ycv_tr <- y_tr_all[tr]
  Xcv_vl <- X_tr_sel[, vl, drop = FALSE]; ycv_vl <- y_tr_all[vl]

  tryCatch({
    m <- cv.glmnet(t(Xcv_tr), ycv_tr, family = "binomial", alpha = 0.01,
                    nfolds = 5, type.measure = "auc")
    p <- predict(m, t(Xcv_vl), s = "lambda.min", type = "response")
    cv_aucs$LR <- c(cv_aucs$LR, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$LR <<- c(cv_aucs$LR, NA))

  tryCatch({
    m <- svm(t(Xcv_tr), as.factor(ycv_tr), kernel = "radial", probability = TRUE)
    p <- attr(predict(m, t(Xcv_vl), probability = TRUE), "probabilities")[, "1"]
    cv_aucs$SVM <- c(cv_aucs$SVM, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$SVM <<- c(cv_aucs$SVM, NA))

  tryCatch({
    m <- randomForest(t(Xcv_tr), as.factor(ycv_tr), ntree = 300)
    p <- predict(m, t(Xcv_vl), type = "prob")[, "1"]
    cv_aucs$RF <- c(cv_aucs$RF, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$RF <<- c(cv_aucs$RF, NA))

  tryCatch({
    m <- xgb.train(
      params = list(objective = "binary:logistic", eval_metric = "auc",
                    max_depth = 3, eta = 0.1),
      data = xgb.DMatrix(t(Xcv_tr), label = ycv_tr), nrounds = 50, verbose = 0)
    p <- predict(m, xgb.DMatrix(t(Xcv_vl)))
    cv_aucs$XGBoost <- c(cv_aucs$XGBoost, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$XGBoost <<- c(cv_aucs$XGBoost, NA))
}

cv_sum <- data.frame(
  Model = names(cv_aucs),
  Mean_AUC = sapply(cv_aucs, mean, na.rm = TRUE),
  SD_AUC   = sapply(cv_aucs, sd, na.rm = TRUE),
  stringsAsFactors = FALSE
)
cat("\n5-Fold CV AUC:\n"); print(cv_sum)

# --- FJMU 独立验证 ---
models <- list()

m_lr <- cv.glmnet(t(X_tr_sel), y_tr_all, family = "binomial",
                    alpha = 0.01, type.measure = "auc")
models$LR <- predict(m_lr, t(X_te_sel), s = "lambda.min", type = "response")

m_svm <- svm(t(X_tr_sel), as.factor(y_tr_all), kernel = "radial", probability = TRUE)
models$SVM <- attr(predict(m_svm, t(X_te_sel), probability = TRUE),
                    "probabilities")[, "1"]

m_rf <- randomForest(t(X_tr_sel), as.factor(y_tr_all), ntree = 500, importance = TRUE)
models$RF <- predict(m_rf, t(X_te_sel), type = "prob")[, "1"]

m_xgb <- xgb.train(
  params = list(objective = "binary:logistic", eval_metric = "auc",
                max_depth = 3, eta = 0.1),
  data = xgb.DMatrix(t(X_tr_sel), label = y_tr_all), nrounds = 100,
  verbose = 0, early_stopping_rounds = 10,
  watchlist = list(train = xgb.DMatrix(t(X_tr_sel), label = y_tr_all)))
models$XGBoost <- predict(m_xgb, xgb.DMatrix(t(X_te_sel)))

cat("\n--- Independent Test (FJMU) AUC ---\n")
for (nm in names(models)) {
  cat(sprintf("  %-8s: AUC = %.3f\n", nm,
              as.numeric(roc(y_te, models[[nm]], direction = "<", quiet = TRUE)$auc)))
}

# CV AUC 柱状图
p_cv <- ggplot(cv_sum, aes(x = reorder(Model, -Mean_AUC), y = Mean_AUC, fill = Model)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, Mean_AUC - SD_AUC),
                     ymax = pmin(1, Mean_AUC + SD_AUC)), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", Mean_AUC)), vjust = -0.8, size = 3.5) +
  scale_fill_manual(values = c("LR" = "#E74C3C", "SVM" = "#3498DB",
                                "RF" = "#27AE60", "XGBoost" = "#F39C12")) +
  labs(title = "Protein Signature: 5-Fold CV AUC", x = "", y = "AUC") +
  theme_bw(base_size = 13) + theme(legend.position = "none") + ylim(0, 1.05)
ggsave("results/08_Protein/ML/CV_AUC_comparison.pdf", p_cv, width = 7, height = 5)

# ROC 曲线
roc_lr  <- roc(y_te, models$LR,     direction = "<", quiet = TRUE)
roc_svm <- roc(y_te, models$SVM,    direction = "<", quiet = TRUE)
roc_rf  <- roc(y_te, models$RF,     direction = "<", quiet = TRUE)
roc_xgb <- roc(y_te, models$XGBoost,direction = "<", quiet = TRUE)

p_roc <- ggroc(list(LR = roc_lr, SVM = roc_svm, RF = roc_rf, XGBoost = roc_xgb),
               linewidth = 0.8) +
  geom_abline(linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("LR" = "#E74C3C", "SVM" = "#3498DB",
                                "RF" = "#27AE60", "XGBoost" = "#F39C12")) +
  annotate("text", x = 0.58, y = 0.30,
           label = sprintf("LR: AUC = %.3f", round(auc(roc_lr), 3)),
           color = "#E74C3C", size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.58, y = 0.23,
           label = sprintf("SVM: AUC = %.3f", round(auc(roc_svm), 3)),
           color = "#3498DB", size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.58, y = 0.16,
           label = sprintf("RF: AUC = %.3f", round(auc(roc_rf), 3)),
           color = "#27AE60", size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.58, y = 0.09,
           label = sprintf("XGBoost: AUC = %.3f", round(auc(roc_xgb), 3)),
           color = "#F39C12", size = 2.8, hjust = 0, fontface = "bold") +
  labs(title = "Protein ROC \u2014 Independent Test (FJMU)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave("results/08_Protein/ML/ROC_FJMU.pdf", p_roc, width = 7, height = 6)

# ================================================================
# 8.12 RF + LR Bootstrap 稳定性
# ================================================================

cat("\n========== 8.12 Bootstrap Stability (RF + LR) ==========\n")

X_full_tr <- t(train_mat)
n_boot <- 200

# RF bootstrap 存储
boot_aucs_rf <- c(); boot_sens_rf <- c(); boot_spec_rf <- c()
boot_feat_freq_rf <- rep(0, length(consensus_prot))
names(boot_feat_freq_rf) <- consensus_prot

# LR bootstrap 存储
boot_aucs_lr <- c(); boot_sens_lr <- c(); boot_spec_lr <- c()

for (b in 1:n_boot) {
  idx <- sample(1:nrow(X_full_tr), nrow(X_full_tr), replace = TRUE)
  oob <- setdiff(1:nrow(X_full_tr), unique(idx))
  if (length(oob) < 5) next
  if (length(unique(y_tr_all[idx])) < 2 || length(unique(y_tr_all[oob])) < 2) next

  X_b <- X_full_tr[idx, consensus_prot, drop = FALSE]
  X_o <- X_full_tr[oob,  consensus_prot, drop = FALSE]
  y_b <- y_tr_all[idx]
  y_o <- y_tr_all[oob]

  # --- RF ---
  tryCatch({
    rf_b <- randomForest(X_b, as.factor(y_b), ntree = 500, importance = TRUE)
    po <- predict(rf_b, X_o, type = "prob")[, "1"]
    pr <- ifelse(po > 0.5, 1, 0)
    boot_aucs_rf <- c(boot_aucs_rf, as.numeric(roc(y_o, po, direction = "<", quiet = TRUE)$auc))
    tp <- sum(pr == 1 & y_o == 1); fn <- sum(pr == 0 & y_o == 1)
    tn <- sum(pr == 0 & y_o == 0); fp <- sum(pr == 1 & y_o == 0)
    boot_sens_rf <- c(boot_sens_rf, tp / max(tp + fn, 1))
    boot_spec_rf <- c(boot_spec_rf, tn / max(tn + fp, 1))
    imp_b <- importance(rf_b)[, "MeanDecreaseGini"]
    top_f <- names(sort(imp_b, decreasing = TRUE))[1:length(consensus_prot)]
    boot_feat_freq_rf[top_f] <- boot_feat_freq_rf[top_f] + 1
  }, error = function(e) {})

  # --- LR (正则化) ---
  tryCatch({
    cv_lr <- cv.glmnet(X_b, y_b, family = "binomial", alpha = 0.01,
                        nfolds = 5, type.measure = "auc")
    po <- predict(cv_lr, X_o, s = "lambda.min", type = "response")[, 1]
    pr <- ifelse(po > 0.5, 1, 0)
    boot_aucs_lr <- c(boot_aucs_lr, as.numeric(roc(y_o, po, direction = "<", quiet = TRUE)$auc))
    tp <- sum(pr == 1 & y_o == 1); fn <- sum(pr == 0 & y_o == 1)
    tn <- sum(pr == 0 & y_o == 0); fp <- sum(pr == 1 & y_o == 0)
    boot_sens_lr <- c(boot_sens_lr, tp / max(tp + fn, 1))
    boot_spec_lr <- c(boot_spec_lr, tn / max(tn + fp, 1))
  }, error = function(e) {})
}

n_ok_rf <- length(boot_aucs_rf)
n_ok_lr <- length(boot_aucs_lr)

# --- 统计汇总 ---
cat(sprintf("\nRF Bootstrap (n=%d):  AUC=%.3f [%.3f, %.3f]  Se=%.3f  Sp=%.3f\n",
            n_ok_rf, mean(boot_aucs_rf), quantile(boot_aucs_rf, .025),
            quantile(boot_aucs_rf, .975), mean(boot_sens_rf), mean(boot_spec_rf)))
cat(sprintf("LR Bootstrap (n=%d):  AUC=%.3f [%.3f, %.3f]  Se=%.3f  Sp=%.3f\n",
            n_ok_lr, mean(boot_aucs_lr), quantile(boot_aucs_lr, .025),
            quantile(boot_aucs_lr, .975), mean(boot_sens_lr), mean(boot_spec_lr)))

# --- 图 8.12A: AUC 分布对比 ---
df_auc_compare <- bind_rows(
  data.frame(AUC = boot_aucs_rf, Model = "RF", stringsAsFactors = FALSE),
  data.frame(AUC = boot_aucs_lr, Model = "LR", stringsAsFactors = FALSE)
)

p_auc_compare <- ggplot(df_auc_compare, aes(x = AUC, fill = Model)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(data = df_auc_compare %>% group_by(Model) %>%
               summarise(m = mean(AUC), .groups = "drop"),
             aes(xintercept = m, color = Model), linewidth = 1.2) +
  scale_fill_manual(values = c("RF" = "#27AE60", "LR" = "#E74C3C")) +
  scale_color_manual(values = c("RF" = "#1E8449", "LR" = "#C0392B")) +
  annotate("text", x = mean(boot_aucs_rf), y = Inf, vjust = 2.5,
           label = sprintf("RF: %.3f [%.3f, %.3f]",
                           mean(boot_aucs_rf), quantile(boot_aucs_rf, .025),
                           quantile(boot_aucs_rf, .975)),
           color = "#1E8449", size = 3.2, hjust = -0.05) +
  annotate("text", x = mean(boot_aucs_lr), y = Inf, vjust = 4.5,
           label = sprintf("LR: %.3f [%.3f, %.3f]",
                           mean(boot_aucs_lr), quantile(boot_aucs_lr, .025),
                           quantile(boot_aucs_lr, .975)),
           color = "#C0392B", size = 3.2, hjust = -0.05) +
  labs(title = sprintf("Bootstrap AUC: RF (n=%d) vs LR (n=%d)", n_ok_rf, n_ok_lr),
       x = "AUC", y = "Count") +
  theme_bw(base_size = 13) + theme(legend.position = "none")
ggsave("results/08_Protein/ML/bootstrap_AUC_RF_vs_LR.pdf", p_auc_compare, width = 8, height = 5)

# --- 图 8.12B: 三指标汇总对比 ---
boot_sum_compare <- data.frame(
  Model  = rep(c("RF", "LR"), each = 3),
  Metric = rep(c("AUC", "Sensitivity", "Specificity"), 2),
  Mean   = c(mean(boot_aucs_rf), mean(boot_sens_rf), mean(boot_spec_rf),
             mean(boot_aucs_lr), mean(boot_sens_lr), mean(boot_spec_lr)),
  Lower  = c(quantile(boot_aucs_rf, .025), quantile(boot_sens_rf, .025),
             quantile(boot_spec_rf, .025),
             quantile(boot_aucs_lr, .025), quantile(boot_sens_lr, .025),
             quantile(boot_spec_lr, .025)),
  Upper  = c(quantile(boot_aucs_rf, .975), quantile(boot_sens_rf, .975),
             quantile(boot_spec_rf, .975),
             quantile(boot_aucs_lr, .975), quantile(boot_sens_lr, .975),
             quantile(boot_spec_lr, .975)),
  stringsAsFactors = FALSE
)
boot_sum_compare$Metric <- factor(boot_sum_compare$Metric,
                                   levels = c("AUC", "Sensitivity", "Specificity"))

p_boot_compare <- ggplot(boot_sum_compare,
                          aes(x = Metric, y = Mean, fill = Model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Mean), y = Upper),
            position = position_dodge(width = 0.7), vjust = -0.6, size = 3) +
  scale_fill_manual(values = c("RF" = "#27AE60", "LR" = "#E74C3C")) +
  labs(title = "Bootstrap Stability: RF vs LR",
       subtitle = "Bar = mean, Error bar = 95% CI",
       x = "", y = "Value") +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.85, 0.85)) +
  ylim(0, max(boot_sum_compare$Upper) * 1.18)
ggsave("results/08_Protein/ML/bootstrap_RF_vs_LR_summary.pdf", p_boot_compare, width = 8, height = 5)

# --- 图 8.12C: LR 单独三指标 ---
boot_sum_lr <- data.frame(
  Metric = factor(c("AUC", "Sensitivity", "Specificity"),
                   levels = c("AUC", "Sensitivity", "Specificity")),
  Mean  = c(mean(boot_aucs_lr), mean(boot_sens_lr), mean(boot_spec_lr)),
  Lower = c(quantile(boot_aucs_lr, .025), quantile(boot_sens_lr, .025),
            quantile(boot_spec_lr, .025)),
  Upper = c(quantile(boot_aucs_lr, .975), quantile(boot_sens_lr, .975),
            quantile(boot_spec_lr, .975))
)

p_boot_lr <- ggplot(boot_sum_lr, aes(x = Metric, y = Mean, fill = Metric)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.5) +
  geom_text(aes(y = Upper, label = sprintf("%.3f", Mean)), vjust = -0.8, size = 3.2) +
  scale_fill_manual(values = c("AUC" = "#E74C3C", "Sensitivity" = "#3498DB",
                                "Specificity" = "#F39C12")) +
  labs(title = sprintf("LR Bootstrap Stability (n=%d)", n_ok_lr),
       x = "", y = "Value") +
  theme_bw(base_size = 13) + theme(legend.position = "none") +
  ylim(0, max(boot_sum_lr$Upper) * 1.18)
ggsave("results/08_Protein/ML/bootstrap_LR_stability.pdf", p_boot_lr, width = 7, height = 5)

# --- 保存结果 ---
write_csv(boot_sum_compare, "results/08_Protein/ML/bootstrap_comparison_RF_LR.csv")

cat("\nBootstrap files:\n")
cat("  bootstrap_AUC_RF_vs_LR.pdf\n")
cat("  bootstrap_RF_vs_LR_summary.pdf\n")
cat("  bootstrap_LR_stability.pdf\n")
cat("  bootstrap_comparison_RF_LR.csv\n")

# ================================================================
# 8.13 共表达网络
# ================================================================

cat("\n========== 8.13 Co-expression Network ==========\n")

prot_cor <- cor(t(train_mat[consensus_prot, ]), method = "spearman")
prot_cor[abs(prot_cor) < 0.4] <- 0
edges <- which(prot_cor != 0 & upper.tri(prot_cor), arr.ind = TRUE)

if (nrow(edges) > 0) {
  g <- graph_from_data_frame(
    data.frame(from = rownames(prot_cor)[edges[, 1]],
               to   = colnames(prot_cor)[edges[, 2]],
               weight = abs(prot_cor[edges])), directed = FALSE)
  V(g)$is_DE <- V(g)$name %in% (de_protein %>% filter(significance == "DE") %>% pull(Protein))
  V(g)$degree <- degree(g)
  write_graph(g, "results/08_Protein/ML/PPI_network.graphml", format = "graphml")
  hub_df <- data.frame(Protein = V(g)$name, Degree = V(g)$degree) %>% arrange(desc(Degree))
  write_csv(hub_df, "results/08_Protein/ML/network_hubs.csv")
  cat("Network:", vcount(g), "nodes,", ecount(g), "edges\n")
}

# ================================================================
# 8.14 保存富集分析输入 + 整合输入
# ================================================================

cat("\n========== 8.14 Saving Enrichment Inputs ==========\n")

write_csv(data.frame(Protein = consensus_prot),
          "results/08_Protein/Enrichment_input_consensus_proteins.csv")
write_csv(de_protein %>% filter(significance == "DE") %>%
            dplyr::select(Protein, logFC, adj.P.Val, direction),
          "results/08_Protein/Enrichment_input_DE_proteins.csv")
write_csv(hub_summary, "results/08_Protein/Enrichment_input_WGCNA_hubs.csv")
write_csv(data.frame(Protein = colnames(datExpr), Module = moduleColors,
                      stringsAsFactors = FALSE),
          "results/08_Protein/Enrichment_input_WGCNA_all_genes.csv")
write_csv(data.frame(Protein = unique(c(consensus_prot, de_protein$Protein))),
          "results/08_Protein/Enrichment_input_all_proteins_for_mapping.csv")

ml_res <- readRDS("data/processed/ml_results.rds")
write_csv(data.frame(Feature = consensus_prot, Type = "Protein"),
          "results/08_Protein/Integration/protein_signature_for_integration.csv")
write_csv(data.frame(Feature = ml_res$consensus_features, Type = "miRNA"),
          "results/08_Protein/Integration/miRNA_signature_for_integration.csv")

cat("Enrichment input files saved. Run 08b_Protein_Enrichment.R separately.\n")

# ================================================================
# 8.15 保存所有结果
# ================================================================

protein_results <- list(
  train_mat = train_mat, test_mat = test_mat,
  combined_meta = combined_meta, fjmu_meta = fjmu_meta,
  shared_proteins = shared_all,
  # 各数据集独立 DEA
  de_j1 = de_j1, de_j2 = de_j2,
  # 合并 DEA
  de_protein = de_protein,
  # WGCNA combined
  wgcna = wgcna_results,
  # WGCNA J1
  wgcna_j1 = wgcna_j1_results,
  # ML
  consensus_proteins = consensus_prot,
  feat_summary = feat_sum, feat_list = feat_list,
  cv_summary = cv_sum, models = models,
  bootstrap = list(
				   rf = list(aucs = boot_aucs_rf, sens = boot_sens_rf,
							 spec = boot_spec_rf, feat_freq = boot_feat_freq_rf),
				   lr  = list(aucs = boot_aucs_lr, sens = boot_sens_lr, spec = boot_spec_lr)
				   ),
  rf_model = m_rf
)
saveRDS(protein_results, "data/processed/protein_results.rds")

cat("\n========================================\n")
cat("Protein Analysis Complete\n")
cat("========================================\n")
cat(sprintf("  Shared proteins (all 3): %d\n", length(shared_all)))
cat(sprintf("  J1 DE proteins: %d\n", sum(de_j1$significance == "DE")))
cat(sprintf("  J2 DE proteins: %d\n", sum(de_j2$significance == "DE")))
cat(sprintf("  Combined DE proteins: %d\n", sum(de_protein$significance == "DE")))
cat(sprintf("  Shared Up   (J1∩J2∩Comb): %d\n", length(shared_up)))
cat(sprintf("  Shared Down (J1∩J2∩Comb): %d\n", length(shared_down)))
cat(sprintf("  WGCNA combined modules: %d, significant: %d\n",
            length(unique(moduleColors)), length(sig_modules)))
cat(sprintf("  WGCNA J1 modules: %d, significant: %d\n",
            length(unique(moduleColors_j1)), length(sig_mods_j1)))
cat(sprintf("  Consensus signature: %d\n", length(consensus_prot)))
cat(sprintf("  Best CV AUC: %s = %.3f\n",
            cv_sum$Model[which.max(cv_sum$Mean_AUC)], max(cv_sum$Mean_AUC)))
cat(sprintf("  FJMU best AUC: %s = %.3f\n",
            names(models)[which.max(sapply(models, function(m)
              as.numeric(roc(y_te, m, direction="<", quiet=TRUE)$auc)))],
            max(sapply(models, function(m)
              as.numeric(roc(y_te, m, direction="<", quiet=TRUE)$auc)))))
cat(sprintf("  Bootstrap AUC: %.3f [%.3f, %.3f]\n",
            mean(boot_aucs_lr), quantile(boot_aucs_lr,.025), quantile(boot_aucs_lr,.975)))
cat("\nRun 08b_Protein_Enrichment.R for enrichment & integration.\n")
