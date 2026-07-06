# ============================================================
# Figure 2: Protein Analysis Publication Figure
# npj Parkinson's Disease style
# ============================================================

library(tidyverse)
library(pROC)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(sva)
library(limma)

prot_res <- readRDS("data/processed/protein_results.rds")

dir.create("figures/publication", recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Constants & Theme
# ================================================================

col_pd   <- "#E74C3C";  col_hc   <- "#3498DB";  col_msa  <- "#27AE60"
col_g1   <- "#E67E22";  col_g2   <- "#27AE60";  col_fjmu <- "#8E44AD"
col_lr   <- "#E74C3C";  col_svm  <- "#3498DB"
col_rf   <- "#27AE60";  col_xgb  <- "#F39C12"

model_colors <- c("LR" = col_lr, "SVM" = col_svm,
                   "RF" = col_rf, "XGBoost" = col_xgb)

theme_j <- theme_bw(base_size = 9) +
  theme(
    text               = element_text(family = "sans"),
    plot.title         = element_text(size = 10, hjust = 0.5,
                                       margin = ggplot2::margin(b = 5)),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 7.5, color = "black"),
    legend.key.size    = unit(0.35, "cm"),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 8, face = "bold"),
    legend.margin      = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin  = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    plot.margin        = ggplot2::margin(3, 5, 3, 3, unit = "pt")
  )

# ================================================================
# Load raw data for PCA panels
# ================================================================

# 各数据集预处理后数据（从 protein_results 中获取）
train_mat <- prot_res$train_mat
test_mat  <- prot_res$test_mat
combined_meta <- prot_res$combined_meta
fjmu_meta     <- prot_res$fjmu_meta

# 读取原始数据用于各数据集独立 PCA
j1_expr <- as.matrix(read.csv("data/protein/JPST003013_protein.csv",
                               row.names = 1, check.names = FALSE))
j2_expr <- as.matrix(read.csv("data/protein/JPST003026_protein.csv",
                               row.names = 1, check.names = FALSE))
fjmu_expr <- as.matrix(read.csv("data/protein/FJMU_protein.csv",
                                  row.names = 1, check.names = FALSE))
j1_meta <- read.csv("data/protein/meta_JPST003013.csv", row.names = 1,
                      check.names = FALSE)
j2_meta <- read.csv("data/protein/meta_JPST003026.csv", row.names = 1,
                      check.names = FALSE)
fjmu_meta_raw <- read.csv("data/protein/FJMU_protein_meta.csv", row.names = 1,
                            check.names = FALSE)

clean_group <- function(g) {
  g <- trimws(as.character(g))
  ifelse(g %in% c("PDND", "PD-MCI", "PDD"), "PD",
  ifelse(g %in% c("MSA-C", "MSA-P"), "MSA",
  ifelse(g == "HC", "HC", g)))
}

j1_meta$Group     <- factor(clean_group(j1_meta$Group),     levels = c("HC", "PD", "MSA"))
j2_meta$Group     <- factor(clean_group(j2_meta$Group),     levels = c("HC", "PD", "MSA"))
fjmu_meta_raw$Group <- factor(clean_group(fjmu_meta_raw$Group), levels = c("HC", "PD", "MSA"))

# ================================================================
# Preprocessing helper
# ================================================================

preprocess_expr <- function(mat, min_detect = 0.5) {
  detected <- rowMeans(!is.na(mat) & mat > 0)
  mat <- mat[detected >= min_detect, ]
  min_val <- min(mat[mat > 0], na.rm = TRUE) / 2
  mat[is.na(mat) | mat <= 0] <- min_val
  mat <- log2(mat)
  return(mat)
}

j1_log2   <- preprocess_expr(j1_expr)
j2_log2   <- preprocess_expr(j2_expr)
fjmu_log2 <- preprocess_expr(fjmu_expr)

# 各数据集独立归一化（用于 Row 1 独立 PCA）
j1_norm   <- normalizeBetweenArrays(j1_log2, method = "quantile")
j2_norm   <- normalizeBetweenArrays(j2_log2, method = "quantile")
fjmu_norm <- normalizeBetweenArrays(fjmu_log2, method = "quantile")

# 三数据集共有蛋白
shared_all <- Reduce(intersect, list(rownames(j1_norm), rownames(j2_norm),
                                      rownames(fjmu_norm)))

# ================================================================
# Panel a–c: PCA per dataset (raw, independent)
# ================================================================

make_pca_prot <- function(mat, meta, ttl, show_legend = FALSE) {
  # 只保留 PD + HC
  keep <- meta$Group %in% c("PD", "HC")
  mat_sub <- mat[, keep]
  grp_sub <- droplevels(meta$Group[keep])

  pca <- prcomp(t(mat_sub), scale. = TRUE)
  ve  <- summary(pca)$importance[2, 1:2] * 100
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = grp_sub)

  p <- ggplot(df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.85) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.6) +
    scale_color_manual(values = c("PD" = col_pd, "HC" = col_hc)) +
    labs(title = ttl,
         x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_j

  if (show_legend) {
    p <- p + theme(legend.position = c(0.95, 0.98),
                    legend.justification = c(1, 1),
                    legend.background = element_rect(fill = "white", color = NA,
                                                      linewidth = 0.3),
                    legend.key = element_blank())
  } else {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

p_a <- make_pca_prot(j1_norm, j1_meta, "JPST003013")
p_b <- make_pca_prot(j2_norm, j2_meta, "JPST003026")
p_c <- make_pca_prot(fjmu_norm, fjmu_meta_raw, "FJMU", show_legend = TRUE)

leg_group <- get_legend(p_c)
p_c <- p_c + theme(legend.position = "none")

row1 <- plot_grid(
  plot_grid(p_a, p_b, p_c, ncol = 3, rel_widths = c(1, 1, 1),
            labels = c("a", "b", "c"), label_size = 12,
            label_fontface = "bold", label_x = 0, label_y = 1),
  leg_group,
  ncol = 1, rel_heights = c(1, 0.06)
)

# ================================================================
# Panel d: J1+J2 ComBat (training set)
# ================================================================

j1_shared <- j1_norm[shared_all, ]
j2_shared <- j2_norm[shared_all, ]

combined_raw_de <- cbind(j1_shared, j2_shared)
batch_de <- c(rep("JPST003013", ncol(j1_shared)),
              rep("JPST003026", ncol(j2_shared)))

cat("Running ComBat for Figure 2d (J1+J2)...\n")
combined_combat_de <- ComBat(dat = combined_raw_de, batch = batch_de,
                               mod = NULL, par.prior = TRUE)
combined_norm_de <- normalizeBetweenArrays(combined_combat_de, method = "quantile")

meta_de <- data.frame(
  Dataset = batch_de,
  Group   = c(j1_meta$Group, j2_meta$Group),
  stringsAsFactors = FALSE
)

# 仅 PD+HC
keep_de <- meta_de$Group %in% c("PD", "HC")
pca_de  <- prcomp(t(combined_norm_de[, keep_de]), scale. = TRUE)
ve_de   <- summary(pca_de)$importance[2, 1:2] * 100
df_de   <- data.frame(PC1 = pca_de$x[, 1], PC2 = pca_de$x[, 2],
                       Dataset = meta_de$Dataset[keep_de],
                       Group   = droplevels(meta_de$Group[keep_de]))

p_d <- ggplot(df_de, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8, alpha = 0.85) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("JPST003013" = col_g1, "JPST003026" = col_g2)) +
  scale_shape_manual(values = c("PD" = 16, "HC" = 17)) +
  labs(title = "Training Set (J1+J2 ComBat)",
       x = sprintf("PC1 (%.1f%%)", ve_de[1]),
       y = sprintf("PC2 (%.1f%%)", ve_de[2])) +
  theme_j +
  theme(legend.position = "none")

# ================================================================
# Panel e: All 3 datasets ComBat
# ================================================================

fjmu_shared <- fjmu_norm[shared_all, ]
combined_all_raw <- cbind(j1_shared, j2_shared, fjmu_shared)
batch_all <- c(rep("JPST003013", ncol(j1_shared)),
               rep("JPST003026", ncol(j2_shared)),
               rep("FJMU",       ncol(fjmu_shared)))

cat("Running ComBat for Figure 2e (all 3 datasets)...\n")
combined_all_combat <- ComBat(dat = combined_all_raw, batch = batch_all,
                                mod = NULL, par.prior = TRUE)
combined_all_norm <- normalizeBetweenArrays(combined_all_combat, method = "quantile")

meta_all <- data.frame(
  Dataset = batch_all,
  Group   = c(j1_meta$Group, j2_meta$Group, fjmu_meta_raw$Group),
  stringsAsFactors = FALSE
)

keep_all <- meta_all$Group %in% c("PD", "HC")
pca_all  <- prcomp(t(combined_all_norm[, keep_all]), scale. = TRUE)
ve_all   <- summary(pca_all)$importance[2, 1:2] * 100
df_all   <- data.frame(PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2],
                        Dataset = meta_all$Dataset[keep_all],
                        Group   = droplevels(meta_all$Group[keep_all]))

p_e <- ggplot(df_all, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8, alpha = 0.85) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("JPST003013" = col_g1,
                                 "JPST003026" = col_g2,
                                 "FJMU"       = col_fjmu)) +
  scale_shape_manual(values = c("PD" = 16, "HC" = 17)) +
  labs(title = "All Datasets (ComBat)",
       x = sprintf("PC1 (%.1f%%)", ve_all[1]),
       y = sprintf("PC2 (%.1f%%)", ve_all[2])) +
  theme_j +
  theme(legend.position = "none")

# 共享 legend (Dataset + Group)
p_legend_src <- ggplot(df_all, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("JPST003013" = col_g1,
                                 "JPST003026" = col_g2,
                                 "FJMU"       = col_fjmu)) +
  scale_shape_manual(values = c("PD" = 16, "HC" = 17)) +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 2),
         linetype = guide_legend(order = 3)) +
  theme_j +
  theme(
    legend.position      = "bottom",
    legend.justification = "left",
    legend.box.just      = "left",
    legend.box           = "horizontal",
    legend.box.margin    = ggplot2::margin(0, 0, 0, 0, unit = "pt")
  )
leg_dataset <- get_legend(p_legend_src)

# ================================================================
# Panel f: Volcano Plot (PD vs HC, combined)
# ================================================================

de_protein <- prot_res$de_protein

# 自动检测列名
pcol <- "Protein"
if (!"Protein" %in% colnames(de_protein)) {
  pcol <- colnames(de_protein)[1]
}

top_label <- de_protein %>%
  dplyr::filter(significance == "DE") %>%
  head(15)

p_f <- ggplot(de_protein, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_text_repel(data = top_label,
                   aes(label = !!sym(pcol)),
                   size = 2.2, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("Up" = col_pd, "Down" = col_hc, "NS" = "#BDC3C7")) +
  labs(title = "DE Proteins: PD vs HC",
       x = "log2(FC)", y = "-log10(FDR)") +
  theme_j +
  theme(legend.position = "none")

# Row 2 组装
de_panel <- plot_grid(p_d, p_e, ncol = 2, rel_widths = c(1, 1),
                       labels = c("d", "e"), label_size = 12,
                       label_fontface = "bold", label_x = 0, label_y = 1)

def_panel <- plot_grid(de_panel, p_f,
                        ncol = 2, rel_widths = c(2.2, 1),
                        labels = c("", "f"), label_size = 12,
                        label_fontface = "bold", label_x = 0, label_y = 1)

row2 <- plot_grid(def_panel, leg_dataset,
                   ncol = 1, rel_heights = c(1, 0.08))

# ================================================================
# Panel g: 5-Fold CV AUC
# ================================================================

cv_sum <- prot_res$cv_summary

cv_df <- cv_sum %>%
  mutate(Model = factor(Model, levels = Model[order(-Mean_AUC)]))

p_g <- ggplot(cv_df, aes(x = Model, y = Mean_AUC, fill = Model)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, Mean_AUC - SD_AUC),
                     ymax = pmin(1, Mean_AUC + SD_AUC)),
                width = 0.2, linewidth = 0.4) +
  geom_text(aes(y = pmin(1, Mean_AUC + SD_AUC),
                 label = sprintf("%.3f", Mean_AUC)),
            vjust = -0.8, size = 2.8) +
  scale_fill_manual(values = model_colors) +
  labs(title = "5-Fold Cross-Validation AUC",
       x = "", y = "AUC") +
  theme_j +
  theme(legend.position = "none") +
  ylim(0, 1.1)

# ================================================================
# Panel h: RF Bootstrap Stability
# ================================================================

boot_data <- prot_res$bootstrap

boot_aucs_rf <- boot_data$rf$aucs
boot_sens_rf <- boot_data$rf$sens
boot_spec_rf <- boot_data$rf$spec
n_ok_rf <- length(boot_aucs_rf)

boot_sum_rf <- data.frame(
  Metric = factor(c("AUC", "Sensitivity", "Specificity"),
                   levels = c("AUC", "Sensitivity", "Specificity")),
  Mean  = c(mean(boot_aucs_rf), mean(boot_sens_rf), mean(boot_spec_rf)),
  Lower = c(quantile(boot_aucs_rf, .025), quantile(boot_sens_rf, .025),
            quantile(boot_spec_rf, .025)),
  Upper = c(quantile(boot_aucs_rf, .975), quantile(boot_sens_rf, .975),
            quantile(boot_spec_rf, .975))
)

p_h <- ggplot(boot_sum_rf, aes(x = Metric, y = Mean, fill = Metric)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.5) +
  geom_text(aes(y = Upper, label = sprintf("%.3f", Mean)),
            vjust = -0.8, size = 2.8) +
  scale_fill_manual(values = c("AUC" = col_rf, "Sensitivity" = col_pd,
                                "Specificity" = col_hc)) +
  labs(title = sprintf("RF Bootstrap Stability (n=%d)", n_ok_rf),
       x = "", y = "Value") +
  theme_j +
  theme(legend.position = "none") +
  ylim(0, max(boot_sum_rf$Upper) * 1.18)


# ================================================================
# Panel i: ROC – Independent Test (FJMU)
# ================================================================

y_te <- ifelse(fjmu_meta$Group == "PD", 1, 0)

models <- prot_res$models

roc_lr  <- roc(y_te, models$LR,      direction = "<", quiet = TRUE)
roc_svm <- roc(y_te, models$SVM,     direction = "<", quiet = TRUE)
roc_rf  <- roc(y_te, models$RF,      direction = "<", quiet = TRUE)
roc_xgb <- roc(y_te, models$XGBoost, direction = "<", quiet = TRUE)

auc_vals <- c(LR     = round(as.numeric(auc(roc_lr)),  3),
              SVM    = round(as.numeric(auc(roc_svm)), 3),
              RF     = round(as.numeric(auc(roc_rf)),  3),
              XGBoost= round(as.numeric(auc(roc_xgb)), 3))

p_i <- ggroc(list(LR = roc_lr, SVM = roc_svm, RF = roc_rf, XGBoost = roc_xgb),
             linewidth = 0.8) +
  #geom_abline(linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("LR" = col_lr, "SVM" = col_svm,
                                "RF" = col_rf, "XGBoost" = col_xgb)) +
  annotate("text", x = 0, y = 0.24,
           label = sprintf("LR: AUC = %.3f",       auc_vals["LR"]),
           color = col_lr,  size = 2.8, hjust = 1, fontface = "bold") +
  annotate("text", x = 0, y = 0.15,
           label = sprintf("RF: AUC = %.3f",       auc_vals["RF"]),
           color = col_rf,  size = 2.8, hjust = 1, fontface = "bold") +
  annotate("text", x = 0, y = 0.08,
           label = sprintf("SVM: AUC = %.3f",      auc_vals["SVM"]),
           color = col_svm, size = 2.8, hjust = 1, fontface = "bold") +
  annotate("text", x = 0, y = 0.01,
           label = sprintf("XGBoost: AUC = %.3f",  auc_vals["XGBoost"]),
           color = col_xgb, size = 2.8, hjust = 1, fontface = "bold") +
  labs(title = "ROC - Independent Test (FJMU)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_j +
  theme(legend.position = "none")

# ================================================================
# Row 3 组装
# ================================================================

row3 <- plot_grid(p_g, p_h, p_i, ncol = 3, rel_widths = c(1, 1, 1.2),
                  labels = c("g", "h", "i"), label_size = 12,
                  label_fontface = "bold", label_x = 0, label_y = 1)

# ================================================================
# Assemble Figure 2
# ================================================================

fig2 <- plot_grid(
  row1, row2, row3,
  ncol   = 1,
  rel_heights = c(1, 1, 1)
)

ggsave("figures/publication/Figure2_Protein_Analysis.pdf",
       fig2, width = 7.2, height = 8.25, units = "in", dpi = 300)

cat("\nFigure 2 saved to figures/publication/\n")
