# ============================================================
# Figure 1: Multi-panel publication figure
# npj Parkinson's Disease style
# ============================================================

library(tidyverse)
library(pROC)
library(randomForest)
library(cowplot)
library(ggplot2)

proc   <- readRDS("data/processed/processed_data.rds")
ml_res <- readRDS("data/processed/ml_results.rds")

dir.create("figures/publication", recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Helper functions and constants
# ================================================================

extract_mature <- function(x) {
  sapply(strsplit(as.character(x), "_"), function(p) {
    if (length(p) >= 2) tail(p, 1) else p[1]
  })
}

# Colors
col_pd  <- "#E74C3C";  col_hc  <- "#3498DB"
col_g1  <- "#E67E22";  col_g2  <- "#27AE60";  col_s1  <- "#8E44AD"
col_lr  <- "#E74C3C";  col_svm <- "#3498DB"
col_rf  <- "#27AE60";  col_xgb <- "#F39C12"

model_colors <- c("LR" = col_lr, "SVM" = col_svm,
                   "RF" = col_rf, "XGBoost" = col_xgb)

# Journal theme
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

# Reconstruct y_test
y_test <- ifelse(proc$s1_meta$Group == "PD", 1, 0)

# ================================================================
# Helper: PCA plot function
# ================================================================

make_pca <- function(cpm_mat, meta, ttl) {
  pca <- prcomp(t(cpm_mat), scale. = TRUE)
  ve  <- summary(pca)$importance[2, 1:2] * 100
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = meta$Group)
  ggplot(df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 1.8, alpha = 0.85) +
    stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.6) +
    scale_color_manual(values = c("PD" = col_pd, "HC" = col_hc)) +
    labs(title = ttl,
         x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_j
}

# ================================================================
# Panel a–c: PCA per dataset
# ================================================================

p_a <- make_pca(proc$g1_cpm, proc$g1_meta, "GSE269775") +
  theme(legend.position = "none")
p_b <- make_pca(proc$g2_cpm, proc$g2_meta, "GSE269776") +
  theme(legend.position = "none")

p_c <- make_pca(proc$s1_cpm, proc$s1_meta, "FJMU") +
  theme(legend.position = c(0.95, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA,
                                          linewidth = 0.3),
        legend.key = element_blank())

leg_group <- get_legend(p_c)

row1 <- plot_grid(
  plot_grid(p_a, p_b, p_c, ncol = 3, rel_widths = c(1, 1, 1),
            labels = c("a", "b", "c"), label_size = 12,
            label_fontface = "bold", label_x = 0, label_y = 1),
  leg_group,
  ncol = 1, rel_heights = c(1, 0.06)
)


# ================================================================
# Panel d: Before batch correction (all 3 datasets, raw)
# ================================================================

common_raw <- Reduce(intersect,
                      list(rownames(proc$g1_cpm), rownames(proc$g2_cpm),
                           rownames(proc$s1_cpm)))
combined_raw <- cbind(proc$g1_cpm[common_raw, ],
                       proc$g2_cpm[common_raw, ],
                       proc$s1_cpm[common_raw, ])
meta_raw <- data.frame(
  Dataset = c(rep("GSE269775", ncol(proc$g1_cpm)),
              rep("GSE269776", ncol(proc$g2_cpm)),
              rep("FJMU",      ncol(proc$s1_cpm))),
  Group   = c(proc$g1_meta$Group, proc$g2_meta$Group, proc$s1_meta$Group),
  stringsAsFactors = FALSE
)

pca_raw <- prcomp(t(combined_raw), scale. = TRUE)
ve_raw  <- summary(pca_raw)$importance[2, 1:2] * 100
df_raw  <- data.frame(PC1 = pca_raw$x[, 1], PC2 = pca_raw$x[, 2],
                       Dataset = meta_raw$Dataset, Group = meta_raw$Group)

p_d <- ggplot(df_raw, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8, alpha = 0.85) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("GSE269775" = col_g1,
                                 "GSE269776" = col_g2,
                                 "FJMU"      = col_s1)) +
  scale_shape_manual(values = c("PD" = 16, "HC" = 17)) +
  labs(title = "Before Batch Correction",
       x = sprintf("PC1 (%.1f%%)", ve_raw[1]),
       y = sprintf("PC2 (%.1f%%)", ve_raw[2])) +
  theme_j +
  theme(legend.position = "none")


# ================================================================
# Panel e: After Batch Correction (3 datasets ComBat)
# ================================================================

# 三数据集共有 miRNA，原始合并后做 ComBat
common_raw_e <- Reduce(intersect,
                         list(rownames(proc$g1_cpm), rownames(proc$g2_cpm),
                              rownames(proc$s1_cpm)))
combined_raw_e <- cbind(proc$g1_cpm[common_raw_e, ],
                          proc$g2_cpm[common_raw_e, ],
                          proc$s1_cpm[common_raw_e, ])

batch_e <- c(rep("GSE269775", ncol(proc$g1_cpm)),
             rep("GSE269776", ncol(proc$g2_cpm)),
             rep("FJMU",      ncol(proc$s1_cpm)))

cat("Running ComBat for Figure 1e (3 datasets)...\n")
combined_combat_e <- sva::ComBat(dat = combined_raw_e, batch = batch_e,
                                   mod = NULL, par.prior = TRUE)

# quantile 归一化
combined_norm_e <- limma::normalizeBetweenArrays(combined_combat_e, method = "quantile")

# PCA
meta_e <- data.frame(
  Dataset = batch_e,
  Group   = c(proc$g1_meta$Group, proc$g2_meta$Group, proc$s1_meta$Group),
  stringsAsFactors = FALSE
)

pca_combat <- prcomp(t(combined_norm_e), scale. = TRUE)
ve_combat  <- summary(pca_combat)$importance[2, 1:2] * 100
df_combat  <- data.frame(PC1 = pca_combat$x[, 1], PC2 = pca_combat$x[, 2],
                          Dataset = meta_e$Dataset,
                          Group   = meta_e$Group)

p_e <- ggplot(df_combat, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8, alpha = 0.85) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("GSE269775" = col_g1,
                                 "GSE269776" = col_g2,
                                 "FJMU"      = col_s1)) +
  scale_shape_manual(values = c("PD" = 16, "HC" = 17)) +
  labs(title = "After Batch Correction",
       x = sprintf("PC1 (%.1f%%)", ve_combat[1]),
       y = sprintf("PC2 (%.1f%%)", ve_combat[2])) +
  theme_j +
  theme(legend.position = "none")


# ================================================================
# 共享 legend（左对齐，横跨 d+e+f 三个 panel）
# ================================================================

p_legend_source <- ggplot(df_raw, aes(PC1, PC2, color = Dataset, shape = Group)) +
  geom_point(size = 1.8) +
  stat_ellipse(aes(linetype = Group), level = 0.95, linewidth = 0.6) +
  scale_color_manual(values = c("GSE269775" = col_g1,
                                 "GSE269776" = col_g2,
                                 "FJMU"      = col_s1)) +
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

leg_batch <- get_legend(p_legend_source)

# ================================================================
# Panel f: Variance Decomposition
# ================================================================

combined_cpm <- proc$g1g2_cpm_combat
meta_comb    <- proc$g1g2_meta

common_ac <- intersect(rownames(combined_cpm), rownames(proc$s1_cpm))
combined_ac <- cbind(combined_cpm[common_ac, ], proc$s1_cpm[common_ac, ])
meta_ac <- data.frame(
  Group   = c(meta_comb$Group, proc$s1_meta$Group),
  Dataset = c(meta_comb$Dataset, proc$s1_meta$Dataset),
  stringsAsFactors = FALSE
)

top_var <- head(order(apply(combined_ac, 1, var), decreasing = TRUE), 500)

var_decomp <- apply(combined_ac[top_var, ], 1, function(y) {
  a  <- anova(lm(y ~ Group + Dataset, data = meta_ac))
  ss <- sum(a[, "Sum Sq"])
  data.frame(Group    = a["Group", "Sum Sq"] / ss,
             Batch    = a["Dataset", "Sum Sq"] / ss,
             Residual = a["Residuals", "Sum Sq"] / ss)
})
var_df <- do.call(rbind, var_decomp)

var_summary <- data.frame(
  Source     = factor(c("Group", "Batch", "Residual"),
                       levels = c("Group", "Batch", "Residual")),
  Proportion = c(mean(var_df$Group), mean(var_df$Batch), mean(var_df$Residual))
)

p_f <- ggplot(var_summary, aes(x = Source, y = Proportion, fill = Source)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", Proportion * 100)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Group" = col_pd, "Batch" = col_xgb,
                                "Residual" = "#95A5A6")) +
  labs(title = "Variance Decomposition", x = "", y = "Proportion") +
  theme_j +
  theme(legend.position = "none") +
  ylim(0, max(var_summary$Proportion) * 1.3)


# ================================================================
# Row 2 组装：d + e + f，legend 横跨三个 panel
# ================================================================

# d 和 e 并排
de_panel <- plot_grid(p_d, p_e, ncol = 2, rel_widths = c(1, 1),
                       labels = c("d", "e"), label_size = 12,
                       label_fontface = "bold", label_x = 0, label_y = 1)

# d+e 并排后与 f 并排
def_panel <- plot_grid(de_panel, p_f,
                        ncol = 2, rel_widths = c(2.2, 1),
                        labels = c("", "f"), label_size = 12,
                        label_fontface = "bold", label_x = 0, label_y = 1)


# ================================================================
# Panel g: 5-Fold CV AUC
# ================================================================

cv_df <- ml_res$cv_summary %>%
  mutate(Model = factor(Model, levels = Model[order(-Mean_AUC)]))

p_g <- ggplot(cv_df, aes(x = Model, y = Mean_AUC, fill = Model)) +
  geom_col(width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, Mean_AUC - SD_AUC),
                     ymax = pmin(1, Mean_AUC + SD_AUC)),
                width = 0.2, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", Mean_AUC)), vjust = -1, size = 2.8) +
  scale_fill_manual(values = model_colors) +
  labs(title = "5-Fold Cross-Validation AUC",
       x = "", y = "AUC") +
  theme_j +
  theme(legend.position = "none") +
  ylim(0, 1.08)

# ================================================================
# Panel h: ROC – Independent Test (FJMU)
# ================================================================

roc_lr  <- roc(y_test, ml_res$models$LR$prob,     direction = "<", quiet = TRUE)
roc_svm <- roc(y_test, ml_res$models$SVM$prob,     direction = "<", quiet = TRUE)
roc_rf  <- roc(y_test, ml_res$models$RF$prob,      direction = "<", quiet = TRUE)
roc_xgb <- roc(y_test, ml_res$models$XGBoost$prob, direction = "<", quiet = TRUE)

auc_vals <- c(LR     = round(as.numeric(auc(roc_lr)),  3),
              SVM    = round(as.numeric(auc(roc_svm)), 3),
              RF     = round(as.numeric(auc(roc_rf)),  3),
              XGBoost= round(as.numeric(auc(roc_xgb)), 3))

p_h <- ggroc(list(LR = roc_lr, SVM = roc_svm, RF = roc_rf, XGBoost = roc_xgb),
             linewidth = 0.8) +
  #geom_abline(linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c("LR" = col_lr, "SVM" = col_svm,
                                "RF" = col_rf, "XGBoost" = col_xgb)) +
  annotate("text", x = 0.38, y = 0.30,
           label = sprintf("RF: AUC = %.3f",      auc_vals["RF"]),
           color = col_rf,  size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.38, y = 0.23,
           label = sprintf("XGBoost: AUC = %.3f", auc_vals["XGBoost"]),
           color = col_xgb, size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.38, y = 0.16,
           label = sprintf("SVM: AUC = %.3f",     auc_vals["SVM"]),
           color = col_svm, size = 2.8, hjust = 0, fontface = "bold") +
  annotate("text", x = 0.38, y = 0.09,
           label = sprintf("LR: AUC = %.3f",      auc_vals["LR"]),
           color = col_lr,  size = 2.8, hjust = 0, fontface = "bold") +
  labs(title = "ROC Curves - Independent Test (FJMU)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_j +
  theme(legend.position = "none")

# ================================================================
# Panel i: RF Bootstrap Summary
# ================================================================

if (!is.null(ml_res$rf_bootstrap)) {
  bd <- ml_res$rf_bootstrap
  boot_aucs <- bd$aucs;  boot_sens <- bd$sensitivities
  boot_spec <- bd$specificities;  n_ok <- bd$n_success
} else {
  cat("Computing RF bootstrap (may take several minutes)...\n")
  set.seed(42)
  X_all <- t(proc$g1g2_cpm_combat)
  s1c   <- intersect(colnames(X_all), rownames(proc$s1_cpm))
  X_all <- X_all[, s1c]
  y_all <- ifelse(proc$g1g2_meta$Group == "PD", 1, 0)
  cf    <- ml_res$consensus_features
  boot_aucs <- c(); boot_sens <- c(); boot_spec <- c()

  for (b in 1:200) {
    idx <- sample(1:nrow(X_all), nrow(X_all), replace = TRUE)
    oob <- setdiff(1:nrow(X_all), unique(idx))
    if (length(oob) < 5) next
    if (length(unique(y_all[idx])) < 2 || length(unique(y_all[oob])) < 2) next
    tryCatch({
      rf_b <- randomForest(X_all[idx, cf, drop = FALSE],
                           as.factor(y_all[idx]), ntree = 500)
      po   <- predict(rf_b, X_all[oob, cf, drop = FALSE], type = "prob")[, "1"]
      pr   <- ifelse(po > 0.5, 1, 0)
      yo   <- y_all[oob]
      boot_aucs <- c(boot_aucs,
                      as.numeric(roc(yo, po, direction = "<", quiet = TRUE)$auc))
      tp <- sum(pr == 1 & yo == 1); fn <- sum(pr == 0 & yo == 1)
      tn <- sum(pr == 0 & yo == 0); fp <- sum(pr == 1 & yo == 0)
      boot_sens <- c(boot_sens, tp / max(tp + fn, 1))
      boot_spec <- c(boot_spec, tn / max(tn + fp, 1))
    }, error = function(e) {})
  }
  n_ok <- length(boot_aucs)
}

boot_sum <- data.frame(
  Metric = factor(c("AUC", "Sensitivity", "Specificity"),
                   levels = c("AUC", "Sensitivity", "Specificity")),
  Mean  = c(mean(boot_aucs), mean(boot_sens), mean(boot_spec)),
  Lower = c(quantile(boot_aucs, .025), quantile(boot_sens, .025),
            quantile(boot_spec, .025)),
  Upper = c(quantile(boot_aucs, .975), quantile(boot_sens, .975),
            quantile(boot_spec, .975))
)

p_i <- ggplot(boot_sum, aes(x = Metric, y = Mean, fill = Metric)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, linewidth = 0.5) +
  geom_text(aes(y = Upper, label = sprintf("%.3f", Mean)),
            vjust = -0.8, size = 2.8) +
  scale_fill_manual(values = c("AUC" = col_rf, "Sensitivity" = col_pd,
                                "Specificity" = col_hc)) +
  labs(title = sprintf("RF Bootstrap Stability (n=%d)", n_ok),
       x = "", y = "Value") +
  theme_j +
  theme(legend.position = "none") +
  ylim(0, max(boot_sum$Upper) * 1.18)

# ================================================================
# Panel j: RF Permutation Importance (mature miRNA names)
# ================================================================

rf_imp <- importance(ml_res$rf_model)
rf_imp_df <- data.frame(
  miRNA = rownames(rf_imp),
  MDA   = rf_imp[, "MeanDecreaseAccuracy"],
  MDG   = rf_imp[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>%
  mutate(mature = extract_mature(miRNA)) %>%
  arrange(desc(MDA))

cat("Computing permutation threshold...\n")
set.seed(42)
X_perm <- t(proc$g1g2_cpm_combat)
s1c_p  <- intersect(colnames(X_perm), rownames(proc$s1_cpm))
X_perm <- X_perm[, s1c_p]
y_perm <- ifelse(proc$g1g2_meta$Group == "PD", 1, 0)

perm_max <- numeric(100)
for (i in 1:100) {
  rf_s <- randomForest(X_perm, as.factor(sample(y_perm)),
                        ntree = 500, importance = TRUE)
  perm_max[i] <- max(importance(rf_s)[, "MeanDecreaseAccuracy"])
}
thresh <- quantile(perm_max, 0.95)

top_n <- min(20, nrow(rf_imp_df))
perm_plot <- rf_imp_df %>%
  head(top_n) %>%
  mutate(confirmed = MDA > thresh)

p_j <- ggplot(perm_plot, aes(x = reorder(mature, MDA), y = MDA, fill = confirmed)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = thresh, linetype = "dashed",
             color = "#E74C3C", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = col_rf, "FALSE" = "#BDC3C7")) +
  labs(title = "RF Permutation Importance",
       subtitle = sprintf("Threshold (95th pct): %.2f", thresh),
       x = "", y = "Mean Decrease Accuracy") +
  theme_j +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6.5, hjust = 1))

# ================================================================
# Assemble Figure 1
# ================================================================

row1 <- plot_grid(
  plot_grid(p_a, p_b, p_c, ncol = 3, rel_widths = c(1, 1, 1),
            labels = c("a", "b", "c"), label_size = 12,
            label_fontface = "bold", label_x = 0, label_y = 1),
  leg_group,
  ncol = 1, rel_heights = c(1, 0.06)
)

# legend 放在三个 panel 整体下方
row2 <- plot_grid(def_panel, leg_batch,
                   ncol = 1, rel_heights = c(1, 0.08))

row3 <- plot_grid(p_g, p_h, ncol = 2, rel_widths = c(1, 1.15),
                  labels = c("g", "h"), label_size = 12,
                  label_fontface = "bold", label_x = 0, label_y = 1)

row4 <- plot_grid(p_i, p_j, ncol = 2, rel_widths = c(1, 1.4),
                  labels = c("i", "j"), label_size = 12,
                  label_fontface = "bold", label_x = 0, label_y = 1)

fig1 <- plot_grid(
  row1, row2, row3, row4,
  ncol   = 1,
  rel_heights = c(1, 1.1, 1, 1.15)
)

ggsave("figures/publication/Figure1_ML_Analysis.pdf",
       fig1, width = 7.2, height = 11, units = "in", dpi = 300)

cat("\nFigure 1 saved to figures/publication/\n")
