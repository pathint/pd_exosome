# ============================================================
# 03_ML.R —— 多方法特征选择 + 分类模型 + Panel优化（完整修复版）
# ============================================================

library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)
library(pROC)

set.seed(42)
proc <- readRDS("data/processed/processed_data.rds")
meta_res <- read_csv("results/02_DEA/meta_analysis_results.csv")

# ================================================================
# 4.1 准备 ML 数据
# ================================================================
# 策略：GSE269775 + GSE269776（ComBat校正后）做训练，FJMU做独立验证

g1g2_expr   <- t(proc$g1g2_cpm_combat)
g1g2_labels <- proc$g1g2_meta$Group
g1g2_dataset<- proc$g1g2_meta$Dataset

# FJMU：使用三数据集共有 miRNA
s1_common   <- intersect(colnames(g1g2_expr), rownames(proc$s1_cpm))
g1g2_expr   <- g1g2_expr[, s1_common]
s1_expr     <- t(proc$s1_cpm[s1_common, ])
s1_labels   <- proc$s1_meta$Group

y_train <- ifelse(g1g2_labels == "PD", 1, 0)
y_test  <- ifelse(s1_labels == "PD", 1, 0)

cat("Feature space:", length(s1_common), "miRNAs\n")
cat("Train:", length(y_train), "(PD:", sum(y_train == 1), "HC:", sum(y_train == 0), ")\n")
cat("Test:", length(y_test), "(PD:", sum(y_test == 1), "HC:", sum(y_test == 0), ")\n")

# ================================================================
# 4.2 多方法特征选择
# ================================================================

cat("\n========== Feature Selection ==========\n")

# --- 方法1: LASSO Logistic Regression ---
cat("1. LASSO...\n")
cv_lasso <- cv.glmnet(
  x = as.matrix(g1g2_expr), y = y_train,
  family = "binomial", alpha = 1,
  nfolds = min(10, sum(y_train == 0), sum(y_train == 1)),
  type.measure = "auc"
)
lasso_coefs   <- coef(cv_lasso, s = "lambda.min")
lasso_features <- rownames(lasso_coefs)[which(lasso_coefs[-1, ] != 0)]
cat("   LASSO selected:", length(lasso_features), "features\n")

# LASSO 系数路径图
pdf("results/03_ML/lasso_coefficient_path.pdf", width = 8, height = 6)
plot(cv_lasso, main = "LASSO Cross-Validation Curve")
abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)
abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2)
legend("topright", legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2)
dev.off()

# --- 方法2: Random Forest Importance ---
cat("2. Random Forest importance...\n")
rf_model <- randomForest(
  x = g1g2_expr, y = as.factor(y_train),
  ntree = 1000, importance = TRUE,
  mtry = floor(sqrt(ncol(g1g2_expr)))
)
rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(
  miRNA = rownames(rf_imp),
  MeanDecreaseAccuracy = rf_imp[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini     = rf_imp[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>% arrange(desc(MeanDecreaseGini))

# Elbow point 选择 top features
gini_vals <- sort(rf_imp_df$MeanDecreaseGini, decreasing = TRUE)
elbow <- which(diff(diff(gini_vals)) > 0)[1]
if (is.na(elbow) || elbow < 10) elbow <- min(50, length(gini_vals))
rf_features <- rf_imp_df$miRNA[1:elbow]
cat("   RF selected:", length(rf_features), "features\n")

# RF importance plot
p_rf <- rf_imp_df %>% head(30) %>%
  ggplot(aes(x = reorder(miRNA, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "#27AE60", alpha = 0.8) +
  coord_flip() +
  labs(title = "Random Forest -- Top 30 Features (Mean Decrease Gini)",
       x = "", y = "Mean Decrease Gini") +
  theme_bw(base_size = 11)
ggsave("results/03_ML/rf_importance.pdf", p_rf, width = 8, height = 8)

# --- 方法3: RF Permutation Importance（替代 Boruta） ---
cat("3. RF Permutation Importance (Boruta replacement)...\n")

n_perm <- 100
perm_max_MDA <- numeric(n_perm)

for (i in 1:n_perm) {
  y_shuffled <- sample(y_train)
  rf_shuf <- randomForest(
    x = g1g2_expr, y = as.factor(y_shuffled),
    ntree = 500, importance = TRUE
  )
  perm_max_MDA[i] <- max(importance(rf_shuf)[, "MeanDecreaseAccuracy"])
}

threshold_MDA <- quantile(perm_max_MDA, 0.95)

boruta_feats <- rf_imp_df %>%
  filter(MeanDecreaseAccuracy > threshold_MDA) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  pull(miRNA)
cat("   Permutation Importance confirmed:", length(boruta_feats),
    "features (MDA >", sprintf("%.2f", threshold_MDA), ")\n")

# 可视化
p_perm <- rf_imp_df %>%
  mutate(confirmed = MeanDecreaseAccuracy > threshold_MDA) %>%
  arrange(desc(MeanDecreaseAccuracy)) %>%
  head(30) %>%
  ggplot(aes(x = reorder(miRNA, MeanDecreaseAccuracy),
             y = MeanDecreaseAccuracy, fill = confirmed)) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = threshold_MDA, linetype = "dashed", color = "red") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#BDC3C7")) +
  labs(title = "RF Permutation Importance",
       subtitle = sprintf("Threshold (95th pct): %.2f", threshold_MDA),
       x = "", y = "Mean Decrease Accuracy") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave("results/03_ML/rf_permutation_importance.pdf", p_perm, width = 8, height = 8)

# Permutation 分布图
p_perm_dist <- ggplot(data.frame(MDA = perm_max_MDA), aes(x = MDA)) +
  geom_histogram(bins = 25, fill = "#3498DB", alpha = 0.7, color = "white") +
  geom_vline(xintercept = threshold_MDA, color = "red", linewidth = 1) +
  annotate("text", x = threshold_MDA + 0.5, y = Inf, vjust = 2,
           label = sprintf("95th pct = %.2f", threshold_MDA), color = "red") +
  labs(title = "Permutation Distribution of Max MDA",
       x = "Max Mean Decrease Accuracy", y = "Count") +
  theme_bw(base_size = 12)
ggsave("results/03_ML/permutation_distribution.pdf", p_perm_dist, width = 7, height = 5)

# --- 方法4: Univariate (meta-analysis) ---
cat("4. Univariate selection (meta-analysis)...\n")
univar_features <- meta_res %>%
  filter(padj_meta < 0.05, abs(beta_meta) > 0.5) %>%
  pull(miRNA) %>%
  intersect(s1_common)
cat("   Univariate selected:", length(univar_features), "features\n")

# --- 方法5: Elastic Net ---
cat("5. Elastic Net...\n")
cv_enet <- cv.glmnet(
  x = as.matrix(g1g2_expr), y = y_train,
  family = "binomial", alpha = 0.5,
  nfolds = min(10, sum(y_train == 0), sum(y_train == 1)),
  type.measure = "auc"
)
enet_coefs    <- coef(cv_enet, s = "lambda.min")
enet_features <- rownames(enet_coefs)[which(enet_coefs[-1, ] != 0)]
cat("   Elastic Net selected:", length(enet_features), "features\n")

# ================================================================
# 4.3 特征交集与共识（修复版）
# ================================================================

cat("\n========== Feature Consensus ==========\n")

feature_list <- list(
  LASSO      = lasso_features,
  RF         = rf_features,
  PermImp    = boruta_feats,
  Univar     = univar_features,
  ElasticNet = enet_features
)

# UpSet 图
feature_upset <- UpSetR::fromList(feature_list)
pdf("results/03_ML/feature_selection_upset.pdf", width = 10, height = 6)
UpSetR::upset(feature_upset, sets = names(feature_list),
              order.by = "freq", mb.ratio = c(0.6, 0.4),
              text.scale = 1.1)
dev.off()

# 共识：至少 >=2 种方法
fc_table <- table(unlist(feature_list))
fc_vec   <- as.integer(fc_table)
names(fc_vec) <- names(fc_table)

consensus_features <- names(fc_vec[fc_vec >= 2])
cat("Consensus features (>=2 methods):", length(consensus_features), "\n")

# 兜底
if (length(consensus_features) < 3) {
  cat("Warning: Few consensus features. Supplementing from meta-analysis top hits.\n")
  meta_top <- meta_res %>% head(20) %>% pull(miRNA)
  consensus_features <- unique(c(consensus_features, meta_top))
}

feat_summary <- data.frame(
  miRNA = consensus_features,
  stringsAsFactors = FALSE
)
feat_summary$n_methods     <- fc_vec[feat_summary$miRNA]
feat_summary$in_LASSO      <- feat_summary$miRNA %in% lasso_features
feat_summary$in_RF         <- feat_summary$miRNA %in% rf_features
feat_summary$in_PermImp    <- feat_summary$miRNA %in% boruta_feats
feat_summary$in_Univar     <- feat_summary$miRNA %in% univar_features
feat_summary$in_ElasticNet <- feat_summary$miRNA %in% enet_features
feat_summary <- feat_summary[order(feat_summary$n_methods, decreasing = TRUE), ]

write_csv(feat_summary, "results/03_ML/consensus_features.csv")

# ================================================================
# 4.4 分类模型训练
# ================================================================

cat("\n========== Model Training ==========\n")

X_tr <- g1g2_expr[, consensus_features, drop = FALSE]
X_te <- s1_expr[, consensus_features, drop = FALSE]

models <- list()

# --- LR ---
cat("Training Logistic Regression...\n")
train_df <- data.frame(y = as.factor(y_train), X_tr)
test_df  <- data.frame(X_te)

lr_model <- glm(y ~ ., data = train_df, family = binomial())
lr_prob  <- predict(lr_model, newdata = test_df, type = "response")
lr_roc   <- roc(y_test, lr_prob, quiet = TRUE)

models$LR <- list(
  auc  = as.numeric(lr_roc$auc),
  prob = lr_prob,
  class = ifelse(lr_prob > 0.5, 1, 0)
)

# --- SVM ---
cat("Training SVM...\n")
svm_model <- svm(x = X_tr, y = as.factor(y_train),
                  kernel = "radial", probability = TRUE)
svm_pred <- predict(svm_model, X_te, probability = TRUE)
svm_prob <- attr(svm_pred, "probabilities")[, "1"]
svm_roc  <- roc(y_test, svm_prob, quiet = TRUE)

models$SVM <- list(
  auc   = as.numeric(svm_roc$auc),
  prob  = svm_prob,
  class = as.numeric(as.character(svm_pred))
)

# --- RF ---
cat("Training Random Forest...\n")
rf_clf  <- randomForest(x = X_tr, y = as.factor(y_train), ntree = 500)
rf_prob <- predict(rf_clf, X_te, type = "prob")[, "1"]
rf_pred <- predict(rf_clf, X_te)
rf_roc  <- roc(y_test, rf_prob, quiet = TRUE)

models$RF <- list(
  auc   = as.numeric(rf_roc$auc),
  prob  = rf_prob,
  class = as.numeric(as.character(rf_pred))
)

# --- XGBoost ---
cat("Training XGBoost...\n")
xgb_tr <- xgb.DMatrix(data = as.matrix(X_tr), label = y_train)
xgb_te <- xgb.DMatrix(data = as.matrix(X_te), label = y_test)

xgb_params <- list(
  objective = "binary:logistic", eval_metric = "auc",
  max_depth = 3, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8
)
xgb_model <- xgb.train(
  params = xgb_params, data = xgb_tr, nrounds = 100,
  verbose = 0, early_stopping_rounds = 10,
  watchlist = list(train = xgb_tr)
)
xgb_prob <- predict(xgb_model, xgb_te)
xgb_roc  <- roc(y_test, xgb_prob, quiet = TRUE)

models$XGBoost <- list(
  auc   = as.numeric(xgb_roc$auc),
  prob  = xgb_prob,
  class = ifelse(xgb_prob > 0.5, 1, 0)
)

cat("\n--- Independent Test (FJMU) AUC ---\n")
for (m in names(models)) {
  cat(sprintf("  %-8s: AUC = %.3f\n", m, models[[m]]$auc))
}

# ================================================================
# 4.5 5-Fold CV（内部验证，G1+G2）
# ================================================================

cat("\n========== 5-Fold Cross-Validation ==========\n")

cv_folds <- createFolds(y_train, k = 5, returnTrain = TRUE)
cv_aucs <- list(LR = c(), SVM = c(), RF = c(), XGBoost = c())

for (fi in seq_along(cv_folds)) {
  tr_idx <- cv_folds[[fi]]
  vl_idx <- setdiff(seq_along(y_train), tr_idx)

  Xcv_tr <- X_tr[tr_idx, , drop = FALSE]; ycv_tr <- y_train[tr_idx]
  Xcv_vl <- X_tr[vl_idx, , drop = FALSE]; ycv_vl <- y_train[vl_idx]

  tryCatch({
    cv_lr <- cv.glmnet(
      x = as.matrix(Xcv_tr), y = ycv_tr,
      family = "binomial", alpha = 0.01,   # 极小弹性网络惩罚
      nfolds = min(5, min(sum(ycv_tr == 0), sum(ycv_tr == 1))),
      type.measure = "auc"
    )
    p <- predict(cv_lr, as.matrix(Xcv_vl), s = "lambda.min", type = "response")
    cv_aucs$LR <- c(cv_aucs$LR, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$LR <<- c(cv_aucs$LR, NA))

  # SVM
  tryCatch({
    m <- svm(Xcv_tr, as.factor(ycv_tr), kernel = "radial", probability = TRUE)
    p <- attr(predict(m, Xcv_vl, probability = TRUE), "probabilities")[, "1"]
    cv_aucs$SVM <- c(cv_aucs$SVM, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$SVM <<- c(cv_aucs$SVM, NA))

  # RF
  tryCatch({
    m <- randomForest(Xcv_tr, as.factor(ycv_tr), ntree = 300)
    p <- predict(m, Xcv_vl, type = "prob")[, "1"]
    cv_aucs$RF <- c(cv_aucs$RF, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$RF <<- c(cv_aucs$RF, NA))

  # XGBoost
  tryCatch({
    m <- xgb.train(
      params = xgb_params,
      data = xgb.DMatrix(as.matrix(Xcv_tr), label = ycv_tr),
      nrounds = 50, verbose = 0
    )
    p <- predict(m, xgb.DMatrix(as.matrix(Xcv_vl)))
    cv_aucs$XGBoost <- c(cv_aucs$XGBoost, as.numeric(roc(ycv_vl, p, quiet = TRUE)$auc))
  }, error = function(e) cv_aucs$XGBoost <<- c(cv_aucs$XGBoost, NA))
}

cv_summary <- data.frame(
  Model   = names(cv_aucs),
  Mean_AUC = sapply(cv_aucs, function(x) mean(x, na.rm = TRUE)),
  SD_AUC   = sapply(cv_aucs, function(x) sd(x, na.rm = TRUE)),
  Min_AUC  = sapply(cv_aucs, function(x) min(x, na.rm = TRUE)),
  Max_AUC  = sapply(cv_aucs, function(x) max(x, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

cat("\n--- 5-Fold CV Results (G1+G2) ---\n")
print(cv_summary)

# ================================================================
# 4.6 ROC 曲线对比
# ================================================================

roc_df <- bind_rows(
  data.frame(FPR = 1 - lr_roc$specificities,  TPR = lr_roc$sensitivities,
             Model = sprintf("LR (AUC=%.3f)", models$LR$auc)),
  data.frame(FPR = 1 - svm_roc$specificities, TPR = svm_roc$sensitivities,
             Model = sprintf("SVM (AUC=%.3f)", models$SVM$auc)),
  data.frame(FPR = 1 - rf_roc$specificities,  TPR = rf_roc$sensitivities,
             Model = sprintf("RF (AUC=%.3f)", models$RF$auc)),
  data.frame(FPR = 1 - xgb_roc$specificities, TPR = xgb_roc$sensitivities,
             Model = sprintf("XGBoost (AUC=%.3f)", models$XGBoost$auc))
)

# 显式指定方向
lr_roc <- roc(y_test, models$LR$prob, direction = "<", quiet = TRUE)
svm_roc <- roc(y_test, models$SVM$prob, direction = "<", quiet = TRUE)
rf_roc <- roc(y_test, models$RF$prob, direction = "<", quiet = TRUE)
xgb_roc <- roc(y_test, models$XGBoost$prob, direction = "<", quiet = TRUE)

# 使用 ggroc 绘制，自动保证单调性
#library(ggroc)
p_roc <- ggroc(list(
  LR = lr_roc,
  SVM = svm_roc,
  RF = rf_roc,
  XGBoost = xgb_roc
), size = 1.2) +
  geom_abline(linetype = "dashed", color = "grey50") +
  labs(title = "ROC Curves -- Independent Test on FJMU",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  scale_color_manual(values = c("LR" = "#E74C3C", "SVM" = "#3498DB",
                                "RF" = "#27AE60", "XGBoost" = "#F39C12")) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.1, 0.8))

# 添加AUC文本（自动从roc对象获取）
p_roc <- p_roc + 
  annotate("text", x = 0.2, y = 0.1, 
           label = paste0("LR AUC = ", round(auc(lr_roc), 3), "\n",
                          "SVM AUC = ", round(auc(svm_roc), 3), "\n",
                          "RF AUC = ", round(auc(rf_roc), 3), "\n",
                          "XGB AUC = ", round(auc(xgb_roc), 3)),
           hjust = 0, size = 4, color = "black")

ggsave("results/03_ML/AUC_ROC_FJMU.pdf", p_roc, width = 8, height = 7)


# CV AUC bar plot
p_cv <- ggplot(cv_summary, aes(x = reorder(Model, -Mean_AUC), y = Mean_AUC, fill = Model)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, Mean_AUC - SD_AUC),
                    ymax = pmin(1, Mean_AUC + SD_AUC)), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", Mean_AUC)), vjust = -0.8, size = 3.5) +
  scale_fill_manual(values = c("LR" = "#E74C3C", "SVM" = "#3498DB",
                               "RF" = "#27AE60", "XGBoost" = "#F39C12")) +
  labs(title = "5-Fold CV AUC on G1+G2", x = "", y = "AUC") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  ylim(0, 1.05)
ggsave("results/03_ML/CV_AUC_comparison.pdf", p_cv, width = 7, height = 5)

# ================================================================
# 4.7 Panel 大小优化
# ================================================================

cat("\n========== Panel Size Optimization ==========\n")

mirna_ranked <- meta_res %>%
  filter(miRNA %in% s1_common) %>%
  arrange(pval_meta) %>%
  pull(miRNA)

panel_sizes <- c(1, 2, 3, 5, 8, 10, 15, 20, 30, 50)
panel_sizes <- panel_sizes[panel_sizes <= length(mirna_ranked)]

panel_aucs <- map_dfr(panel_sizes, function(k) {
  features_k <- mirna_ranked[1:k]
  X_k <- g1g2_expr[, features_k, drop = FALSE]

  fold_aucs <- c()
  for (fi in seq_along(cv_folds)) {
    tr <- cv_folds[[fi]]
    vl <- setdiff(1:nrow(X_k), tr)

    tryCatch({
      m <- glm(as.factor(y_train[tr]) ~ .,
               data = data.frame(y = as.factor(y_train[tr]), X_k[tr, , drop = FALSE]),
               family = binomial())
      p <- predict(m, data.frame(X_k[vl, , drop = FALSE]), type = "response")
      fold_aucs <- c(fold_aucs, as.numeric(roc(y_train[vl], p, quiet = TRUE)$auc))
    }, error = function(e) {
      fold_aucs <<- c(fold_aucs, NA)
    })
  }
  data.frame(n = k, auc = mean(fold_aucs, na.rm = TRUE),
             sd = sd(fold_aucs, na.rm = TRUE), stringsAsFactors = FALSE)
})

optimal_n <- panel_aucs$n[which.max(panel_aucs$auc)]
cat("Optimal panel size:", optimal_n, "miRNAs (CV AUC =",
    sprintf("%.3f", max(panel_aucs$auc, na.rm = TRUE)), ")\n")

p_panel <- ggplot(panel_aucs, aes(x = n, y = auc)) +
  geom_ribbon(aes(ymin = pmax(0, auc - sd), ymax = pmin(1, auc + sd)),
              fill = "#3498DB", alpha = 0.15) +
  geom_line(color = "#2C3E50", linewidth = 1) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_vline(xintercept = optimal_n, linetype = "dashed", color = "#E74C3C") +
  annotate("text", x = optimal_n + 2, y = min(panel_aucs$auc, na.rm = TRUE),
           label = sprintf("Optimal: %d miRNAs\nAUC = %.3f", optimal_n,
                           max(panel_aucs$auc, na.rm = TRUE)),
           hjust = 0, color = "#E74C3C", size = 3.5) +
  labs(title = "Panel Size Optimization (5-Fold CV, Logistic Regression)",
       x = "Number of miRNAs in Panel", y = "Mean AUC (CV)") +
  theme_bw(base_size = 13)
ggsave("results/03_ML/panel_size_optimization.pdf", p_panel, width = 8, height = 6)

write_csv(panel_aucs, "results/03_ML/panel_size_aucs.csv")

# ================================================================
# 4.8 Bootstrap 稳定性评估
# ================================================================

cat("\n========== RF Bootstrap Stability ==========\n")

n_boot <- 200
boot_aucs      <- c()
boot_sensitivities <- c()
boot_specificities <- c()
boot_feature_freq  <- rep(0, length(consensus_features))
names(boot_feature_freq) <- consensus_features

panel_feats <- consensus_features[1:min(optimal_n, length(consensus_features))]

set.seed(42)

for (b in 1:n_boot) {
  boot_idx <- sample(1:nrow(g1g2_expr), nrow(g1g2_expr), replace = TRUE)
  oob_idx  <- setdiff(1:nrow(g1g2_expr), unique(boot_idx))

  if (length(oob_idx) < 5) next

  X_boot <- g1g2_expr[boot_idx, consensus_features, drop = FALSE]
  X_oob  <- g1g2_expr[oob_idx,  consensus_features, drop = FALSE]
  y_boot <- y_train[boot_idx]
  y_oob  <- y_train[oob_idx]

  # 至少两个类别都有样本
  if (length(unique(y_boot)) < 2 || length(unique(y_oob)) < 2) next

  tryCatch({
    rf_boot <- randomForest(
      x = X_boot, y = as.factor(y_boot),
      ntree = 500, importance = TRUE
    )

    # OOB 预测
    p_oob    <- predict(rf_boot, X_oob, type = "prob")[, "1"]
    pred_oob <- ifelse(p_oob > 0.5, 1, 0)

    # AUC
    auc_val   <- as.numeric(roc(y_oob, p_oob, direction = "<", quiet = TRUE)$auc)
    boot_aucs <- c(boot_aucs, auc_val)

    # Sensitivity / Specificity
    tp <- sum(pred_oob == 1 & y_oob == 1)
    fn <- sum(pred_oob == 0 & y_oob == 1)
    tn <- sum(pred_oob == 0 & y_oob == 0)
    fp <- sum(pred_oob == 1 & y_oob == 0)
    boot_sensitivities <- c(boot_sensitivities, tp / max(tp + fn, 1))
    boot_specificities <- c(boot_specificities, tn / max(tn + fp, 1))

    # 特征选择频率
    imp_boot  <- importance(rf_boot)[, "MeanDecreaseGini"]
    top_feats <- names(sort(imp_boot, decreasing = TRUE))[1:length(panel_feats)]
    boot_feature_freq[top_feats] <- boot_feature_freq[top_feats] + 1

  }, error = function(e) {})
}

n_success <- length(boot_aucs)
cat(sprintf("Successful bootstrap iterations: %d / %d\n", n_success, n_boot))

if (n_success >= 10) {
  # --- 统计汇总 ---
  cat(sprintf("\nRF Bootstrap AUC:          %.3f +/- %.3f", mean(boot_aucs), sd(boot_aucs)))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", quantile(boot_aucs, 0.025), quantile(boot_aucs, 0.975)))
  cat(sprintf("RF Bootstrap Sensitivity:  %.3f +/- %.3f", mean(boot_sensitivities), sd(boot_sensitivities)))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", quantile(boot_sensitivities, 0.025), quantile(boot_sensitivities, 0.975)))
  cat(sprintf("RF Bootstrap Specificity:  %.3f +/- %.3f", mean(boot_specificities), sd(boot_specificities)))
  cat(sprintf("  95%% CI: [%.3f, %.3f]\n", quantile(boot_specificities, 0.025), quantile(boot_specificities, 0.975)))

  # --- 图 4.8A: AUC 分布 ---
  p_boot_auc <- ggplot(data.frame(AUC = boot_aucs), aes(x = AUC)) +
    geom_histogram(bins = 30, fill = "#27AE60", alpha = 0.7, color = "white") +
    geom_vline(xintercept = mean(boot_aucs), color = "#E74C3C", linewidth = 1.2) +
    geom_vline(xintercept = quantile(boot_aucs, c(0.025, 0.975)),
               color = "#E74C3C", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = mean(boot_aucs), y = Inf, vjust = 2,
             label = sprintf("Mean = %.3f\n95%% CI [%.3f, %.3f]",
                             mean(boot_aucs), quantile(boot_aucs, 0.025),
                             quantile(boot_aucs, 0.975)),
             hjust = -0.1, color = "#E74C3C", size = 3.5) +
    labs(title = sprintf("RF Bootstrap AUC Distribution (n=%d)", n_success),
         x = "AUC", y = "Count") +
    theme_bw(base_size = 13)
  ggsave("results/03_ML/RF_bootstrap_AUC.pdf", p_boot_auc, width = 8, height = 5)

  # --- 图 4.8B: Sensitivity + Specificity 分布 ---
  boot_metrics <- data.frame(
    value  = c(boot_sensitivities, boot_specificities),
    metric = rep(c("Sensitivity", "Specificity"), each = n_success)
  )

  p_boot_se_sp <- ggplot(boot_metrics, aes(x = value, fill = metric)) +
    geom_histogram(bins = 25, alpha = 0.6, position = "identity", color = "white") +
    geom_vline(data = boot_metrics %>% group_by(metric) %>%
                 summarise(m = mean(value), .groups = "drop"),
               aes(xintercept = m, color = metric), linewidth = 1.2) +
    scale_fill_manual(values = c("Sensitivity" = "#E74C3C", "Specificity" = "#3498DB")) +
    scale_color_manual(values = c("Sensitivity" = "#C0392B", "Specificity" = "#2471A3")) +
    facet_wrap(~ metric, ncol = 1) +
    labs(title = "RF Bootstrap Sensitivity & Specificity",
         x = "Metric Value", y = "Count") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave("results/03_ML/RF_bootstrap_sensitivity_specificity.pdf",
         p_boot_se_sp, width = 8, height = 8)

  # --- 图 4.8C: 三指标汇总柱状图 ---
  summary_df <- data.frame(
    Metric = c("AUC", "Sensitivity", "Specificity"),
    Mean   = c(mean(boot_aucs), mean(boot_sensitivities), mean(boot_specificities)),
    Lower  = c(quantile(boot_aucs, 0.025), quantile(boot_sensitivities, 0.025),
               quantile(boot_specificities, 0.025)),
    Upper  = c(quantile(boot_aucs, 0.975), quantile(boot_sensitivities, 0.975),
               quantile(boot_specificities, 0.975))
  )

  p_summary <- ggplot(summary_df, aes(x = Metric, y = Mean, fill = Metric)) +
    geom_col(width = 0.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", Mean)), vjust = -0.8, size = 4) +
    scale_fill_manual(values = c("AUC" = "#27AE60", "Sensitivity" = "#E74C3C",
                                  "Specificity" = "#3498DB")) +
    labs(title = "RF Bootstrap Performance Summary",
         subtitle = sprintf("Bar = mean, Error bar = 95%% CI (n=%d)", n_success),
         x = "", y = "Value") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none") +
    ylim(0, 1.05)
  ggsave("results/03_ML/RF_bootstrap_summary.pdf", p_summary, width = 7, height = 5)

  # --- 特征选择频率 ---
  feat_freq_df <- data.frame(
    miRNA     = names(boot_feature_freq),
    frequency = boot_feature_freq / n_success,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(frequency))

  p_feat_freq <- feat_freq_df %>%
    head(20) %>%
    ggplot(aes(x = reorder(miRNA, frequency), y = frequency)) +
    geom_col(fill = "#27AE60", alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", frequency * 100)),
              hjust = -0.1, size = 3.2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    coord_flip() +
    labs(title = "RF Bootstrap: Feature Selection Frequency",
         subtitle = "Proportion of bootstraps where feature ranked in top importance",
         x = "", y = "Selection Frequency") +
    theme_bw(base_size = 11) +
    ylim(0, max(feat_freq_df$frequency) * 1.15)
  ggsave("results/03_ML/RF_bootstrap_feature_frequency.pdf", p_feat_freq, width = 8, height = 7)

  write_csv(feat_freq_df, "results/03_ML/RF_bootstrap_feature_frequency.csv")
}

# --- 保存 bootstrap 结果到 ml_results ---
bootstrap_results <- list(
  aucs              = boot_aucs,
  sensitivities     = boot_sensitivities,
  specificities     = boot_specificities,
  feature_frequency = boot_feature_freq,
  n_success         = n_success
)
# ================================================================
# 4.9 特征重要性汇总
# ================================================================

cat("\n========== Feature Importance Summary ==========\n")

# LASSO 系数
lasso_coef_df <- data.frame(
  miRNA = lasso_features,
  lasso_coef = as.numeric(lasso_coefs[lasso_features, ]),
  stringsAsFactors = FALSE
)

# RF importance
rf_imp_sel <- rf_imp_df %>%
  filter(miRNA %in% consensus_features) %>%
  select(miRNA, RF_MDA = MeanDecreaseAccuracy, RF_MDG = MeanDecreaseGini)

# Elastic Net 系数
enet_coef_df <- data.frame(
  miRNA = enet_features,
  enet_coef = as.numeric(enet_coefs[enet_features, ]),
  stringsAsFactors = FALSE
)

# XGBoost feature importance
xgb_imp <- xgb.importance(model = xgb_model)
xgb_imp_df <- data.frame(
  miRNA = xgb_imp$Feature,
  xgb_gain = xgb_imp$Gain,
  stringsAsFactors = FALSE
)

# 合并
feat_importance <- feat_summary %>%
  select(miRNA, n_methods) %>%
  left_join(lasso_coef_df, by = "miRNA") %>%
  left_join(rf_imp_sel, by = "miRNA") %>%
  left_join(enet_coef_df, by = "miRNA") %>%
  left_join(xgb_imp_df, by = "miRNA")

write_csv(feat_importance, "results/03_ML/feature_importance_summary.csv")

# ================================================================
# 4.10 保存 ML 结果
# ================================================================

ml_results <- list(
  consensus_features = consensus_features,
  feature_table      = feat_summary,
  feat_importance    = feat_importance,
  models             = models,
  cv_summary         = cv_summary,
  cv_aucs            = cv_aucs,
  panel_aucs         = panel_aucs,
  optimal_n          = optimal_n,
  bootstrap_aucs     = boot_aucs,
  lasso_cv           = cv_lasso,
  rf_model           = rf_model,
  xgb_model          = xgb_model,
  rf_bootstrap       = bootstrap_results
)

saveRDS(ml_results, "data/processed/ml_results.rds")

cat("\n========== ML Summary ==========\n")
cat("Consensus features:", length(consensus_features), "\n")
cat("Optimal panel size:", optimal_n, "\n")
cat("Models trained:", paste(names(models), collapse = ", "), "\n")
cat(sprintf("Best CV AUC: %s = %.3f\n",
            cv_summary$Model[which.max(cv_summary$Mean_AUC)],
            max(cv_summary$Mean_AUC)))
cat(sprintf("Independent test AUC: %s = %.3f\n",
            names(models)[which.max(sapply(models, `[[`, "auc"))],
            max(sapply(models, `[[`, "auc"))))
cat("\nML module complete.\n")
