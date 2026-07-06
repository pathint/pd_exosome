# ============================================================
# Figure 4: Multi-omics Integration & Mechanism
# npj Parkinson's Disease style
# Left: 3 rows of panels | Right: heatmap spanning full height
# ============================================================

library(tidyverse)
library(cowplot)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridis)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(grid)

dir.create("figures/publication", recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Constants & Theme
# ================================================================

col_pd  <- "#E74C3C";  col_hc  <- "#3498DB"
col_cns <- "#E74C3C";  col_other <- "#8C8C8C"

theme_j <- theme_bw(base_size = 9) +
  theme(
    text               = element_text(family = "sans"),
    plot.title         = element_text(size = 10, hjust = 0.5,
                                       margin = ggplot2::margin(b = 5)),
    axis.title         = element_text(size = 9),
    axis.text          = element_text(size = 7.5, color = "black"),
    legend.key.size    = unit(0.35, "cm"),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 8),
    legend.margin      = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin  = ggplot2::margin(0, 0, 0, 0, unit = "pt"),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    plot.margin        = ggplot2::margin(3, 5, 3, 3, unit = "pt")
  )

# ================================================================
# Panel a: FJMU mRNA PCA
# ================================================================

cat("\n========== Figure 4 Panel a: FJMU mRNA PCA ==========\n")

fjmu_mrna <- as.matrix(read.csv("data/raw/FJMU_mRNA_counts.csv",
                                  row.names = 1, check.names = FALSE))
fjmu_mrna_meta <- read.csv("data/raw/FJMU_mRNA_meta.csv", row.names = 1,
                            check.names = FALSE)

keep_a <- fjmu_mrna_meta$Group %in% c("PD", "HC")
mat_a <- fjmu_mrna[, keep_a]
grp_a <- droplevels(factor(fjmu_mrna_meta$Group[keep_a], levels = c("HC", "PD")))

pca_a <- prcomp(t(mat_a), scale. = TRUE)
ve_a  <- summary(pca_a)$importance[2, 1:2] * 100
df_a  <- data.frame(PC1 = pca_a$x[, 1], PC2 = pca_a$x[, 2], Group = grp_a)

p_a <- ggplot(df_a, aes(PC1, PC2, color = Group)) +
  geom_point(size = 1.8, alpha = 0.85) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.6) +
  scale_color_manual(values = c("PD" = col_pd, "HC" = col_hc)) +
  labs(title = "FJMU mRNA PCA",
       x = sprintf("PC1 (%.1f%%)", ve_a[1]),
       y = sprintf("PC2 (%.1f%%)", ve_a[2])) +
  theme_j +
  theme(legend.position = c(0.95, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA, size = 0.3))

# ================================================================
# Panel b: KEGG Barplot
# ================================================================

cat("\n========== Figure 4 Panel b: KEGG ==========\n")

kegg_res <- readRDS("data/raw_mRNA/mRNA_kegg_results.rds")

plot_kegg <- fortify(kegg_res, showCategory = 5) %>%
  mutate(Log_P_Adjust = -log10(p.adjust)) %>%
  arrange(desc(Log_P_Adjust))

plot_kegg$Description <- factor(plot_kegg$Description,
                                  levels = rev(plot_kegg$Description))



p_b <- ggplot(plot_kegg, aes(x = Description, y = Log_P_Adjust)) +
  #geom_col(width = 0.55, fill = "#A61E22", alpha = 0.9) +
  geom_col(width = 0.55, fill = "#a80326", alpha = 0.9) +
  geom_text(aes(label = Description, y = 0),
            hjust = 0, nudge_x = 0.05, size = 2.5, fontface = "bold", color = "white") +
  coord_flip() +
  labs(x = NULL, y = expression(-log[10](FDR)),
       title = "KEGG Pathway Enrichment") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic(base_size = 8) +
  theme(
    panel.grid    = element_blank(),
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.text.x   = element_text(size = 7, colour = "black"),
    axis.title.y  = element_blank(),
    axis.title.x  = element_text(size = 8),
    axis.line     = element_line(colour = "black", size = 0.4),
    axis.ticks    = element_line(colour = "black", size = 0.4),
    plot.title    = element_text(size = 9, hjust = 0.5),
    legend.position = "none",
    plot.margin   = ggplot2::margin(5, 5, 5, 5)
  )


# ================================================================
# Panel c: GO BP/CC Barplot
# ================================================================

cat("\n========== Figure 4 Panel c: GO ==========\n")

go_raw <- readRDS("data/raw_mRNA/mRNA_go_results.rds")
go_df <- go_raw@result

plot_go <- go_df %>%
  dplyr::filter(ONTOLOGY %in% c("BP", "CC")) %>%
  mutate(Log_P_Adjust = -log10(p.adjust)) %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Log_P_Adjust), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()

plot_go$ONTOLOGY <- factor(plot_go$ONTOLOGY, levels = c("BP", "CC"))

plot_go <- plot_go %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  ungroup()

p_c <- ggplot(plot_go, aes(x = Description, y = Log_P_Adjust)) +
  geom_col(width = 0.55, fill = "#a80326", alpha = 0.9) +
  geom_text(aes(label = Description, y = 0),
            hjust = 0, nudge_x = 0.05, size = 2.5, fontface = "bold", color = "white") +
  coord_flip() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y",  switch = "y") +
  labs(x = NULL, y = expression(-log[10](FDR)),
       title = "GO Enrichment") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    #strip.text       = element_text(face = "bold", size = 8),
	strip.text.y.left = element_text(face = "bold", size = 8, angle = 0),
	strip.placement   = "outside",
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_text(size = 7, colour = "black"),
    axis.title.y     = element_blank(),
    axis.title.x     = element_text(size = 8),
    axis.line        = element_line(size = 0.4, colour = "black"),
    axis.ticks       = element_line(size = 0.4, colour = "black"),
    panel.grid       = element_blank(),
    plot.title       = element_text(size = 9, hjust = 0.5),
    legend.position  = "none",
    plot.margin      = ggplot2::margin(5, 5, 5, 5)
  )


# ================================================================
# Panel d: EV-origin Cell Composition（ME18 样式）
# ================================================================

cat("\n========== Figure 4 Panel d: EV-origin ==========\n")

cibersort <- readRDS("data/raw_mRNA/CIBERSORT_mRNA_results.rds")

# 逐个 celltype 计算 wilcox p 值
cells <- unique(cibersort$celltype)
wilcox_labels <- data.frame()

for (ct in cells) {
  sub <- cibersort %>% filter(celltype == ct)
  hc_vals <- sub %>% filter(Group == "HC") %>% pull(composition)
  pd_vals <- sub %>% filter(Group == "PD") %>% pull(composition)
  if (length(hc_vals) >= 2 && length(pd_vals) >= 2) {
    wt <- wilcox.test(hc_vals, pd_vals, exact = FALSE)
    p <- wt$p.value
    label <- ifelse(p < 0.001, sprintf("p = %.1e", p),
                    ifelse(p < 0.05, sprintf("p = %.3f", p), "NS"))
    sig <- p < 0.05
  } else {
    label <- "NA"
    sig <- FALSE
  }
  y_pos <- max(sub$composition, na.rm = TRUE)
  wilcox_labels <- rbind(wilcox_labels,
                          data.frame(celltype = ct, label = label, sig = sig,
                                     y_pos = y_pos, stringsAsFactors = FALSE))
}

wilcox_sig <- wilcox_labels %>% filter(sig == TRUE)
wilcox_sig$x_pos <- match(wilcox_sig$celltype, cells)

p_d <- ggplot(cibersort, aes(x = celltype, y = composition, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, position = position_dodge(0.7)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 1., alpha = 1, shape = 16) +
  {if (nrow(wilcox_sig) > 0) list(
    geom_segment(data = wilcox_sig,
                 aes(x = x_pos - 0.18, xend = x_pos - 0.18,
                     y = y_pos, yend = y_pos + 0.02),
                 color = "black", linewidth = 0.4),
    geom_segment(data = wilcox_sig,
                 aes(x = x_pos + 0.18, xend = x_pos + 0.18,
                     y = y_pos, yend = y_pos + 0.02),
                 color = "black", linewidth = 0.4),
    geom_segment(data = wilcox_sig,
                 aes(x = x_pos - 0.18, xend = x_pos + 0.18,
                     y = y_pos + 0.02, yend = y_pos + 0.02),
                 color = "black", linewidth = 0.4),
    geom_text(data = wilcox_sig,
              aes(x = x_pos, y = y_pos + 0.06, label = label),
              color = "black", size = 2.5)
  )} +
  scale_color_manual(values = c("HC" = col_hc, "PD" = col_pd)) +
  labs(title = "Cell Composition (EV-origin)",
       y = "Proportion", x = "") +
  theme_j +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA, size = 0.3)
  ) +
  ylim(NA, max(cibersort$composition, na.rm = TRUE) * 1.25)
# ================================================================
# Panel e–f: UMAP Density（HC/PD 标题无边框）
# ================================================================

cat("\n========== Figure 4 Panels e-f: UMAP Density ==========\n")

cd4_df   <- readRDS("data/raw_mRNA/cd4_df.rds")
Bcell_df <- readRDS("data/raw_mRNA/Bcell_df.rds")

p_e <- ggplot(cd4_df, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                  alpha = 0.7, colour = NA, linewidth = 0.1) +
  facet_wrap(~ Group) +
  scale_fill_viridis_c(option = "viridis", name = "Density") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, face = "bold"),
    legend.key.height = unit(0.8 / 2.54, "cm"),
    legend.key.width  = unit(0.8 / 2.54, "cm"),
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 6),
    plot.title        = element_text(size = 9, hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 7)
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "CD4+ T cell")

p_f <- ggplot(Bcell_df, aes(x = umap_1, y = umap_2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                  alpha = 0.7, colour = NA, linewidth = 0.1) +
  facet_wrap(~ Group) +
  scale_fill_viridis_c(option = "viridis", name = "Density") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, face = "bold"),
    legend.key.height = unit(0.8 / 2.54, "cm"),
    legend.key.width  = unit(0.8 / 2.54, "cm"),
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 6),
    plot.title        = element_text(size = 9, hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 7)
  ) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "B cell")

# ================================================================
# Multi-omics Heatmap（独立保存 + 转为 grob 用于主图）
# ================================================================

cat("\n========== Figure 4 Heatmap ==========\n")

mo_proc <- readRDS("data/processed/multiomics_processed.rds")

fjmu_mirna <- read_csv("results/04_MultiOmics/FJMU_DE_miRNA.csv", show_col_types = FALSE) %>%
  filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(20)
fjmu_mrna_de <- read_csv("results/04_MultiOmics/FJMU_DE_mRNA.csv", show_col_types = FALSE) %>%
  filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(20)
fjmu_prot <- read_csv("results/04_MultiOmics/FJMU_DE_protein.csv", show_col_types = FALSE) %>%
  filter(significance != "NS") %>% arrange(adj.P.Val) %>% head(15)

# 添加 Figure 3 的 4 个 miRNA
extra_mirnas <- c("hsa-miR-130a-3p", "hsa-miR-598-3p",
                   "hsa-miR-101-3p",  "hsa-miR-146a-5p")
mirna_all <- mo_proc$v_mirna$E
extra_in_mat <- intersect(extra_mirnas, rownames(mirna_all))
extra_not_in_de <- setdiff(extra_in_mat, fjmu_mirna$miRNA)
if (length(extra_not_in_de) > 0) {
  cat(sprintf("  Adding %d Figure 3 miRNAs\n", length(extra_not_in_de)))
  fjmu_mirna <- bind_rows(fjmu_mirna,
                           data.frame(miRNA = extra_not_in_de, stringsAsFactors = FALSE))
}

# Z-score
mirna_heat <- t(scale(t(mirna_all[fjmu_mirna$miRNA, ])))
mrna_heat  <- t(scale(t(mo_proc$v_mrna$E[fjmu_mrna_de$gene, ])))
prot_heat  <- t(scale(t(mo_proc$prot_mat[fjmu_prot$protein, ])))
mirna_heat[is.nan(mirna_heat)] <- 0
mrna_heat[is.nan(mrna_heat)] <- 0
prot_heat[is.nan(prot_heat)] <- 0

# 统一样本
sample_order <- mo_proc$s1_meta_4o$ID[order(mo_proc$s1_meta_4o$Group)]
sample_order <- sample_order[sample_order %in% colnames(mirna_heat)]
sample_order <- sample_order[sample_order %in% colnames(mrna_heat)]
sample_order <- sample_order[sample_order %in% colnames(prot_heat)]
mirna_heat <- mirna_heat[, sample_order]
mrna_heat  <- mrna_heat[, sample_order]
prot_heat  <- prot_heat[, sample_order]

# 提取 mature miRNA name（最后一个 "_" 之后的部分）
rownames(mirna_heat) <- sapply(strsplit(rownames(mirna_heat), "_"), function(x) {
							     if (length(x) >= 2) tail(x, 1) else x[1]
						   })

cat(sprintf("  miRNA: %d, mRNA: %d, Protein: %d, Samples: %d\n",
            nrow(mirna_heat), nrow(mrna_heat), nrow(prot_heat), length(sample_order)))

# 顶部注释
anno_col <- HeatmapAnnotation(
  Group = mo_proc$s1_meta_4o$Group[match(sample_order, mo_proc$s1_meta_4o$ID)],
  col = list(Group = c(PD = col_pd, HC = col_hc)),
  show_legend = TRUE, show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    Group = list(title = "Group", title_gp = gpar(fontsize = 8),
                 labels_gp = gpar(fontsize = 7)))
)

ht_mirna <- Heatmap(mirna_heat, name = "miRNA",
                     col = colorRamp2(c(-2, 0, 2), c("#3498DB", "white", "#E74C3C")),
                     cluster_columns = TRUE, cluster_rows = TRUE,
                     show_column_names = TRUE, show_row_names = TRUE,
					 show_row_dend = FALSE, show_column_dend = FALSE,
                     row_names_gp = gpar(fontsize = 7),
                     column_names_gp = gpar(fontsize = 6),
                     #column_title = "Exosomal miRNA",
                     #column_title_gp = gpar(fontsize = 9, fontface = "bold"),
					 row_title = "Exosomal miRNA",
					 row_title_side = "left",
					 row_title_rot = 90,
					 row_title_gp = gpar(fontsize = 9),
                     top_annotation = anno_col,
                     #width = unit(8, "cm"),
                     heatmap_legend_param = list(title = "Z-score",
                                                  title_gp = gpar(fontsize = 8),
                                                  labels_gp = gpar(fontsize = 7),
                                                  legend_height = unit(2, "cm")))

ht_mrna <- Heatmap(mrna_heat, name = "mRNA",
                    col = colorRamp2(c(-2, 0, 2), c("#27AE60", "white", "#F39C12")),
                    cluster_columns = TRUE, cluster_rows = TRUE,
                    show_column_names = TRUE, show_row_names = TRUE,
					 show_row_dend = FALSE, show_column_dend = FALSE,
                    row_names_gp = gpar(fontsize = 6),
                    #column_names_gp = gpar(fontsize = 6),
                    #column_title = "mRNA",
                    #column_title_gp = gpar(fontsize = 9, fontface = "bold"),
					 row_title = "Exosomal mRNA",
					 row_title_side = "left",
					 row_title_rot = 90,
					 row_title_gp = gpar(fontsize = 9),
                    #width = unit(8, "cm"),
                    heatmap_legend_param = list(title = "Z-score",
                                                title_gp = gpar(fontsize = 8),
                                                labels_gp = gpar(fontsize = 7),
                                                legend_height = unit(2, "cm")))

ht_prot <- Heatmap(prot_heat, name = "Protein",
                    col = colorRamp2(c(-2, 0, 2), c("#8E44AD", "white", "#E67E22")),
                    cluster_columns = TRUE, cluster_rows = TRUE,
                    show_column_names = FALSE, show_row_names = TRUE,
					 show_row_dend = FALSE, show_column_dend = FALSE,
                    row_names_gp = gpar(fontsize = 7),
                    #column_names_gp = gpar(fontsize = 6),
                    #column_title = "Protein",
                    #column_title_gp = gpar(fontsize = 9, fontface = "bold"),
					 row_title = "Exosomal Protein",
					 row_title_side = "left",
					 row_title_rot = 90,
					 row_title_gp = gpar(fontsize = 9),
                    #width = unit(8, "cm"),
                    heatmap_legend_param = list(title = "Z-score",
                                                title_gp = gpar(fontsize = 8),
                                                labels_gp = gpar(fontsize = 7),
                                                legend_height = unit(2, "cm")))

ht_list <- ht_mirna %v% ht_mrna %v% ht_prot

# 独立保存热图
heatmap_height <- max(12, (nrow(mirna_heat) + nrow(mrna_heat) + nrow(prot_heat)) * 0.25)
pdf("figures/publication/Figure4b_heatmap.pdf",
    width = 2.5, height = 6.5, useDingbats = FALSE)
draw(ht_list, gap = unit(3, "mm"),
     column_title = " FJMU \nMulti-Omics DE Heatmap (Z-score)",
     column_title_gp = gpar(fontsize = 11))
dev.off()

# 将热图转为 grob 用于主图嵌入
hm_grob <- grid.grabExpr(
  draw(ht_list, gap = unit(3, "mm"),
       column_title = "FJMU\nMulti-omics DE Heatmap",
       column_title_gp = gpar(fontsize = 10)),
  width = 2.5, height = 6.5
)

cat("  Heatmap saved and captured as grob\n")

# ================================================================
# Assemble Figure 4（左三行 + 右侧热图贯穿）
# ================================================================

cat("\n========== Assembling Figure 4 ==========\n")

# --- 左侧三行 ---
left_row1 <- plot_grid(p_a, p_d, ncol = 2, rel_widths = c(0.8, 1.3),
					    rel_heights = c(0.8,1),
                        labels = c("a", "b"), label_size = 12,
                        label_fontface = "bold", label_x = 0, label_y = 1)

left_row2 <- plot_grid(p_b, p_c, ncol = 2, rel_widths = c(1, 1),
                        labels = c("d", "e"), label_size = 12,
                        label_fontface = "bold", label_x = 0, label_y = 1)

left_row3 <- plot_grid(p_e, p_f, ncol = 2, rel_widths = c(1, 1),
                        labels = c("f", "g"), label_size = 12,
                        label_fontface = "bold", label_x = 0, label_y = 1)

left_panel <- plot_grid(left_row1, left_row2, left_row3,
                          ncol = 1, rel_heights = c(1, 1, 1))

fig4 <- plot_grid(left_panel, hm_grob,
                   ncol = 2, rel_widths = c(2, 1),
                   labels = c("", "c"), label_size = 12,
                   label_fontface = "bold", label_x = 0, label_y = 1)


ggsave("figures/publication/Figure4_MultiOmics_Mechanism.pdf",
       fig4, width = 7.2, height = 6.5, units = "in", dpi = 300)

cat("\nFigure 4 saved to figures/publication/\n")
cat("  Figure4_MultiOmics_Mechanism.pdf (left 2/3, right 1/3 blank)\n")
cat("  Figure4b_heatmap.pdf (independent heatmap, to be composited)\n")
