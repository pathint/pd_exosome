# ============================================================
# Figure 3: miRNA Tissue Specificity + Protein WGCNA
# npj Parkinson's Disease style (updated)
# ============================================================

library(tidyverse)
library(cowplot)
library(ggplot2)
library(enrichplot)
library(ggpubr)


dir.create("figures/publication", recursive = TRUE, showWarnings = FALSE)

# ================================================================
# Constants & Theme
# ================================================================

col_pd   <- "#E74C3C";  col_hc   <- "#3498DB"
col_cns  <- "#E74C3C";  col_other <- "#8C8C8C"

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
# Panel a: miRNA Tissue Specificity (4 miRNAs combined, single panel)
# ================================================================

cat("\n========== Figure 3 Panel a: miRNA Tissue Specificity ==========\n")

read_tissue_data <- function(filepath) {
  if (grepl("\\.txt$", filepath)) {
    df <- read.delim(filepath, sep = "\t", stringsAsFactors = FALSE)
  } else {
    df <- read.csv(filepath, sep = ",", stringsAsFactors = FALSE)
  }
  colnames(df) <- trimws(colnames(df))
  return(df)
}

mirna_files <- c(
  "results/03_ML/miRNA/hsa-miR-130a-3p_Mean_Top20.csv",
  "results/03_ML/miRNA/hsa-miR-598-3p_Mean_Top20.csv",
  "results/03_ML/miRNA/hsa-miR-101-3p_Mean_Top20.csv",
  "results/03_ML/miRNA/hsa-miR-146a-5p_Mean_Top20.csv"
)


mirna_labels <- c("miR-130a-3p", "miR-598-3p", "miR-101-3p", "miR-146a-5p")

tissue_plots <- list()

for (i in seq_along(mirna_files)) {
  cat(sprintf("  Loading: %s\n", basename(mirna_files[i])))
  df <- read_tissue_data(mirna_files[i])
  df$Expression <- as.numeric(df$Expression)
  df <- df[!is.na(df$Expression), ]
  df$Label <- paste0(df$Organ_system, " | ", df$Tissue)

  # 同名 Label 取最大值
  df <- df %>%
    group_by(Label, CNS) %>%
    summarise(Expression = max(Expression, na.rm = TRUE), .groups = "drop")

  # 降序排列：最高在上
  df <- df %>% arrange(desc(Expression)) %>% head(10)
  df$Label <- factor(df$Label, levels = rev(df$Label))

  # 颜色
  #df$fill_col <- ifelse(df$CNS == "CNS", col_cns, col_other)
  #print(df)

  # 颜色：CNS = 红色，外周血 = 蓝色，其他 = 灰色
  df$fill_col <- case_when(
    df$CNS == "CNS" ~ col_cns,
    grepl("immune|blood|marrow", df$Label, ignore.case = TRUE) ~ "#3498DB",
    TRUE ~ col_other
  )

  p <- ggplot(df, aes(x = Label, y = Expression, fill = fill_col)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = Label, y = 0),
              hjust = 0, nudge_y = 0.3, size = 2.2, color = "white", fontface = "bold") +
    coord_flip() +
    scale_fill_identity() +
    labs(title = mirna_labels[i], x = "", y = "Mean Expression (RPMM)") +
    theme_j +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title   = element_text(size = 9, face = "bold")
    )
  tissue_plots[[i]] <- p
}

# 2x2 组合
#p_a <- plot_grid(plotlist = tissue_plots, ncol = 2)

# 简洁图例
p_leg <- ggplot(data.frame(
    Category = factor(c("CNS", "Blood/Immune", "Other"),
                       levels = c("CNS", "Blood/Immune", "Other")),
    x = c(1, 2, 3)),
  aes(x = x, y = 1, fill = Category)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("CNS" = col_cns, "Blood/Immune" = "#3498DB",
                                "Other" = col_other)) +
  theme_j +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

leg_tissue <- get_legend(p_leg)

p_a <- plot_grid(
  plot_grid(plotlist = tissue_plots, ncol = 2),
  leg_tissue,
  ncol = 1, rel_heights = c(1, 0.06)
)
# ================================================================
# Panel b: ME18 Eigengene Boxplot
# ================================================================

cat("\n========== Figure 3 Panel b: ME18 Eigengene ==========\n")

eig_df <- read.csv("results/08_Protein/Enrichment/eigMM18_df.csv",
                    stringsAsFactors = FALSE)
colnames(eig_df) <- trimws(colnames(eig_df))

if ("Group" %in% colnames(eig_df)) {
  eig_df$Group <- factor(eig_df$Group, levels = c("HC", "PD"))
} else if ("Disease" %in% colnames(eig_df)) {
  eig_df$Group <- factor(ifelse(eig_df$Disease == 1, "PD", "HC"),
                          levels = c("HC", "PD"))
}

me_col <- grep("^ME|^me", colnames(eig_df), value = TRUE)
if (length(me_col) == 0) me_col <- "ME"

# Wilcox p 值
wt <- wilcox.test(eig_df[[me_col]] ~ eig_df$Group)
p_val <- wt$p.value
p_label <- ifelse(p_val < 0.001, sprintf("p = %.1e", p_val),
                   sprintf("p = %.3f", p_val))

y_max <- max(eig_df[[me_col]], na.rm = TRUE)
y_min <- min(eig_df[[me_col]], na.rm = TRUE)
y_range <- y_max - y_min

p_b <- ggplot(eig_df, aes(x = Group, y = .data[[me_col]], color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  # bracket + p value
  annotate("segment", x = 1, xend = 1,
           y = y_max + y_range * 0.05, yend = y_max + y_range * 0.1,
           linewidth = 0.4, color = "black") +
  annotate("segment", x = 2, xend = 2,
           y = y_max + y_range * 0.05, yend = y_max + y_range * 0.1,
           linewidth = 0.4, color = "black") +
  annotate("segment", x = 1, xend = 2,
           y = y_max + y_range * 0.1, yend = y_max + y_range * 0.1,
           linewidth = 0.4, color = "black") +
  annotate("text", x = 1.5, y = y_max + y_range * 0.15,
           label = p_label, size = 3, color = "black") +
  scale_color_manual(values = c("HC" = col_hc, "PD" = col_pd)) +
  labs(title = "ME18 Eigengene",
       x = "", y = "Eigengene") +
  theme_j +
  theme(legend.position = "none") +
  ylim(y_min - y_range * 0.1, y_max + y_range * 0.25)


# ================================================================
# Panel c: GSEA Plot
# ================================================================

cat("\n========== Figure 3 Panel c: GSEA ==========\n")

kegg_kk1 <- readRDS("results/08_Protein/Enrichment/protein_gsea_results.rds")

target_id   <- kegg_kk1$ID[1]
target_desc <- kegg_kk1$Description[1]
target_nes  <- kegg_kk1$NES[1]
target_fdr  <- kegg_kk1$p.adjust[1]

cat(sprintf("  Pathway: %s\n", target_desc))
cat(sprintf("  NES = %.3f, FDR = %.2e\n", target_nes, target_fdr))

p_gsea <- gseaplot2(kegg_kk1,
                     geneSetID = target_id,
                     color = "red",
                     base_size = 8,
                     rel_heights = c(1.5, 0.5, 1),
                     subplots = c(1, 2),
                     ES_geom = "line",
                     pvalue_table = FALSE)

p_c <- p_gsea +
  patchwork::plot_annotation(
    title    = target_desc,
    subtitle = sprintf("NES = %.3f | FDR = %.2e", target_nes, target_fdr),
    theme    = theme(
      plot.title    = element_text(size = 9, face = "bold", hjust = 1),
      plot.subtitle = element_text(size = 8, hjust = 1, color = "#555555"),
      plot.margin   = ggplot2::margin(25, 5, 12, 5, unit = "pt")
    )
  )

# ================================================================
# Row 2 组装 + 最终组装（调整比例避免裁切）
# ================================================================

row2 <- plot_grid(p_b, p_c,
                   ncol = 2, rel_widths = c(1, 2.3),
                   labels = c("b", "c"), label_size = 12,
                   label_fontface = "bold", label_x = 0, label_y = 1)

fig3 <- plot_grid(p_a, row2,
                   ncol = 1, rel_heights = c(1.3, 1.2),
                   labels = c("a", NA), label_size = 12,
                   label_fontface = "bold", label_x = 0, label_y = 1)

ggsave("figures/publication/Figure3_MiRNA_Tissue_WGCNA.pdf",
       fig3, width = 7.2, height = 7.5, units = "in", dpi = 300)
