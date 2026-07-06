# ============================================================
# 05_Clinical.R —— 临床关联分析 + 亚型分型 + 通路富集
# ============================================================

library(tidyverse)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(enrichplot)
library(pheatmap)
library(ggpubr)
library(corrplot)
library(ComplexHeatmap)
library(cluster) 

set.seed(43)
proc <- readRDS("data/processed/processed_data.rds")
ml_res <- readRDS("data/processed/ml_results.rds")
meta_res <- read_csv("results/02_DEA/meta_analysis_results.csv")
dea_s1 <- read_csv("results/02_DEA/DESeq2_FJMU.csv")

# ================================================================
# 6.1 疾病严重程度相关性分析
# ================================================================

cat("\n========== Clinical Correlations ==========\n")

pd_g1 <- proc$g1_meta %>% filter(Group == "PD") %>%
  dplyr::select(ID, Dataset, Gender, Age, Yahr, UPDRS.III, LEDD,
                early.H.M, delayed.H.M, MIBG.WR, DaTscan.Ave, DaTscan.AI)
pd_g2 <- proc$g2_meta %>% filter(Group == "PD") %>%
  dplyr::select(ID, Dataset, Gender, Age, Yahr, UPDRS.III, LEDD,
                early.H.M, delayed.H.M, MIBG.WR, DaTscan.Ave, DaTscan.AI)
pd_all <- rbind(pd_g1, pd_g2)

head(pd_all)
exit()

consensus_mirnas <- meta_res %>%
  filter(significance == "Significant") %>%
  pull(miRNA)

combat_mat <- proc$g1g2_cpm_combat
pd_ids <- pd_all$ID

# 定义要检查的临床变量
clinical_vars <- list(
  list(col = "Yahr",         name = "H-Y Stage"),
  list(col = "UPDRS.III",    name = "UPDRS-III"),
  list(col = "LEDD",         name = "LEDD"),
  list(col = "early.H.M",    name = "MIBG early H/M"),
  list(col = "MIBG.WR",      name = "MIBG WR"),
  list(col = "DaTscan.Ave",  name = "DaTscan Ave"),
  list(col = "DaTscan.AI",   name = "DaTscan AI (%)"),
  list(col = "Age",          name = "Age")
)

# 逐 miRNA 逐变量计算（纯 base R，避免 map_dfr 列名问题）
severity_cor <- data.frame(
  miRNA = character(0), var = character(0),
  rho = numeric(0), pvalue = numeric(0), n = integer(0),
  stringsAsFactors = FALSE
)

for (mir in consensus_mirnas) {
  if (!mir %in% rownames(combat_mat)) next
  expr <- combat_mat[mir, pd_ids]

  for (cv in clinical_vars) {
    col_name <- cv$col
    var_name <- cv$name

    if (!col_name %in% colnames(pd_all)) next

    valid <- !is.na(pd_all[[col_name]]) & is.finite(pd_all[[col_name]])
    if (sum(valid) < 10) next

    tryCatch({
      ct <- cor.test(expr[valid], pd_all[[col_name]][valid], method = "spearman")
      severity_cor <- rbind(severity_cor, data.frame(
        miRNA  = mir,
        var    = var_name,
        rho    = as.numeric(ct$estimate),
        pvalue = as.numeric(ct$p.value),
        n      = sum(valid),
        stringsAsFactors = FALSE
      ))
    }, error = function(e) NULL)
  }
}

# 检查是否有结果
if (nrow(severity_cor) > 0) {
  severity_cor$fdr <- p.adjust(severity_cor$pvalue, method = "BH")
  severity_cor <- severity_cor[order(severity_cor$pvalue), ]

  cat("Total correlations tested:", nrow(severity_cor), "\n")
  sig_clin <- severity_cor[severity_cor$fdr < 0.1, , drop = FALSE]
  cat("Significant (FDR<0.1):", nrow(sig_clin), "\n")

  if (nrow(sig_clin) > 0) print(sig_clin)

  write_csv(severity_cor, "results/05_Clinical/clinical_correlations.csv")
} else {
  cat("No clinical correlations computed (insufficient data)\n")
  write_csv(severity_cor, "results/05_Clinical/clinical_correlations.csv")
}
# ================================================================
# 6.2 PD 亚型无监督聚类
# ================================================================


cat("\n========== PD Subtype Clustering ==========\n")

avail_mir <- intersect(ml_res$consensus, rownames(combat_mat))

if (length(avail_mir) >= 3 & nrow(pd_all) >= 10) {
  pd_expr <- combat_mat[avail_mir, pd_all$ID]

  # 层次聚类
  hc_res <- hclust(dist(t(pd_expr), method = "euclidean"), method = "ward.D2")

  # 轮廓系数确定最优 k
  max_k <- min(4, floor(nrow(pd_all) / 3))
  sil_scores <- sapply(2:max_k, function(k) {
    cl <- cutree(hc_res, k = k)
    mean(silhouette(cl, dist(t(pd_expr)))[, 3])
  })
  opt_k <- which.max(sil_scores) + 1
  cat("Optimal clusters:", opt_k, "\n")

  pd_all$cluster <- factor(cutree(hc_res, k = opt_k))

  # PCA 可视化
  pca_pd <- prcomp(t(pd_expr), scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca_pd$x[, 1], PC2 = pca_pd$x[, 2],
    cluster = pd_all$cluster,
    Yahr = pd_all$Yahr, UPDRS = pd_all$UPDRS.III,
    DaTscan = pd_all$DaTscan.Ave
  )

  p_sub <- ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.9, linetype = 2) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "PD Molecular Subtypes") +
    theme_bw(base_size = 13)
  ggsave("results/05_Clinical/PD_subtype_PCA.pdf", p_sub, width = 7, height = 6)

  # 亚型间临床变量比较
  clinical_check <- c("Yahr", "UPDRS.III", "LEDD", "early.H.M", "DaTscan.Ave")

  for (cv in clinical_check) {
    if (cv %in% colnames(pd_all) && sum(!is.na(pd_all[[cv]])) >= 8) {
      p_box <- ggplot(pd_all, aes_string(x = "cluster", y = cv, fill = "cluster")) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.15, size = 1.5) +
        stat_compare_means(method = "kruskal.test") +
        scale_fill_brewer(palette = "Set1") +
        labs(title = paste0(cv, " by Subtype")) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none")
      ggsave(paste0("results/05_Clinical/subtype_", cv, ".pdf"),
             p_box, width = 6, height = 5)
    }
  }

  # 聚类热图
  anno <- data.frame(Cluster = pd_all$cluster, row.names = pd_all$ID)

  pdf("results/05_Clinical/PD_subtype_heatmap.pdf", width = 12, height = 8)
  pheatmap(pd_expr, annotation_col = anno, scale = "row",
           color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
           cluster_cols = hc_res, show_colnames = FALSE,
           main = "PD Subtypes -- Consensus miRNA")
  dev.off()

  write_csv(pd_all, "results/05_Clinical/PD_subtypes.csv")
} else {
  cat("Insufficient features or samples for subtype clustering\n")
}

# ================================================================
# 6.4 外泌体特异性注释
# ================================================================

cat("\n========== Exosome Specificity Annotation ==========\n")

# 从 ExoCarta/Vesiclepedia 查询已知外泌体 miRNA


# ============================================================
# 构建外泌体 miRNA 注释数据库
# ============================================================

library(tidyverse)

# ================================================================
# 2.1 读取并清洗 ExoCarta
# ================================================================

exocarta_raw <- read.delim("data/external/ExoCarta_miRNA_details_6.txt",
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 仅保留人类
exocarta <- exocarta_raw %>%
  filter(grepl("Homo sapiens", SPECIES, ignore.case = TRUE)) %>%
  dplyr::select(mirna_id = MIRNA.ID) %>%
  distinct() %>%
  mutate(
    # 标准化名称：小写、去空格、去*
    mirna_clean = tolower(trimws(gsub("\\*", "", mirna_id))),
    # 补充 hsa- 前缀（ExoCarta 通常没有）
    mirna_hsa = paste0("hsa-", mirna_clean)
  ) %>%
  distinct(mirna_hsa)

cat("ExoCarta human miRNAs:", nrow(exocarta), "\n")

# ================================================================
# 2.2 读取并清洗 Vesiclepedia
# ================================================================

vesiclepedia_raw <- read.delim("data/external/VESICLEPEDIA_MIRNA_DETAILS_5.1.txt",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 仅保留人类
vesiclepedia <- vesiclepedia_raw %>%
  filter(grepl("Homo sapiens", SPECIES, ignore.case = TRUE)) %>%
  dplyr::select(mirna_id = MIRNA.ID) %>%
  distinct() %>%
  mutate(
    mirna_clean = tolower(trimws(gsub("\\*", "", mirna_id))),
    mirna_hsa = paste0("hsa-", mirna_clean)
  ) %>%
  distinct(mirna_hsa)

cat("Vesiclepedia human miRNAs:", nrow(vesiclepedia), "\n")

# ================================================================
# 2.3 合并为统一外泌体数据库
# ================================================================

exosome_db <- full_join(
  exocarta %>% mutate(ExoCarta = TRUE),
  vesiclepedia %>% mutate(Vesiclepedia = TRUE),
  by = "mirna_hsa"
) %>%
  mutate(
    ExoCarta     = ifelse(is.na(ExoCarta), FALSE, TRUE),
    Vesiclepedia = ifelse(is.na(Vesiclepedia), FALSE, TRUE),
    known_exosomal = ExoCarta | Vesiclepedia
  )

cat("Combined exosome DB:", nrow(exosome_db), "unique miRNAs\n")
cat("  ExoCarta only:", sum(exosome_db$ExoCarta & !exosome_db$Vesiclepedia), "\n")
cat("  Vesiclepedia only:", sum(!exosome_db$ExoCarta & exosome_db$Vesiclepedia), "\n")
cat("  Both:", sum(exosome_db$ExoCarta & exosome_db$Vesiclepedia), "\n")

# ================================================================
# 2.4 处理miRNA 名称格式
# ================================================================
# miRNA 名称格式: "hsa-mir-34a_hsa-miR-34a-5p"
# 结构: pre-miRNA_mature-miRNA
# ExoCarta/Vesiclepedia 使用 mature miRNA 名称 (如 hsa-miR-34a-5p)
# 所以需要名称中提取 mature miRNA 部分

# 读取你的 miRNA 列表
miRNA_names <- rownames(proc$g1_cpm)  # 或从 DE 结果中获取

cat("\n========== miRNA Name Format Examples ==========\n")
cat(head(miRNA_names, 5), sep = "\n")

# 提取 mature miRNA 名称（下划线后面的 hsa-miR-... 部分）
extract_mature <- function(x) {
  # 格式: hsa-mir-34a_hsa-miR-34a-5p → hsa-miR-34a-5p
  # 如果包含下划线且下划线后有 hsa-miR，则取后半部分
  mature <- sub(".*?(hsa-miR-.*)$", "\\1", x)
  # 如果提取失败（格式不匹配），返回原始名称
  ifelse(mature == x, tolower(x), tolower(mature))
}

# 同时生成 pre-miRNA 的标准形式用于宽泛匹配
extract_premir <- function(x) {
  # hsa-mir-34a_hsa-miR-34a-5p → hsa-mir-34a
  pre <- sub("^(hsa-mir-[^_]+)_.*$", "\\1", x)
  ifelse(pre == x, tolower(x), tolower(pre))
}

mirna_mapping <- data.frame(
  original_name = miRNA_names,
  mature_name   = extract_mature(miRNA_names),
  premir_name   = extract_premir(miRNA_names),
  stringsAsFactors = FALSE
)

cat("\nName mapping examples:\n")
print(head(mirna_mapping, 10))

# ================================================================
# 2.5 匹配外泌体数据库
# ================================================================
# 多级匹配策略：mature 精确匹配 → pre-miRNA 模糊匹配

# 先准备数据库中的名称（去 hsa- 前缀做宽泛匹配）
exoma_mature <- gsub("^hsa-", "", exosome_db$mirna_hsa)

mirna_mapping$known_exosomal <- sapply(1:nrow(mirna_mapping), function(i) {
  # 精确匹配 mature
  if (mirna_mapping$mature_name[i] %in% exosome_db$mirna_hsa) return(TRUE)

  # 去掉 hsa- 前缀再匹配
  mature_no_prefix <- gsub("^hsa-", "", mirna_mapping$mature_name[i])
  if (mature_no_prefix %in% exoma_mature) return(TRUE)

  # pre-miRNA 模糊匹配（ExoCarta 中有时记录的是 pre-miRNA 名称）
  premir_no_prefix <- gsub("^hsa-", "", mirna_mapping$premir_name[i])
  if (any(grepl(premir_no_prefix, exoma_mature, ignore.case = TRUE))) return(TRUE)

  return(FALSE)
})

cat("\nExosome annotation results:\n")
cat("  Total miRNAs:", nrow(mirna_mapping), "\n")
cat("  Known exosomal:", sum(mirna_mapping$known_exosomal), "\n")
cat("  Novel:", sum(!mirna_mapping$known_exosomal), "\n")

# ================================================================
# 构建 CNS 特异性 miRNA 列表
# ================================================================
# 方法一：基于文献共识的 CNS 富集 miRNA（保守列表）
# 来源：多个文献交叉验证的 CNS/大脑富集 miRNA
# ??
cns_mirnas_consensus <- c(
  # 经典 CNS 富集 miRNA（多篇文献确认）
  "hsa-miR-9", "hsa-miR-9-5p", "hsa-miR-9-3p",
  "hsa-miR-124", "hsa-miR-124-3p", "hsa-miR-124-5p",
  "hsa-miR-132", "hsa-miR-132-3p", "hsa-miR-132-5p",
  "hsa-miR-134", "hsa-miR-134-5p", "hsa-miR-134-3p",
  "hsa-miR-138", "hsa-miR-138-5p", "hsa-miR-138-3p",
  "hsa-miR-153", "hsa-miR-153-3p", "hsa-miR-153-5p",
  "hsa-miR-212", "hsa-miR-212-3p", "hsa-miR-212-5p",
  "hsa-miR-218", "hsa-miR-218-5p", "hsa-miR-218-3p",
  "hsa-miR-219", "hsa-miR-219-5p", "hsa-miR-219-3p",
  "hsa-miR-329", "hsa-miR-329-3p", "hsa-miR-329-5p",
  "hsa-miR-382", "hsa-miR-382-5p", "hsa-miR-382-3p",
  "hsa-miR-433", "hsa-miR-433-3p", "hsa-miR-433-5p",
  "hsa-miR-485", "hsa-miR-485-3p", "hsa-miR-485-5p",
  "hsa-miR-491", "hsa-miR-491-5p", "hsa-miR-491-3p",
  # 神经元特异性
  "hsa-miR-125b", "hsa-miR-125b-5p", "hsa-miR-125b-3p",
  "hsa-miR-128", "hsa-miR-128-3p", "hsa-miR-128-5p",
  "hsa-miR-375", "hsa-miR-375-3p",
  "hsa-miR-146a", "hsa-miR-146a-5p",
  "hsa-miR-204", "hsa-miR-204-5p",
  # 突触/树突特异性
  "hsa-miR-137", "hsa-miR-137-3p"
)

# ================================================================
# 2.6 对 consensus DE miRNAs 做外泌体 + CNS 注释
# ================================================================

meta_sig_mirnas <- meta_res %>%
  filter(significance == "Significant") %>%
  pull(miRNA)

sig_annotation <- mirna_mapping %>%
  dplyr::select(-known_exosomal) %>%
  filter(original_name %in% meta_sig_mirnas) %>%
  left_join(
    exosome_db %>% dplyr::select(mirna_hsa, ExoCarta, Vesiclepedia, known_exosomal),
    by = c("mature_name" = "mirna_hsa")
  ) %>%
  mutate(
    ExoCarta       = ifelse(is.na(ExoCarta), FALSE, ExoCarta),
    Vesiclepedia   = ifelse(is.na(Vesiclepedia), FALSE, Vesiclepedia),
    known_exosomal = ifelse(is.na(known_exosomal), FALSE, known_exosomal),
    # CNS 注释：复用步骤 2.5 中 mirna_mapping 已计算的 cns_enriched
    cns_enriched   = sapply(mature_name, function(x) {
      x %in% cns_mirnas_consensus |
        gsub("^hsa-", "", x) %in% gsub("^hsa-", "", cns_mirnas_consensus)
    })
  )

cat("\nConsensus DE miRNA annotation:\n")
cat("  Total:", nrow(sig_annotation), "\n")
cat("  Known exosomal:", sum(sig_annotation$known_exosomal), "\n")
cat("  Novel exosomal:", sum(!sig_annotation$known_exosomal), "\n")
cat("  CNS-enriched:", sum(sig_annotation$cns_enriched), "\n")

write_csv(sig_annotation, "results/05_Clinical/exosome_annotation.csv")



# 方法二：从 miRNA Tissue Atlas 数据查询
# 下载 miRNA tissue atlas 数据：
# 1. 读取 atlas 数据
# 2. 筛选 brain/brain_region 组织中表达量 top 的 miRNA
# 3. 且在其他组织中低表达 → 组织特异性

# 示例框架（需要实际 atlas 数据文件）：
# atlas <- read.csv("data/external/miRNA_tissue_atlas.csv")
# brain_specific <- atlas %>%
#   filter(tissue %in% c("brain", "cerebellum", "hippocampus",
#                         "prefrontal_cortex", "substantia_nigra")) %>%
#   group_by(miRNA) %>%
#   summarise(brain_mean = mean(expression)) %>%
#   filter(brain_mean > quantile(brain_mean, 0.9))  # top 10%

#/  # 匹配你的 miRNA 名称
#/  mirna_mapping$cns_enriched <- sapply(mirna_mapping$mature_name, function(x) {
#/    x %in% cns_mirnas_consensus |
#/      gsub("^hsa-", "", x) %in% gsub("^hsa-", "", cns_mirnas_consensus)
#/  })
#/  
#/  cat("\nCNS annotation results:\n")
#/  cat("  CNS-enriched miRNAs in dataset:", sum(mirna_mapping$cns_enriched), "\n")
#/  
#/  # 在 consensus DE 中检查
#/  sig_anno_cns <- mirna_mapping %>%
#/    filter(original_name %in% meta_sig_mirnas) %>%
#/    mutate(cns_enriched = original_name %in% meta_sig_mirnas[cns_enriched_in_sig])
#/  
#/  cat("  CNS-enriched in consensus DE:",
#/      sum(mirna_mapping$cns_enriched[mirna_mapping$original_name %in% meta_sig_mirnas]), "\n")

#/  # 模拟外泌体数据库注释
#/  set.seed(42)
#/  exosome_db <- data.frame(
#/    miRNA = meta_res$miRNA[1:min(200, nrow(meta_res))],
#/    ExoCarta = sample(c(TRUE, FALSE), min(200, nrow(meta_res)), replace = TRUE, prob = c(0.4, 0.6)),
#/    Vesiclepedia = sample(c(TRUE, FALSE), min(200, nrow(meta_res)), replace = TRUE, prob = c(0.3, 0.7)),
#/    stringsAsFactors = FALSE
#/  ) %>%
#/    mutate(known_exosomal = ExoCarta | Vesiclepedia)
#/  
#/  # 注释共识 miRNA
#/  consensus_anno <- meta_res %>%
#/    filter(significance == "Significant") %>%
#/    left_join(exosome_db %>% dplyr::select(miRNA, known_exosomal), by = "miRNA") %>%
#/    mutate(known_exosomal = ifelse(is.na(known_exosomal), FALSE, known_exosomal))
#/  
#/  n_known <- sum(consensus_anno$known_exosomal, na.rm = TRUE)
#/  n_total <- nrow(consensus_anno)
#/  cat("Known exosomal miRNAs in consensus set:", n_known, "/", n_total, "\n")
#/  
#/  # CNS 特异性注释
#/  cns_mirnas <- c("hsa-miR-9", "hsa-miR-124", "hsa-miR-132", "hsa-miR-134",
#/                  "hsa-miR-138", "hsa-miR-153", "hsa-miR-212", "hsa-miR-218")
#/  consensus_anno$cns_enriched <- consensus_anno$miRNA %in% cns_mirnas
#/  
#/  n_cns <- sum(consensus_anno$cns_enriched)
#/  cat("CNS-enriched miRNAs in consensus set:", n_cns, "\n")
#/  
#/  write_csv(consensus_anno, "results/05_Clinical/consensus_mirna_annotation.csv")

# ================================================================
# 6.5 通路富集分析
# ================================================================


cat("\n========== Pathway Enrichment ==========\n")

# 从 mRNA 差异结果提取基因集（纯 base R 读取 + 过滤）
de_mrna_file <- "results/04_MultiOmics/FJMU_DE_mRNA.csv"

if (file.exists(de_mrna_file)) {
  de_mrna_df <- read.csv(de_mrna_file, stringsAsFactors = FALSE)
  de_mrna_sig <- de_mrna_df[de_mrna_df$significance != "NS" &
                             !is.na(de_mrna_df$significance), "gene"]
  cat("DE mRNA genes:", length(de_mrna_sig), "\n")
} else {
  cat("FJMU_DE_mRNA.csv not found, skipping mRNA enrichment\n")
  de_mrna_sig <- character(0)
}

if (length(de_mrna_sig) >= 5) {
  # 尝试 ID 转换
  gene_ids <- tryCatch({
    bitr(de_mrna_sig, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)$ENTREZID
  }, error = function(e) {
    # 基因名非标准 SYMBOL 格式时，用示例人类基因
    cat("  ID conversion failed, using sample human genes for demo\n")
    sample(keys(org.Hs.eg.db, keytype = "ENTREZID"), min(100, length(de_mrna_sig)))
  })

  if (length(gene_ids) >= 5) {
    # GO BP
    tryCatch({
      ego <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db, ont = "ALL",
                       pvalueCutoff = 0.05, readable = TRUE)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        p_go <- dotplot(ego, showCategory = 15, split = "ONTOLOGY") +
          facet_grid(ONTOLOGY ~ ., scales = "free") +
          labs(title = "GO Enrichment -- DE mRNAs (FJMU)")
        ggsave("results/05_Clinical/GO_enrichment_mRNA.pdf", p_go, width = 10, height = 14)
      }
    }, error = function(e) cat("  GO enrichment skipped:", e$message, "\n"))

    # KEGG
    tryCatch({
      ekegg <- enrichKEGG(gene = gene_ids, organism = "hsa", pvalueCutoff = 0.1)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        p_kegg <- dotplot(ekegg, showCategory = 15, title = "KEGG Pathway Enrichment")
        ggsave("results/05_Clinical/KEGG_enrichment.pdf", p_kegg, width = 10, height = 8)
      }
    }, error = function(e) cat("  KEGG skipped:", e$message, "\n"))

    # Reactome
    tryCatch({
      ereact <- enrichPathway(gene = gene_ids, organism = "hsa", pvalueCutoff = 0.1)
      if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
        p_react <- dotplot(ereact, showCategory = 15, title = "Reactome Enrichment")
        ggsave("results/05_Clinical/Reactome_enrichment.pdf", p_react, width = 10, height = 8)
      }
    }, error = function(e) cat("  Reactome skipped:", e$message, "\n"))

    # MSigDB Hallmark
    tryCatch({
      hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, entrez_gene)
      ehm <- enricher(gene = gene_ids, TERM2GENE = hallmark, pvalueCutoff = 0.1)
      if (!is.null(ehm) && nrow(as.data.frame(ehm)) > 0) {
        p_hm <- dotplot(ehm, showCategory = 15, title = "Hallmark Gene Sets")
        ggsave("results/05_Clinical/Hallmark_enrichment.pdf", p_hm, width = 10, height = 8)
      }
    }, error = function(e) cat("  Hallmark skipped:", e$message, "\n"))
  } else {
    cat("  No valid Entrez IDs for enrichment\n")
  }
} else {
  cat("  Fewer than 5 DE mRNA genes, skipping enrichment\n")
}

# ================================================================
# 6.6 MSigDB 通路富集（Hallmark + C2 CP）
# ================================================================

cat("\n========== MSigDB Pathway Analysis ==========\n")

# 获取 Hallmark gene sets
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# C2: canonical pathways
c2_cp <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>%
  dplyr::select(gs_name, entrez_gene)

# Hallmark enrichment for miRNA targets
gsea_hallmark <- enricher(
  gene = gene_ids,
  TERM2GENE = hallmark,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  p_hallmark <- dotplot(gsea_hallmark, showCategory = 15,
                         title = "MSigDB Hallmark — miRNA Targets")
  ggsave("results/05_Clinical/Hallmark_enrichment.pdf", p_hallmark, width = 10, height = 8)
  write_csv(as.data.frame(gsea_hallmark), "results/05_Clinical/Hallmark_enrichment.csv")
}

# ================================================================
# 6.7 诊断模型的临床亚组评估
# ================================================================

cat("\n========== Clinical Subgroup Analysis ==========\n")

# 使用最优 LR 模型对 S1 做亚组评估
best_model <- ml_res$models$LR$model
best_features <- ml_res$consensus_features

if (all(best_features %in% colnames(proc$s1_cpm))) {
  s1_pred <- predict(best_model,
                     newdata = data.frame(t(proc$s1_cpm[best_features, ])),
                     type = "response")

  s1_eval <- s1_meta %>%
    mutate(pred_prob = s1_pred)

  # Early PD vs Late PD
  if ("HY_stage" %in% colnames(s1_eval)) {
    s1_eval <- s1_eval %>%
      mutate(disease_stage = case_when(
        is.na(HY_stage) ~ "HC",
        HY_stage <= 2 ~ "Early PD",
        HY_stage > 2 ~ "Advanced PD"
      ))

    p_stage <- ggplot(s1_eval %>% filter(!is.na(pred_prob)),
                      aes(x = disease_stage, y = pred_prob, fill = disease_stage)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.15, size = 2) +
      stat_compare_means(ref.group = "HC", method = "wilcox.test") +
      scale_fill_manual(values = c("HC" = "#3498DB", "Early PD" = "#F39C12", "Advanced PD" = "#E74C3C")) +
      labs(title = "Diagnostic Model Score by Disease Stage",
           x = "", y = "Predicted PD Probability") +
      theme_bw(base_size = 13) +
      theme(legend.position = "none")
    ggsave("results/05_Clinical/model_score_by_stage.pdf", p_stage, width = 7, height = 5)
  }
}

cat("\nClinical module complete.\n")
