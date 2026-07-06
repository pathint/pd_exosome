# ============================================================
# 08b_Protein_Enrichment.R（更新版）
# 适配 08_Protein_Analysis.R 新输出
# 蛋白 ID 为 Gene Symbol → keytype = "SYMBOL"
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

select <- dplyr::select

dir.create("results/08_Protein/Enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/08_Protein/Integration", showWarnings = FALSE)

# ================================================================
# 1. 加载输入文件
# ================================================================

cat("Loading input files...\n")

consensus_prot <- read_csv("results/08_Protein/Enrichment_input_consensus_proteins.csv",
                            show_col_types = FALSE) %>% pull(Protein)

de_protein <- read_csv("results/08_Protein/DE_proteins_PD_vs_HC.csv",
                        show_col_types = FALSE)

de_j1 <- read_csv("results/08_Protein/DE_JPST003013_independent.csv",
                    show_col_types = FALSE)

de_j2 <- read_csv("results/08_Protein/DE_JPST003026_independent.csv",
                    show_col_types = FALSE)

hub_summary <- read_csv("results/08_Protein/Enrichment_input_WGCNA_hubs.csv",
                          show_col_types = FALSE)

module_gene_list <- read_csv("results/08_Protein/Enrichment_input_WGCNA_all_genes.csv",
                              show_col_types = FALSE)

# J1 WGCNA
j1_hub_summary <- read_csv("results/08_Protein/WGCNA_J1/hub_genes_all_modules.csv",
                             show_col_types = FALSE)

j1_module_genes <- read_csv("results/08_Protein/WGCNA_J1/all_module_genes.csv",
                              show_col_types = FALSE)

# ================================================================
# 2. Gene Symbol → ENTREZID 映射
# ================================================================

cat("\n========== ID Mapping (Gene Symbol → ENTREZID) ==========\n")

all_symbols <- unique(c(consensus_prot,
                         de_protein$Protein,
                         de_j1$Protein,
                         de_j2$Protein,
                         hub_summary$Protein,
                         j1_hub_summary$Protein))

# 清洗：去除可能的 UniProt 残留格式
all_symbols_clean <- trimws(all_symbols)

entrez_map <- mapIds(org.Hs.eg.db,
                     keys = all_symbols_clean,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

sym_map <- mapIds(org.Hs.eg.db,
                  keys = all_symbols_clean,
                  column = "SYMBOL",
                  keytype = "SYMBOL",
                  multiVals = "first")

mapping_df <- data.frame(
  Protein  = all_symbols,
  Gene_Symbol = unname(sym_map[all_symbols]),
  ENTREZID = unname(entrez_map[all_symbols]),
  stringsAsFactors = FALSE
) %>% distinct(Protein, .keep_all = TRUE)

# 映射失败的尝试小写/大小写变体
failed <- mapping_df %>% filter(is.na(ENTREZID)) %>% pull(Protein)
if (length(failed) > 0) {
  cat("Retrying", length(failed), "unmapped proteins...\n")
  retry_entrez <- mapIds(org.Hs.eg.db, keys = toupper(failed),
                          column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  for (p in failed) {
    if (!is.na(retry_entrez[p])) {
      mapping_df$ENTREZID[mapping_df$Protein == p] <- retry_entrez[p]
    }
  }
}

write_csv(mapping_df, "results/08_Protein/Enrichment/protein_mapping.csv")
cat("Mapped:", sum(!is.na(mapping_df$ENTREZID)), "/", nrow(mapping_df), "\n")

# ================================================================
# 3. 辅助函数
# ================================================================

run_go_kegg <- function(gene_list, label, p_cut = 0.05) {
  entrez <- mapping_df %>%
    dplyr::filter(Protein %in% gene_list, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>% unique()

  cat(sprintf("  %s: %d genes with ENTREZID\n", label, length(entrez)))
  if (length(entrez) < 3) return(list(go = NULL, kegg = NULL))

  ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = p_cut,
                   qvalueCutoff = 0.2, readable = TRUE)

  entrez_clean <- as.character(as.integer(entrez))
  entrez_clean <- entrez_clean[!is.na(entrez_clean)]

  ekegg <- NULL
  tryCatch({
    ekegg <- enrichKEGG(gene = entrez_clean, organism = "hsa",
                         keyType = "ncbi-geneid", pvalueCutoff = p_cut)
    if (nrow(as.data.frame(ekegg)) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }
  }, error = function(e) cat(sprintf("  KEGG failed for %s: %s\n", label, e$message)))

  return(list(go = ego, kegg = ekegg))
}

save_and_plot <- function(result, label) {
  if (!is.null(result$go) && nrow(as.data.frame(result$go)) > 0) {
    write_csv(as.data.frame(result$go),
              paste0("results/08_Protein/Enrichment/", label, "_GO_BP.csv"))
    p <- dotplot(result$go, showCategory = 15) +
      ggtitle(paste0(label, " GO:BP")) + theme_bw(base_size = 11)
    ggsave(paste0("results/08_Protein/Enrichment/", label, "_GO_BP_dotplot.pdf"),
           p, width = 9, height = 7)
    cat(sprintf("  %s GO:BP: %d terms\n", label, nrow(as.data.frame(result$go))))
  }
  if (!is.null(result$kegg) && nrow(as.data.frame(result$kegg)) > 0) {
    write_csv(as.data.frame(result$kegg),
              paste0("results/08_Protein/Enrichment/", label, "_KEGG.csv"))
    p <- dotplot(result$kegg, showCategory = 15) +
      ggtitle(paste0(label, " KEGG")) + theme_bw(base_size = 11)
    ggsave(paste0("results/08_Protein/Enrichment/", label, "_KEGG_dotplot.pdf"),
           p, width = 9, height = 7)
    cat(sprintf("  %s KEGG: %d terms\n", label, nrow(as.data.frame(result$kegg))))
  }
}

# ================================================================
# 4. ORA 富集分析
# ================================================================

cat("\n========== ORA Enrichment ==========\n")

# --- 4.1 J1 DE 蛋白（主要信号，40 个 DE）---
cat("\n--- J1 DE Proteins ---\n")
j1_de_up   <- de_j1 %>% dplyr::filter(direction == "Up") %>% pull(Protein)
j1_de_down <- de_j1 %>% dplyr::filter(direction == "Down") %>% pull(Protein)
j1_de_all  <- de_j1 %>% dplyr::filter(significance == "DE") %>% pull(Protein)

save_and_plot(run_go_kegg(j1_de_all, "J1_DE_all"), "J1_DE_all")
save_and_plot(run_go_kegg(j1_de_up,   "J1_DE_up"),  "J1_DE_up")
save_and_plot(run_go_kegg(j1_de_down, "J1_DE_down"), "J1_DE_down")

# --- 4.2 J2 DE 蛋白（5 个 DE）---
cat("\n--- J2 DE Proteins ---\n")
j2_de_all <- de_j2 %>% dplyr::filter(significance == "DE") %>% pull(Protein)
save_and_plot(run_go_kegg(j2_de_all, "J2_DE_all"), "J2_DE_all")

# --- 4.3 Combined DE 蛋白（仅 2 个，可能不够做 ORA）---
cat("\n--- Combined DE Proteins ---\n")
comb_de_all <- de_protein %>% dplyr::filter(significance == "DE") %>% pull(Protein)
if (length(comb_de_all) >= 3) {
  save_and_plot(run_go_kegg(comb_de_all, "Combined_DE_all"), "Combined_DE_all")
} else {
  cat(sprintf("  Only %d combined DE proteins, ORA not feasible\n", length(comb_de_all)))
}

# --- 4.4 J1∩J2 交集 DE ---
cat("\n--- J1 ∩ J2 DE ---\n")
j1j2_up   <- intersect(
  de_j1 %>% dplyr::filter(direction == "Up") %>% pull(Protein),
  de_j2 %>% dplyr::filter(direction == "Up") %>% pull(Protein))
j1j2_down <- intersect(
  de_j1 %>% dplyr::filter(direction == "Down") %>% pull(Protein),
  de_j2 %>% dplyr::filter(direction == "Down") %>% pull(Protein))
j1j2_all  <- c(j1j2_up, j1j2_down)

cat(sprintf("  J1∩J2 Up: %d, Down: %d\n", length(j1j2_up), length(j1j2_down)))

if (length(j1j2_all) >= 3) {
  save_and_plot(run_go_kegg(j1j2_all, "J1_J2_intersect_DE"), "J1_J2_intersect_DE")
}

# --- 4.5 Consensus Signature（63 个蛋白）---
cat("\n--- Consensus Signature ---\n")
save_and_plot(run_go_kegg(consensus_prot, "Consensus_Signature"), "Consensus_Signature")

# --- 4.6 WGCNA Combined 显著模块 ---
cat("\n--- WGCNA Combined Modules ---\n")
if (length(hub_summary$Module) > 0) {
  unique_mods <- unique(hub_summary$Module)
  for (mod in unique_mods) {
    mod_genes <- hub_summary %>% dplyr::filter(Module == mod) %>% pull(Protein)
    save_and_plot(run_go_kegg(mod_genes, paste0("WGCNA_combined_", mod)),
                  paste0("WGCNA_combined_", mod))
  }
}

# --- 4.7 WGCNA J1 所有模块（逐模块富集）---
cat("\n--- WGCNA J1 Modules ---\n")
j1_unique_mods <- unique(j1_module_genes$Module)

for (mod in j1_unique_mods) {
  mod_genes <- j1_module_genes %>%
    dplyr::filter(Module == mod) %>% pull(Protein)
  if (length(mod_genes) < 10) next

  save_and_plot(run_go_kegg(mod_genes, paste0("WGCNA_J1_", mod)),
                paste0("WGCNA_J1_", mod))
}

# --- 4.8 WGCNA J1 Hub Genes ---
cat("\n--- WGCNA J1 Hub Genes ---\n")
j1_hub_all <- unique(j1_hub_summary$Protein)
if (length(j1_hub_all) >= 3) {
  save_and_plot(run_go_kegg(j1_hub_all, "WGCNA_J1_hubs"), "WGCNA_J1_hubs")
}

# ================================================================
# 5. GSEA（ranked by logFC，适合信号弱的场景）
# ================================================================

cat("\n========== GSEA Analysis ==========\n")

run_gsea <- function(de_df, label) {
  de_mapped <- de_df %>%
    dplyr::filter(!is.na(Protein)) %>%
    left_join(mapping_df %>% dplyr::select(Protein, ENTREZID), by = "Protein") %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  if (nrow(de_mapped) < 10) {
    cat(sprintf("  %s: only %d mapped genes, GSEA skipped\n", label, nrow(de_mapped)))
    return(NULL)
  }

  gene_list <- setNames(de_mapped$logFC, de_mapped$ENTREZID)
  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", ont = "BP",
                    pvalueCutoff = 0.05, verbose = FALSE)

  gsea_kegg <- NULL
  tryCatch({
    gsea_kegg <- gseKEGG(geneList = gene_list, organism = "hsa",
                           keyType = "ncbi-geneid", pvalueCutoff = 0.05,
                           verbose = FALSE)
  }, error = function(e) {})

  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    write_csv(as.data.frame(gsea_go),
              paste0("results/08_Protein/Enrichment/GSEA_", label, "_GO_BP.csv"))
    p <- dotplot(gsea_go, showCategory = 15) +
      ggtitle(paste0("GSEA GO:BP - ", label)) + theme_bw(base_size = 11)
    ggsave(paste0("results/08_Protein/Enrichment/GSEA_", label, "_GO_BP_dotplot.pdf"),
           p, width = 9, height = 7)
    cat(sprintf("  GSEA %s GO:BP: %d terms\n", label, nrow(as.data.frame(gsea_go))))
  }
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    write_csv(as.data.frame(gsea_kegg),
              paste0("results/08_Protein/Enrichment/GSEA_", label, "_KEGG.csv"))
    p <- dotplot(gsea_kegg, showCategory = 15) +
      ggtitle(paste0("GSEA KEGG - ", label)) + theme_bw(base_size = 11)
    ggsave(paste0("results/08_Protein/Enrichment/GSEA_", label, "_KEGG_dotplot.pdf"),
           p, width = 9, height = 7)
    cat(sprintf("  GSEA %s KEGG: %d terms\n", label, nrow(as.data.frame(gsea_kegg))))
  }

  return(list(go = gsea_go, kegg = gsea_kegg))
}

# GSEA on J1（全蛋白排序，不需要 DE 阈值）
cat("\n--- GSEA: J1 (all proteins ranked) ---\n")
gsea_j1 <- run_gsea(de_j1, "J1")

# GSEA on J2
cat("\n--- GSEA: J2 (all proteins ranked) ---\n")
gsea_j2 <- run_gsea(de_j2, "J2")

# GSEA on Combined
cat("\n--- GSEA: Combined (all proteins ranked) ---\n")
gsea_comb <- run_gsea(de_protein, "Combined")

# ================================================================
# 6. 蛋白-miRNA 通路层面整合
# ================================================================

cat("\n========== Protein-miRNA Pathway Integration ==========\n")

mirna_go_path   <- "results/04_Pathway/GO_BP_enrichment.csv"
mirna_kegg_path <- "results/04_Pathway/KEGG_enrichment.csv"

# 读取最佳蛋白 GO 结果（优先 J1，信号最强）
best_go <- NULL
best_go_label <- ""
for (lbl in c("J1_DE_all", "Consensus_Signature")) {
  fpath <- paste0("results/08_Protein/Enrichment/", lbl, "_GO_BP.csv")
  if (file.exists(fpath)) {
    tmp <- read_csv(fpath, show_col_types = FALSE)
    if (nrow(tmp) > 0) {
      best_go <- tmp
      best_go_label <- lbl
      break
    }
  }
}

shared_go <- data.frame()
if (!is.null(best_go) && file.exists(mirna_go_path)) {
  mirna_go <- read_csv(mirna_go_path, show_col_types = FALSE)
  shared_ids <- intersect(best_go$ID, mirna_go$ID)
  cat(sprintf("Shared GO:BP terms (protein [%s] & miRNA): %d\n",
              best_go_label, length(shared_ids)))
  if (length(shared_ids) > 0) {
    shared_go <- bind_rows(
      best_go %>% dplyr::filter(ID %in% shared_ids) %>% mutate(Source = "Protein"),
      mirna_go %>% dplyr::filter(ID %in% shared_ids) %>% mutate(Source = "miRNA"))
    write_csv(shared_go, "results/08_Protein/Integration/shared_GO_terms.csv")
  }
}

best_kegg <- NULL
best_kegg_label <- ""
for (lbl in c("J1_DE_all", "Consensus_Signature")) {
  fpath <- paste0("results/08_Protein/Enrichment/", lbl, "_KEGG.csv")
  if (file.exists(fpath)) {
    tmp <- read_csv(fpath, show_col_types = FALSE)
    if (nrow(tmp) > 0) {
      best_kegg <- tmp
      best_kegg_label <- lbl
      break
    }
  }
}

shared_kegg <- data.frame()
if (!is.null(best_kegg) && file.exists(mirna_kegg_path)) {
  mirna_kegg <- read_csv(mirna_kegg_path, show_col_types = FALSE)
  shared_ids <- intersect(best_kegg$ID, mirna_kegg$ID)
  cat(sprintf("Shared KEGG (protein [%s] & miRNA): %d\n",
              best_kegg_label, length(shared_ids)))
  if (length(shared_ids) > 0) {
    shared_kegg <- bind_rows(
      best_kegg %>% dplyr::filter(ID %in% shared_ids) %>% mutate(Source = "Protein"),
      mirna_kegg %>% dplyr::filter(ID %in% shared_ids) %>% mutate(Source = "miRNA"))
    write_csv(shared_kegg, "results/08_Protein/Integration/shared_KEGG_pathways.csv")
  }
}

# 共有通路可视化
if (nrow(shared_go) > 0) {
  p_shared <- ggplot(shared_go %>% distinct(ID, .keep_all = TRUE),
                      aes(x = reorder(Description, -log10(pvalue)),
                          y = -log10(pvalue), fill = Source)) +
    geom_col(position = "dodge", width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Protein" = "#E74C3C", "miRNA" = "#3498DB")) +
    labs(title = sprintf("Shared GO:BP (Protein [%s] & miRNA)", best_go_label),
         x = "", y = "-log10(p-value)") +
    theme_bw(base_size = 11)
  ggsave("results/08_Protein/Integration/shared_GO_barplot.pdf", p_shared,
         width = 10, height = max(5, nrow(shared_go %>% distinct(ID)) * 0.3))
}

# ================================================================
# 7. 汇总
# ================================================================

cat("\n========================================\n")
cat("Protein Enrichment Analysis Complete\n")
cat("========================================\n")

# 统计所有输出
enrich_files <- list.files("results/08_Protein/Enrichment", pattern = "\\.csv$", full.names = TRUE)
cat(sprintf("  Total enrichment result files: %d\n", length(enrich_files)))
cat(sprintf("  Shared GO:BP (protein & miRNA): %d\n", nrow(shared_go)))
cat(sprintf("  Shared KEGG (protein & miRNA): %d\n", nrow(shared_kegg)))
cat("\nAll results in: results/08_Protein/Enrichment/\n")
