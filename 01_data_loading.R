# ============================================================
# 01_data_loading.R —— 数据加载
# ============================================================

meta <- read.table("data/raw/meta_sample_info.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

rownames(meta) <- meta$ID
head(meta)
dim(meta)

meta2020 <- subset(meta, Dataset == "GSE269775")
dim(meta2020)

meta2021 <- subset(meta, Dataset == "GSE269776")
dim(meta2021)

metaFJMU <- subset(meta, Dataset == "FJMU")
dim(metaFJMU)

miRNA2021 <- read.table("data/raw/GSE269775_2020_miRNA_expression_matrix.csv", header = TRUE, row.names = 1, sep = ",")
colnames(miRNA2021) <- rownames(meta2021)
head(miRNA2021)
dim(miRNA2021)


miRNA2020 <- read.table("data/raw/GSE269776_2021_miRNA_expression_matrix.csv", header = TRUE, row.names = 1, sep = ",")
colnames(miRNA2020) <- rownames(meta2020)
head(miRNA2020)
dim(miRNA2020)

keep <- colSums(miRNA2020) > 5e4
miRNA2020 <- miRNA2020[,keep]
meta2020  <- meta2020[colnames(miRNA2020),]

excl <- c("HC152",
		  "HC158",
		  "HC159",
		  "HC160",
		  "HC161",
		  "HC163",
		  "HC164",
		  "HC165",
		  "HC167",
		  "HC168",
		  "HC169",
		  "HC172",
		  "HC170",
		  "PD143")
keep <- !(colnames(miRNA2020) %in% excl)

miRNA2020 <- miRNA2020[,keep]
meta2020  <- meta2020[colnames(miRNA2020),]

keep <- colSums(miRNA2021) > 5e4
miRNA2021 <- miRNA2021[,keep]
meta2021  <- meta2021[colnames(miRNA2021),]

miRNA_FJMU <- read.table("data/raw/FJMU_miRNA_counts.csv", header = TRUE, row.names = 1, sep = ",")
table(colnames(miRNA_FJMU) == rownames(metaFJMU))
head(miRNA_FJMU)
dim(miRNA_FJMU)

mRNA_FJMU <- read.table("data/raw/FJMU_mRNA_counts.csv", header = TRUE, row.names = 1, sep = ",")
table(colnames(mRNA_FJMU) %in% rownames(metaFJMU))
head(mRNA_FJMU)
dim(mRNA_FJMU)

lncRNA_FJMU <- read.table("data/raw/FJMU_lncRNA_counts.csv", header = TRUE, row.names = 1, sep = ",")
colnames(lncRNA_FJMU)
lncRNA_info <- lncRNA_FJMU[,c(1:6, 20:22)] 
lncRNA_FJMU <- lncRNA_FJMU[,c(7:19)]

table(colnames(lncRNA_FJMU) %in% rownames(metaFJMU))
head(lncRNA_FJMU)
dim(lncRNA_FJMU)

table(lncRNA_info[,9] > 0)

protein_FJMU <- read.table("data/raw/FJMU_protein_expr.csv", header = TRUE, row.names = 1, sep = ",")
dim(protein_FJMU)
table(colnames(protein_FJMU) %in% rownames(metaFJMU))
head(protein_FJMU)


proc <- list(
  # 公共数据集1
  g1_counts = miRNA2020,
  g1_meta   = meta2020,
  # 公共数据集2
  g2_counts = miRNA2021,
  g2_meta   = meta2021,
  # 自测数据
  s1_mirna_counts  = miRNA_FJMU,
  s1_mrna_counts   = mRNA_FJMU,
  s1_lncrna_counts = lncRNA_FJMU,
  s1_lncrna_info   = lncRNA_info,
  s1_protein       = protein_FJMU,
  s1_meta          = metaFJMU,
  # 完整meta
  meta_all = meta
)


# 检查各数据集 miRNA 名称交集
g1_mirnas <- rownames(proc$g1_counts)
g2_mirnas <- rownames(proc$g2_counts)
s1_mirnas <- rownames(proc$s1_mirna_counts)

common_all <- Reduce(intersect, list(g1_mirnas, g2_mirnas, s1_mirnas))
common_g1g2 <- intersect(g1_mirnas, g2_mirnas)

cat("G1 miRNAs:", length(g1_mirnas), "\n")
cat("G2 miRNAs:", length(g2_mirnas), "\n")
cat("S1 miRNAs:", length(s1_mirnas), "\n")
cat("G1 ∩ G2:",   length(common_g1g2), "\n")
cat("G1 ∩ G2 ∩ S1:", length(common_all), "\n")

# 检查 FJMU 多组学样本一致性
fjmu_mirna_samp <- colnames(proc$s1_mirna_counts)
fjmu_mrna_samp  <- colnames(proc$s1_mrna_counts)
fjmu_prot_samp  <- colnames(proc$s1_protein)
fjmu_lnc_samp   <- colnames(proc$s1_lncrna_counts)

cat("\nFJMU sample overlap:\n")
cat("  miRNA ∩ mRNA:", length(intersect(fjmu_mirna_samp, fjmu_mrna_samp)), "\n")
cat("  miRNA ∩ protein:", length(intersect(fjmu_mirna_samp, fjmu_prot_samp)), "\n")
cat("  miRNA ∩ lncRNA:", length(intersect(fjmu_mirna_samp, fjmu_lnc_samp)), "\n")
cat("  All 4:", length(Reduce(intersect, list(fjmu_mirna_samp, fjmu_mrna_samp,
                                                fjmu_prot_samp, fjmu_lnc_samp))), "\n")

# 检查临床变量可用性
pd_meta <- subset(proc$g1_meta, Group == "PD")
cat("\nClinical data availability (GSE269775 PD):\n")
cat("  Age:", sum(!is.na(pd_meta$Age)), "/", nrow(pd_meta), "\n")
cat("  Gender:", sum(!is.na(pd_meta$Gender)), "/", nrow(pd_meta), "\n")
cat("  Yahr:", sum(pd_meta$Yahr != "-" & !is.na(pd_meta$Yahr)), "/", nrow(pd_meta), "\n")
cat("  UPDRS-III:", sum(pd_meta$UPDRS.III != "-" & !is.na(pd_meta$UPDRS.III)), "/", nrow(pd_meta), "\n")
cat("  LEDD:", sum(pd_meta$LEDD > 0, na.rm = TRUE), "/", nrow(pd_meta), "\n")
cat("  MIBG early H/M:", sum(!is.na(pd_meta$early.H.M) & pd_meta$early.H.M != "-"), "/", nrow(pd_meta), "\n")
cat("  DaTscan Ave:", sum(!is.na(pd_meta$DaTscan.Ave) & pd_meta$DaTscan.Ave != "-"), "/", nrow(pd_meta), "\n")

# 保存
saveRDS(proc, "data/processed/proc_raw.rds")

common_all <- read.csv("common_mirnas.csv")
common_all <- common_all$x

common_g1g2 <- read.csv("g1g2_common.csv")
common_g1g2 <- common_g1g2$x

write.csv(metaFJMU, file = "FJMU_meta_filtered.csv")
write.csv(meta2020, file = "GSE269775_meta_filtered.csv")
write.csv(meta2021, file = "GSE269776_meta_filtered.csv")

table(rownames(meta2020) == colnames(miRNA2020))

write.csv(miRNA2020[common_g1g2,], file = "GSE269775_counts_g1g2.csv")
write.csv(miRNA2021[common_g1g2,], file = "GSE269776_counts_g1g2.csv")


write.csv(miRNA_FJMU[common_all,], file = "FJMU_miRNA_counts_common.csv")
write.csv(miRNA2020[common_all,], file = "GSE269775_counts_common.csv")
write.csv(miRNA2021[common_all,], file = "GSE269776_counts_common.csv")


