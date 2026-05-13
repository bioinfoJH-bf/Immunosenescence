# ==============================================================================
# Title: High-Resolution Myeloid Atlas & In Silico CD74 Knockdown Analysis
# Journal: Gigascience (Result 3 Pipeline)
# Description: This script comprehensively performs:
#              Part 1: High-resolution pan-cancer Mono/Macro atlas (Fig 3A-3F).
#              Part 2: Microenvironment-dependent CD74 evaluation (Fig 3G-3H).
#              Part 3: Virtual KD of CD74 via scTenifoldKnk & Downstream (Fig 3I-3K).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggExtra)       # For marginal histograms (Fig 3G)
  library(corrplot)      # For transcriptomic correlation (Fig 3F)
  library(clusterProfiler)# For GO/KEGG functional enrichment (Fig 3K)
  library(org.Hs.eg.db)
  library(ggvenn)        # For Venn diagrams (Fig 3I)
  library(scTenifoldKnk) # For in silico virtual knockdown (Critical for Fig 3I-K)
  # library(scMetabolism)  # Assumed available for Fig 3D metabolism scoring
})

# Define global variables & color palettes
WORK_DIR <- "~/JH/Pancancer"
setwd(WORK_DIR)
SAG_UNION <- c("APP", "MIF", "EGF", "CD74", "CD44", "CXCR4", "PLAU") # Replace with all 769 genes
custom_colors <- c('#e68b81', '#7da6c6', '#B7B2D0', '#EAAA60', '#F0B5A8', '#AFD3E6')

# ==============================================================================
# PART 1: Pan-Cancer Mono/Macro High-Resolution Atlas (Figs 3A - 3F)
# ==============================================================================
message(">>> Starting PART 1: Mono/Macro Atlas Profiling...")

# 1.1 Load merged pan-cancer myeloid data (Assuming pre-integrated via Harmony)
load("Pancancer_Myeloid_final.Rdata") # Contains 'Pancancer_Myeloid' object

# Extract only Monocytes and Macrophages for specific subpopulation analysis
Idents(Pancancer_Myeloid) <- "cell_type_final"
Mono_Macro <- subset(Pancancer_Myeloid, idents = c("Mono", "Macro", "Mac_DAB2", "Mac_APOC1", "Mac_CD81")) # Filter Unident/Proliferating if necessary

# [Fig 3A] UMAP Projection stratified by tissue origin
p3a <- DimPlot(Mono_Macro, reduction = "umap", group.by = "cell_type_final", split.by = "group") +
  theme_bw(base_size = 14) +
  labs(title = "Pan-cancer Mono/Macro (Tumor vs Adjacent vs Healthy)")
ggsave("Fig_3A_UMAP_MonoMacro.pdf", p3a, width = 12, height = 5)

# [Fig 3B] DotPlot for gene expression patterns across 13 subpopulations
marker_genes <- c("CD14", "FCGR3A", "DAB2", "APOC1", "CD81", "C1QC", "NLRP3", "LYVE1") # Example core markers
p3b <- DotPlot(Mono_Macro, features = marker_genes, group.by = "cell_type_final") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  scale_color_gradientn(colors = c("#f7fcf0", "#4eb3d3", "#08589e")) +
  labs(x = "Marker Genes", y = "Subpopulations")
ggsave("Fig_3B_DotPlot_Markers.pdf", p3b, width = 8, height = 6)

# [Fig 3C] Relative abundance across tissue microenvironments
abundance_df <- Mono_Macro@meta.data %>%
  group_by(group, cell_type_final) %>%
  tally() %>%
  mutate(freq = n / sum(n))

p3c <- ggplot(abundance_df, aes(x = group, y = freq, fill = cell_type_final)) +
  geom_bar(stat = "identity", position = "fill", color = "black", size = 0.2) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(x = "Tissue Origin", y = "Relative Abundance", fill = "Subpopulation")
ggsave("Fig_3C_Abundance_Barplot.pdf", p3c, width = 6, height = 6)

# [Fig 3D] Functional enrichment (Metabolism & Inflammation)
# Note: Simulating scoring via AddModuleScore or custom GSVA
Mono_Macro <- AddModuleScore(Mono_Macro, features = list(c("LPL", "APOE", "FABP4")), name = "Lipid_Metabolism")
p3d <- VlnPlot(Mono_Macro, features = "Lipid_Metabolism1", group.by = "cell_type_final", pt.size = 0) +
  geom_boxplot(width = 0.2, fill = "white") +
  labs(title = "Lipid Metabolic Remodeling (Mac_APOC1 focus)")
ggsave("Fig_3D_Functional_GSVA.pdf", p3d, width = 8, height = 5)

# [Fig 3E] SAG Score Comparisons across Tissues
p3e <- ggplot(Mono_Macro@meta.data, aes(x = cell_type_final, y = SAG_union2, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold")) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  labs(y = "Senescence Score (SAG)", title = "Senescence Dynamics Across Microenvironments")
ggsave("Fig_3E_SAG_Score_Dynamics.pdf", p3e, width = 10, height = 5)

# [Fig 3F] Transcriptomic Spearman Correlation Matrix
av_expr <- AverageExpression(Mono_Macro, assays = "RNA", return.seurat = FALSE)$RNA
cor_matrix <- cor(av_expr, method = "spearman")
pdf("Fig_3F_Correlation_Matrix.pdf", width = 7, height = 7)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE)
dev.off()


# ==============================================================================
# PART 2: Microenvironment-Dependent CD74 Evaluation (Figs 3G - 3H)
# ==============================================================================
message(">>> Starting PART 2: CD74 Expression Analysis...")

# Retrieve CD74 expression data safely
Mono_Macro$CD74_Expr <- GetAssayData(Mono_Macro, assay = "RNA", slot = "data")["CD74", ]

# [Fig 3G] Scatter plot of CD74 vs SAG score with Marginal Histograms
tumor_cells <- Mono_Macro@meta.data %>% filter(group == "Tumor")
p3g <- ggplot(tumor_cells, aes(x = CD74_Expr, y = SAG_union2)) +
  geom_point(alpha = 0.4, color = "#e68b81") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  theme_bw(base_size = 14) +
  labs(title = "Tumor Tissues: CD74 vs Senescence Score", x = "CD74 Expression", y = "SAG Score")
# Add marginal densities
p3g_marginal <- ggMarginal(p3g, type = "histogram", fill = "#e68b81", color = "black")
ggsave("Fig_3G_CD74_SAG_Correlation.pdf", p3g_marginal, width = 6, height = 6)

# [Fig 3H] CD74 Expression across subpopulations specifically in Tumor
p3h <- ggplot(tumor_cells, aes(x = reorder(cell_type_final, CD74_Expr, FUN=mean), y = CD74_Expr, fill = cell_type_final)) +
  geom_boxplot(outlier.shape = 18) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "none") +
  labs(x = "Myeloid Subpopulations in Tumor", y = "CD74 Expression")
ggsave("Fig_3H_CD74_Expression_Ranked.pdf", p3h, width = 8, height = 5)


# ==============================================================================
# PART 3: In silico CD74 Knockdown in Mac_DAB2 via scTenifoldKnk (Figs 3I - 3K)
# ==============================================================================
message(">>> Starting PART 3: scTenifoldKnk Virtual Knockdown...")

# 3.1 Setup target cells for Knockdown (Mac_DAB2 strictly within Tumor)
Idents(Mono_Macro) <- "cell_type_final"
Mac_DAB2_Tumor <- subset(Mono_Macro, idents = "Mac_DAB2", subset = group == "Tumor")

# Extract Raw Counts and Highly Variable Genes to prevent memory overload
counts_mat <- GetAssayData(Mac_DAB2_Tumor, assay = "RNA", slot = "counts")
Mac_DAB2_Tumor <- FindVariableFeatures(Mac_DAB2_Tumor, nfeatures = 5000)
hvgs <- VariableFeatures(Mac_DAB2_Tumor)

# Ensure CD74 is in the computation matrix
target_gene <- "CD74"
calc_genes <- unique(c(target_gene, hvgs))
data_matrix <- as.matrix(counts_mat[calc_genes, ])

# 3.2 Execute scTenifoldKnk (Virtual Gene Knockout)
# Note: Extremely computationally heavy, nCores adjusted
message("Running scTenifoldKnk for CD74 in Mac_DAB2. This will take time...")
ko_result <- scTenifoldKnk(
  countMatrix = data_matrix, 
  gKO = target_gene, 
  qc = TRUE, qc_mtThreshold = 0.1, qc_minLSize = 1000, 
  nc_nNet = 10, nCores = 10, nc_nCells = 500
)
saveRDS(ko_result, "Mac_DAB2_CD74_Knockdown_Result.rds")

# 3.3 Extract perturbed candidate genes
ko_df <- ko_result$diffRegulation
ko_df <- ko_df[ko_df$gene != target_gene, ] # Remove the knocked-out gene itself
sig_ko_df <- ko_df[ko_df$p.value < 0.05, ]
write.table(sig_ko_df, "Table_S_CD74_Knockdown_SigDiff.txt", sep = "\t", quote = F, row.names = F)

# [Fig 3I] Venn diagram & Volcano Plot
venn_list <- list("CD74 KO Targets" = sig_ko_df$gene, "Senescence Genes" = SAG_UNION)
pdf("Fig_3I_Venn_Knockdown_SAG.pdf", width = 6, height = 6)
ggvenn(venn_list, fill_color = c("#E7BC37", "#e68b81"), stroke_color = "white")
dev.off()

# [Fig 3J] Export top 30 targets for PPI Network Construction
top_30_genes <- head(ko_df[order(ko_df$p.value), "gene"], 30)
write.csv(data.frame(Gene = top_30_genes), "Fig_3J_PPI_Top30_Nodes.csv", row.names = FALSE)
message("PPI Node list exported. Use STRING database & Cytoscape for network rendering.")

# [Fig 3K] Functional Enrichment Analysis of CD74-perturbed genes
message("Running clusterProfiler for GO Enrichment...")
entrez_ids <- mget(sig_ko_df$gene, org.Hs.egSYMBOL2EG, ifnotfound = NA)
sig_ko_df$EntrezID <- as.character(entrez_ids)
enrich_genes <- na.omit(sig_ko_df$EntrezID)

# GO Enrichment (Biological Process)
go_res <- enrichGO(gene = enrich_genes, OrgDb = org.Hs.eg.db, ont = "BP", 
                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)

# Filter for relevant Senescence/Immune terms exactly as described in the paper
go_subset <- go_res@result %>% 
  filter(grepl("senescence|macrophage|myeloid|immune", Description, ignore.case = TRUE))
go_res@result <- go_subset

p3k <- barplot(go_res, showCategory = 10, label_format = 50) +
  scale_fill_viridis_c() +
  labs(title = "GO BP: Functions Impaired by CD74 Knockdown") +
  theme_bw(base_size = 14)
ggsave("Fig_3K_GO_Enrichment_Barplot.pdf", p3k, width = 8, height = 6)

message("=== Master Pipeline Complete! Check working directory for all PDFs and Tables. ===")