# ==============================================================================
# Title: Master Pipeline for Result 5 - MSAG.SIG Subtyping & Clinical Validation
# Journal: Gigascience
# Description: This script executes the translation of single-cell signatures 
#              to bulk TCGA cohorts. It covers Univariate Cox (Fig 5A), Sankey 
#              Enrichment (Fig 5B), Consensus Clustering (Fig 5C-D), PCA (Fig 5E),
#              Survival (Fig 5G), DEGs (Fig 5H-I), and TME analysis (Fig 5K-L).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Global Setup & Library Loading
# ------------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ConsensusClusterPlus)
  library(scatterplot3d)
  library(ComplexHeatmap)
  library(limma)
  library(ggpubr)
  library(RColorBrewer)
})

# Define global variables & Working Directories
WORK_DIR <- "~/JH/Pancancer_TCGA"
if(!dir.exists(WORK_DIR)) dir.create(WORK_DIR, recursive = TRUE)
setwd(WORK_DIR)

# The 59 initial candidates from single-cell analysis (49 key genes + 10 core TFs)
# Replace with actual gene list
CANDIDATE_59_GENES <- c("APP", "MIF", "CD74", "CD44", "ATF3", "NR1H3", "CXCR4") # Placeholder

# ------------------------------------------------------------------------------
# PART 1: Univariate Cox Regression & Enrichment (Fig 5A, 5B)
# ------------------------------------------------------------------------------
run_msag_identification <- function(expr_matrix, clin_data) {
  message(">>> [Part 1] Performing Univariate Cox Regression (Fig 5A)...")
  
  # Merge expression with survival data (assuming clin_data has 'OS.time' and 'OS')
  merged_data <- merge(clin_data, t(expr_matrix), by = "row.names")
  
  # Univariate Cox loop
  cox_results <- data.frame()
  for (gene in intersect(CANDIDATE_59_GENES, colnames(merged_data))) {
    formula <- as.formula(paste("Surv(OS.time, OS) ~", gene))
    fit <- coxph(formula, data = merged_data)
    summary_fit <- summary(fit)
    
    cox_results <- rbind(cox_results, data.frame(
      Gene = gene,
      HR = summary_fit$conf.int[1],
      Lower_95 = summary_fit$conf.int[3],
      Upper_95 = summary_fit$conf.int[4],
      P_value = summary_fit$waldtest[3]
    ))
  }
  
  # Filter significant genes (The 35 MSAG.SIG)
  msag_sig <- cox_results %>% filter(P_value < 0.05)
  write.csv(msag_sig, "Table_MSAG_SIG_35_Genes.csv", row.names = FALSE)
  
  # [Fig 5A] Forest Plot
  p_forest <- ggforest(coxph(Surv(OS.time, OS) ~ ., data = merged_data[, c("OS.time", "OS", msag_sig$Gene)]), data = merged_data)
  ggsave("Fig_5A_Cox_Forest.pdf", p_forest, width = 10, height = 12)
  
  # [Fig 5B] Functional Enrichment (Sankey-Bar Concept)
  message(">>> [Part 1] Functional Enrichment for 35 MSAG Genes (Fig 5B)...")
  entrez_ids <- mget(msag_sig$Gene, org.Hs.egSYMBOL2EG, ifnotfound = NA)
  kegg_res <- enrichKEGG(gene = na.omit(as.character(entrez_ids)), organism = "hsa", pvalueCutoff = 0.05)
  
  # Extract top pathways for plotting (simplified representation of your custom plot)
  p_bar <- barplot(kegg_res, showCategory = 10) + scale_fill_viridis_c() + labs(title="KEGG Pathways (PD-1, NF-kB, etc.)")
  ggsave("Fig_5B_KEGG_Bar.pdf", p_bar, width=8, height=6)
  
  return(msag_sig$Gene)
}

# ------------------------------------------------------------------------------
# PART 2: Consensus Clustering & PCA (Fig 5C - 5E, 5F)
# ------------------------------------------------------------------------------
run_subtyping_and_pca <- function(expr_matrix, msag_genes) {
  message(">>> [Part 2] Unsupervised Consensus Clustering (Fig 5C-D)...")
  
  # Subset matrix to MSAG.SIG genes only
  mat_sub <- expr_matrix[intersect(rownames(expr_matrix), msag_genes), ]
  mat_sub <- sweep(mat_sub, 1, apply(mat_sub, 1, median)) # Median centering
  
  # Run ConsensusClusterPlus
  results <- ConsensusClusterPlus(
    as.matrix(mat_sub), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
    title = "Consensus_Clustering", clusterAlg = "hc", distance = "pearson", plot = "pdf"
  )
  
  # Extract K=2 (Group 1 & Group 2)
  cluster_assignments <- results[[2]]$consensusClass
  subtypes_df <- data.frame(Sample = names(cluster_assignments), Subtype = paste0("Group", cluster_assignments))
  
  # [Fig 5E] 3D PCA
  pca_res <- prcomp(t(mat_sub), scale. = TRUE)
  colors <- ifelse(subtypes_df$Subtype == "Group1", "#7da6c6", "#D97373")
  
  pdf("Fig_5E_3D_PCA.pdf", width = 6, height = 6)
  s3d <- scatterplot3d(pca_res$x[,1:3], color = colors, pch = 16, main = "3D PCA (Group 1 vs Group 2)",
                       xlab = "PC1", ylab = "PC2", zlab = "PC3", angle = 45)
  legend("topright", legend = c("Group 1", "Group 2"), col = c("#7da6c6", "#D97373"), pch = 16)
  dev.off()
  
  # [Fig 5F] Overall MSAG Score Comparison
  msag_scores <- colMeans(mat_sub)
  score_df <- merge(subtypes_df, data.frame(Sample = names(msag_scores), Score = msag_scores), by = "Sample")
  
  p_score <- ggplot(score_df, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_boxplot() + scale_fill_manual(values = c("#7da6c6", "#D97373")) +
    stat_compare_means(method = "wilcox.test") + theme_bw() + labs(title = "Overall MSAG Score")
  ggsave("Fig_5F_MSAG_Score_Boxplot.pdf", p_score, width = 5, height = 5)
  
  return(subtypes_df)
}

# ------------------------------------------------------------------------------
# PART 3: Survival, DEGs & Clinical Characteristics (Fig 5G - 5J)
# ------------------------------------------------------------------------------
run_clinical_deg_evaluation <- function(expr_matrix, clin_data, subtypes_df) {
  message(">>> [Part 3] Survival & Clinical Feature Evaluation...")
  
  # Merge data
  full_clin <- merge(clin_data, subtypes_df, by.x = "row.names", by.y = "Sample")
  
  # [Fig 5G] Kaplan-Meier Survival Analysis
  fit_km <- survfit(Surv(OS.time, OS) ~ Subtype, data = full_clin)
  p_km <- ggsurvplot(fit_km, data = full_clin, pval = TRUE, conf.int = TRUE,
                     palette = c("#7da6c6", "#D97373"), title = "Overall Survival by Subtype")
  pdf("Fig_5G_KM_Survival.pdf", width = 6, height = 6, onefile = FALSE); print(p_km); dev.off()
  
  # [Fig 5H] DEGs (Group 2 vs Group 1)
  design <- model.matrix(~ 0 + factor(full_clin$Subtype))
  colnames(design) <- c("Group1", "Group2")
  fit <- lmFit(expr_matrix[, full_clin$Row.names], design)
  contrast <- makeContrasts(Group2 - Group1, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  deg_res <- topTable(fit2, coef = 1, number = Inf)
  
  # Volcano Plot
  deg_res$Significance <- "NS"
  deg_res$Significance[deg_res$logFC > 1 & deg_res$adj.P.Val < 0.05] <- "Up in Group 2"
  deg_res$Significance[deg_res$logFC < -1 & deg_res$adj.P.Val < 0.05] <- "Down in Group 2"
  
  p_volcano <- ggplot(deg_res, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.8) + scale_color_manual(values = c("blue", "grey", "red")) +
    theme_bw() + labs(title = "DEGs: Group 2 vs Group 1")
  ggsave("Fig_5H_Volcano.pdf", p_volcano, width = 6, height = 5)
  
  # [Fig 5J] Clinical Heatmap (Age, Gender, Stage, Grade)
  clin_anno <- full_clin[, c("Subtype", "Age", "Gender", "Stage")]
  clin_anno <- clin_anno[order(clin_anno$Subtype), ]
  ha <- HeatmapAnnotation(df = clin_anno, col = list(
    Subtype = c("Group1" = "#7da6c6", "Group2" = "#D97373"),
    Stage = c("I" = "#FEE0D2", "II" = "#FC9272", "III" = "#DE2D26", "IV" = "#99000D")
  ))
  
  mat_plot <- t(scale(t(expr_matrix[CANDIDATE_59_GENES[1:10], rownames(clin_anno)]))) # Top 10 genes for visual
  pdf("Fig_5J_Clinical_Heatmap.pdf", width=8, height=6)
  print(Heatmap(mat_plot, top_annotation = ha, show_column_names = FALSE, cluster_columns = FALSE))
  dev.off()
}

# ------------------------------------------------------------------------------
# PART 4: Tumor Microenvironment & Immune Infiltration (Fig 5K - 5L)
# ------------------------------------------------------------------------------
run_tme_analysis <- function(cibersort_res, estimate_res, subtypes_df) {
  message(">>> [Part 4] TME Profiling (CIBERSORT & ESTIMATE)...")
  
  # Merge all TME info
  tme_df <- merge(subtypes_df, cibersort_res, by.x = "Sample", by.y = "row.names")
  tme_df <- merge(tme_df, estimate_res, by.x = "Sample", by.y = "row.names")
  
  # [Fig 5K] CIBERSORT Violin plots
  ciber_melt <- tme_df %>% pivot_longer(cols = starts_with("B_cell"):starts_with("Neutrophils"), names_to = "Immune_Cell", values_to = "Fraction")
  
  p_ciber <- ggplot(ciber_melt, aes(x = Immune_Cell, y = Fraction, fill = Subtype)) +
    geom_violin(position = position_dodge(0.8), scale = "width") +
    geom_boxplot(width=0.1, position = position_dodge(0.8), outlier.shape = NA) +
    scale_fill_manual(values = c("#7da6c6", "#D97373")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
    stat_compare_means(aes(group = Subtype), label = "p.signif") + labs(title="CIBERSORT Infiltration")
  ggsave("Fig_5K_CIBERSORT_Violin.pdf", p_ciber, width = 12, height = 6)
  
  # [Fig 5L] ESTIMATE Scores
  estim_melt <- tme_df %>% pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"), names_to = "ScoreType", values_to = "Score")
  
  p_estim <- ggplot(estim_melt, aes(x = ScoreType, y = Score, fill = Subtype)) +
    geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill="white") +
    scale_fill_manual(values = c("#7da6c6", "#D97373")) + theme_bw() +
    stat_compare_means(method = "wilcox.test") + labs(title = "TME ESTIMATE Scores")
  ggsave("Fig_5L_ESTIMATE_Scores.pdf", p_estim, width = 6, height = 5)
}

# ==============================================================================
# EXECUTION BLOCK
# ==============================================================================
# In a real run, replace these randomly generated placeholder matrices with your actual TCGA load:
# load("TCGA_Pancancer_Expression.Rdata") -> expr_mat
# load("TCGA_Clinical.Rdata") -> clin_mat
# load("Cibersort_Estimate.Rdata") -> cib_res, est_res

message("Master Pipeline Setup Complete. Ready to execute functions with TCGA Data.")
# msag_35_genes <- run_msag_identification(expr_mat, clin_mat)
# subtypes <- run_subtyping_and_pca(expr_mat, msag_35_genes)
# run_clinical_deg_evaluation(expr_mat, clin_mat, subtypes)
# run_tme_analysis(cib_res, est_res, subtypes)