# ==============================================================================
# Title: Master Pipeline for Result 2 - Pan-Cancer TME Senescence Interactions
# Journal: Gigascience
# Description: This unified script executes the three core modules of Result 2:
#              Part 1: Malignant Epithelial Cell Identification via InferCNV
#              Part 2: Senescence-Associated Cell Communication via CellChat
#              Part 3: Ligand-Receptor Expression Profiling in Myeloid Subsets
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Global Setup & Library Loading
# ------------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(infercnv)
  library(CellChat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggvenn)
  library(patchwork)
  library(RColorBrewer)
})

# Define global variables
GENE_ORDER_FILE <- "~/JH/inferCNV/hg38_gencode_v27.txt"
SAG_UNION2 <- c("APP", "MIF", "EGF", "GRN", "PLAU", "AGT", "NGF", "CD74", "CD44") # Replace with complete 769 SAG vector
CUSTOM_COLORS <- c('Macro'='#7da6c6', 'Mono'='#e68b81', 'DC'='#B7B2D0', 'Mast'='#EAAA60')

# ------------------------------------------------------------------------------
# PART 1: Malignant Epithelial Cell Identification (InferCNV)
# ------------------------------------------------------------------------------
run_cnv_clustering <- function(sc, cancer_type, out_dir) {
  message(sprintf("\n[Part 1] Starting CNV Analysis for %s...", cancer_type))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Subset to relevant lineages
  Idents(sc) <- "cell_type"
  sc_sub <- subset(sc, idents = c("Epithelium", "Lymphocyte", "Myeloid"))
  
  # Prepare Annotation
  gene_order <- read.table(GENE_ORDER_FILE, sep = '\t', row.names = 1)
  gene_file <- gene_order[intersect(rownames(sc_sub), rownames(gene_order)), ] %>% na.omit()
  
  # Run InferCNV
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = sc_sub@assays$RNA@counts,
    annotations_file  = as.matrix(sc_sub@active.ident),
    gene_order_file   = gene_file, delim = '\t',
    ref_group_names   = c("Lymphocyte", "Myeloid")
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj, cutoff = 0.1, out_dir = out_dir, cluster_by_groups = TRUE,
    denoise = TRUE, num_threads = 10, HMM = TRUE, analysis_mode = "subclusters",
    leiden_resolution = 0.0001, output_format = "pdf"
  )
  
  # K-Means Clustering for Epithelial Malignancy (Epithelium 1-7)
  expr <- readRDS(file.path(out_dir, "run.final.infercnv_obj"))@expr.data
  test_loc <- readRDS(file.path(out_dir, "run.final.infercnv_obj"))@observation_grouped_cell_indices$Epithelium
  expr_obs <- expr[, colnames(expr)[test_loc], drop = FALSE]
  
  cnv_score_df <- data.frame(CB = colnames(expr_obs), CNV_score = colMeans((expr_obs - 1)^2))
  
  set.seed(20210418)
  kmeans_res <- kmeans(t(expr_obs), centers = 7)
  cnv_score_df$kmeans_cluster <- as.factor(kmeans_res$cluster)
  
  # Order by mean CNV score
  cluster_means <- cnv_score_df %>% group_by(kmeans_cluster) %>% summarise(mean = mean(CNV_score)) %>% arrange(mean)
  mapping <- setNames(paste0("Epithelium", 1:7), cluster_means$kmeans_cluster)
  cnv_score_df$Epithelial_Subgroup <- factor(mapping[as.character(cnv_score_df$kmeans_cluster)], levels = paste0("Epithelium", 1:7))
  
  # Save Violin Plot (Fig 2B)
  p <- ggplot(cnv_score_df, aes(x = Epithelial_Subgroup, y = CNV_score, fill = Epithelial_Subgroup)) +
    geom_violin(color = NA, alpha = 0.8) + scale_fill_brewer(palette = "Set3") + theme_bw(base_size=14) +
    stat_compare_means(method = 'anova', label.y.npc = 0.95) + theme(legend.position = "none") +
    labs(title = paste0(cancer_type, " - Malignant Subgroups"), x = NULL, y = "CNV Score")
  ggsave(file.path(out_dir, paste0(cancer_type, "_Fig2B_CNV_Subgroups.pdf")), plot = p, width = 7, height = 5)
  
  return(cnv_score_df)
}

# ------------------------------------------------------------------------------
# PART 2: Intercellular Communication Analysis (CellChat)
# ------------------------------------------------------------------------------
run_cellchat_analysis <- function(sc, cnv_meta, cancer_type, out_dir) {
  message(sprintf("\n[Part 2] Starting CellChat Analysis for %s...", cancer_type))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  setwd(out_dir)
  
  # Integrate CNV Subgroups into Seurat Meta
  rownames(cnv_meta) <- cnv_meta$CB
  sc <- AddMetaData(sc, metadata = cnv_meta[, "Epithelial_Subgroup", drop = FALSE], col.name = "Epithelium_Malignancy")
  sc$CellChat_Labels <- ifelse(is.na(sc$Epithelium_Malignancy), as.character(sc$cell_type), as.character(sc$Epithelium_Malignancy))
  Idents(sc) <- "CellChat_Labels"
  
  # Build CellChat Object
  cellchat <- createCellChat(object = sc, group.by = "CellChat_Labels")
  cellchat@DB <- CellChatDB.human 
  cellchat <- subsetData(cellchat) %>% identifyOverExpressedGenes() %>% identifyOverExpressedInteractions()
  cellchat <- computeCommunProb(cellchat, type = "triMean") %>% filterCommunication(min.cells = 10) %>% computeCommunProbPathway() %>% aggregateNet()
  
  # Find Senescence-Associated Pathways (APP, MIF etc.)
  active_pathways <- cellchat@netP$pathways
  senescence_pathways <- intersect(active_pathways, SAG_UNION2)
  
  # Plot Venn Diagram (Fig 2E)
  pdf(paste0(cancer_type, "_Fig2E_Venn_SAG.pdf"), width = 8, height = 5)
  print(ggvenn(list("CellChat Pathways" = active_pathways, "Senescence Genes" = SAG_UNION2), fill_color = c("#5E81AC", "#EBCB8B")))
  dev.off()
  
  # Plot Target Sender/Receiver Scatter (Fig 2H)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  for (pathway in c("APP", "MIF")) {
    if (pathway %in% active_pathways) {
      p_path <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathway, title = paste('Role in', pathway))
      p_all <- netAnalysis_signalingRole_scatter(cellchat, title = 'All pathways')
      ggsave(paste0(cancer_type, "_Fig2H_Scatter_", pathway, ".pdf"), plot = p_path + p_all, width = 10, height = 5)
    }
  }
  return(cellchat)
}

# ------------------------------------------------------------------------------
# PART 3: Gene Expression Visualization (Receptor Distributions)
# ------------------------------------------------------------------------------
plot_receptor_expression <- function(sc_myeloid, cancer_type, out_dir, target_genes = c("CD74", "CD44")) {
  message(sprintf("\n[Part 3] Plotting Receptor Expression for %s...", cancer_type))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  for (gene in target_genes) {
    if (gene %in% rownames(sc_myeloid)) {
      # Extract Gene Expression
      sc_myeloid@meta.data[[gene]] <- GetAssayData(sc_myeloid, assay = "RNA", slot = "data")[gene, ]
      
      # Sort subsets by mean expression of the specific gene
      ordered_subtypes <- sc_myeloid@meta.data %>% 
        group_by(cell_subtype) %>% 
        summarise(mean_exp = mean(!!sym(gene), na.rm = TRUE)) %>% 
        arrange(mean_exp) %>% pull(cell_subtype)
      
      sc_myeloid@meta.data$cell_subtype <- factor(sc_myeloid@meta.data$cell_subtype, levels = ordered_subtypes)
      
      # Generate Boxplot with Error Bars (Fig 2K/2L Style)
      p_box <- ggplot(sc_myeloid@meta.data, aes(x = cell_subtype, y = !!sym(gene))) +  
        stat_boxplot(aes(color = cell_subtype), geom = "errorbar", width = 0.3, size = 0.6) + 
        geom_boxplot(aes(fill = cell_subtype, color = cell_subtype), outlier.shape = 18, size = 0.6, alpha=0.8) +  
        scale_fill_manual(values = CUSTOM_COLORS) +  
        scale_color_manual(values = CUSTOM_COLORS) +
        theme_bw() +  
        theme(
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1),
          axis.title = element_text(face = "bold", size = 12),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none"
        ) +
        labs(title = paste0(cancer_type, " - ", gene, " Expression in Myeloid Cells"), x = "Myeloid Subtypes", y = paste(gene, "Expression Level"))
      
      ggsave(file.path(out_dir, paste0(cancer_type, "_Fig2L_", gene, "_Expression.pdf")), plot = p_box, width = 6, height = 5)
    }
  }
}

# ------------------------------------------------------------------------------
# 4. Master Execution Loop (Run all cancers)
# ------------------------------------------------------------------------------
# Map your cancer types to their processed Seurat object paths
cancer_datasets <- list(
  "ccRCC"  = "~/JH/GSE207493_RAW/sc_processed.Rdata",
  "Breast" = "~/JH/GSE176078_RAW/(Breast)sc(质控后).Rdata"
  # Add CRC, Gastric, NSCLC, Ova, Prostate, RCC here...
)

# Run Pan-Cancer Master Pipeline
for (cancer in names(cancer_datasets)) {
  dataset_path <- cancer_datasets[[cancer]]
  if (!file.exists(dataset_path)) next
  
  # Load Data
  load_env <- new.env()
  load(dataset_path, envir = load_env)
  sc_obj <- load_env$sc 
  
  # Define Output Directories
  out_cnv      <- paste0("~/JH/Result2_Output/", cancer, "/1_CNV")
  out_cellchat <- paste0("~/JH/Result2_Output/", cancer, "/2_CellChat")
  out_expr     <- paste0("~/JH/Result2_Output/", cancer, "/3_Expression")
  
  # Step 1: InferCNV Malignancy Clustering
  cnv_metadata <- run_cnv_clustering(sc = sc_obj, cancer_type = cancer, out_dir = out_cnv)
  
  # Step 2: CellChat Communication Analysis
  cellchat_res <- run_cellchat_analysis(sc = sc_obj, cnv_meta = cnv_metadata, cancer_type = cancer, out_dir = out_cellchat)
  
  # Step 3: Expression Visualization (Filter Myeloid Cells first)
  if ("cell_type" %in% colnames(sc_obj@meta.data)) {
    sc_myeloid <- subset(sc_obj, subset = cell_type %in% c("Myeloid", "Mono", "Macro", "DC"))
    plot_receptor_expression(sc_myeloid = sc_myeloid, cancer_type = cancer, out_dir = out_expr, target_genes = c("CD74", "CD44", "APP", "MIF"))
  }
}
message("=== Master Pipeline Complete! All outputs have been successfully generated. ===")