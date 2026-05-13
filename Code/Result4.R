# ==============================================================================
# Title: Master Pipeline for Result 4 - Mono/Macro Trajectories & Regulons
# Journal: Gigascience 
# Description: This pipeline comprehensively performs:
#              Part 1: Monocle3 Trajectory Inference (Fig 4A-I)
#              Part 2: Senescence-Trajectory Genes Intersection & GO (Fig 4J-L)
#              Part 3: Co-expression Modules & Dynamics (Fig 4M-O)
#              Part 4: SCENIC Regulon Inference Networks (Fig 4P-R)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages & Global Setup
# ------------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(monocle3)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(ggvenn)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(networkD3)
  library(SCENIC)
  library(pheatmap)
})

WORK_DIR <- "~/JH/Pancancer"
setwd(WORK_DIR)
# Assume SAG_union2 is a character vector of the 769 senescence genes
# load("SAG_union2.Rdata")

# ------------------------------------------------------------------------------
# 1. Data Cleaning & Preparation
# ------------------------------------------------------------------------------
message(">>> Step 1: Data Preparation...")
load("Mono_Macro.Rdata") # Load your Seurat object

# Remove Unidentified and Proliferating (MKI67) cells per biological logic
Mono_Macro <- subset(Mono_Macro, cell_anno != "C9_Mac_Unident" & cell_anno != "C10_Mac_MKI67")

# Merge Mono_CD14 and Mono_CD16 to serve as a unified trajectory root
Mono_Macro$cell_anno_new <- Mono_Macro$cell_anno
Mono_Macro$cell_anno_new[Mono_Macro$cell_anno %in% c("C0_Mono_CD14", "C8_Mono_CD16")] <- "C0_Mono_CD14/C8_Mono_CD16"

# ------------------------------------------------------------------------------
# 2. Part 1: Monocle3 Trajectory Inference Module (Fig 4A - 4I)
# ------------------------------------------------------------------------------
message(">>> Step 2: Monocle3 Trajectory Inference...")

# Helper function to auto-detect the root node corresponding to the Monocyte start
get_earliest_principal_node <- function(cds, time_bin="C0_Mono_CD14/C8_Mono_CD16"){
  cell_ids <- which(colData(cds)[, "cell_anno_new"] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}

# Core function to process monocle3 for any tissue group
run_monocle_pipeline <- function(seurat_obj, group_name) {
  message(sprintf("Processing Monocle3 for: %s", group_name))
  out_dir <- paste0("Result4_Output/", group_name)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Subset and build CDS
  sub_obj <- subset(seurat_obj, group == group_name)
  expr_matrix <- GetAssayData(sub_obj, assay = 'RNA', slot = 'counts')
  cell_meta <- sub_obj@meta.data
  gene_meta <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))
  
  cds <- new_cell_data_set(expr_matrix, cell_metadata = cell_meta, gene_metadata = gene_meta)
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  
  # Order cells dynamically setting the Monocyte cluster as root
  root_node <- get_earliest_principal_node(cds)
  cds <- order_cells(cds, root_pr_nodes = root_node)
  
  # [Fig 4A/C/E] Trajectory UMAP Plot
  p_traj <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) +
    theme_bw(base_size = 14) + labs(title = paste0(group_name, " - Pseudotime Trajectory"))
  ggsave(file.path(out_dir, paste0("Fig4_Trajectory_", group_name, ".pdf")), p_traj, width = 6, height = 5)
  
  # Extract pseudotime data
  sub_obj$Pseudotime <- pseudotime(cds)
  meta_df <- sub_obj@meta.data %>% filter(is.finite(Pseudotime))
  
  # [Fig 4B/D/F] Subpopulation Pseudotime Ranking
  rank_df <- meta_df %>% group_by(cell_anno) %>% summarise(mean_pt = mean(Pseudotime)) %>% arrange(mean_pt)
  meta_df$cell_anno <- factor(meta_df$cell_anno, levels = rank_df$cell_anno)
  p_rank <- ggplot(meta_df, aes(x = cell_anno, y = Pseudotime, fill = cell_anno)) +
    geom_boxplot(outlier.shape = NA) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "none") +
    labs(title = paste0(group_name, " - Mean Pseudotime Ranking"))
  ggsave(file.path(out_dir, paste0("Fig4_Pseudotime_Rank_", group_name, ".pdf")), p_rank, width = 7, height = 5)
  
  # [Fig 4G-I] Spearman correlation: Pseudotime vs SAG score
  p_cor <- ggplot(meta_df, aes(x = Pseudotime, y = SAG_union2)) +
    geom_point(alpha = 0.5, color = "#e68b81") +
    geom_smooth(method = "lm", color = "black") +
    stat_cor(method = "spearman", label.x.npc = "left") +
    theme_bw(base_size = 14) + labs(title = paste0(group_name, " - SAG vs Pseudotime"))
  ggsave(file.path(out_dir, paste0("Fig4_Correlation_", group_name, ".pdf")), p_cor, width = 5, height = 5)
  
  # Compute Trajectory-Associated Genes (graph_test)
  diff_test <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
  sig_traj_genes <- row.names(subset(diff_test, q_value < 0.05))
  
  return(list(cds = cds, seurat = sub_obj, traj_genes = sig_traj_genes))
}

# Run for all three tissues
res_tumor <- run_monocle_pipeline(Mono_Macro, "Tumor")
res_adj   <- run_monocle_pipeline(Mono_Macro, "Adjacent")
res_heal  <- run_monocle_pipeline(Mono_Macro, "Healthy")


# ------------------------------------------------------------------------------
# 3. Part 2: Intersection & Functional Enrichment (Fig 4J - 4L)
# ------------------------------------------------------------------------------
message(">>> Step 3: Intersection and Enrichment Analysis...")

# [Fig 4J] Intersect with SAG_union2 for each tissue
driver_tumor <- intersect(res_tumor$traj_genes, SAG_union2)
driver_adj   <- intersect(res_adj$traj_genes, SAG_union2)
driver_heal  <- intersect(res_heal$traj_genes, SAG_union2)

# [Fig 4K] Find core senescence key genes (Shared across contexts)
core_key_genes <- intersect(driver_heal, union(driver_tumor, driver_adj)) # 49 genes

pdf("Result4_Output/Fig4K_Venn_CoreGenes.pdf", width=6, height=5)
ggvenn(list(Healthy = driver_heal, Tumor_Adj = union(driver_tumor, driver_adj)), fill_color = c("#7da6c6", "#e68b81"))
dev.off()

# [Fig 4L] Functional enrichment of 49 Core Genes
entrez_ids <- mget(core_key_genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
enrich_res <- enrichGO(gene = na.omit(as.character(entrez_ids)), OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)

pdf("Result4_Output/Fig4L_Enrichment_CoreGenes.pdf", width=8, height=6)
print(barplot(enrich_res, showCategory = 10, title = "Core Senescence Drivers (GO BP)"))
dev.off()


# ------------------------------------------------------------------------------
# 4. Part 3: Co-expression Modules & Dynamics (Fig 4M - 4O)
# ------------------------------------------------------------------------------
message(">>> Step 4: Gene Co-expression Modules...")

# Example using Tumor tissue for module analysis (Matches text logic)
cds_tumor <- res_tumor$cds
gene_module_df <- find_gene_modules(cds_tumor[res_tumor$traj_genes,], resolution = 1e-2)

# [Fig 4M] Module Heatmap across Subpopulations
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_tumor)), cell_group=colData(cds_tumor)$cell_anno)
agg_mat <- aggregate_gene_expression(cds_tumor, gene_module_df, cell_group_df)
pdf("Result4_Output/Fig4M_Module_Heatmap_Tumor.pdf", width = 8, height = 7)
pheatmap(agg_mat, scale = "row", cluster_cols = TRUE, main = "Module Expression in Tumor")
dev.off()

# [Fig 4N] Sankey Diagram (Distribution of key genes across modules)
key_module_df <- gene_module_df %>% filter(id %in% core_key_genes)
nodes <- data.frame(name = unique(c(as.character(key_module_df$module), key_module_df$id)))
key_module_df$source <- match(key_module_df$module, nodes$name) - 1 
key_module_df$target <- match(key_module_df$id, nodes$name) - 1
key_module_df$value <- 1

sankey_plot <- sankeyNetwork(Links = key_module_df, Nodes = nodes, Source = "source",
                             Target = "target", Value = "value", NodeID = "name", fontSize = 12, nodeWidth = 30)
saveNetwork(sankey_plot, "Result4_Output/Fig4N_Sankey_Modules.html")


# ------------------------------------------------------------------------------
# 5. Part 4: Transcriptional Regulators Inference via SCENIC (Fig 4P - 4R)
# ------------------------------------------------------------------------------
message(">>> Step 5: SCENIC Regulon Inference...")

# FUNCTION: Run SCENIC safely in ISOLATED directories to prevent 'int/' overwriting
run_isolated_scenic <- function(seurat_obj, group_name) {
  scenic_dir <- paste0("Result4_Output/SCENIC_", group_name)
  if(!dir.exists(scenic_dir)) dir.create(scenic_dir, recursive = TRUE)
  setwd(scenic_dir) # CRITICAL: Change working directory
  
  exprMat <- as.matrix(GetAssayData(seurat_obj, assay="RNA", slot="counts"))
  cellInfo <- seurat_obj@meta.data
  
  # Initialize settings (Ensure your rcisTarget databases are in a reachable path)
  scenicOptions <- initializeScenic(org="hgnc", dbDir="~/JH/cisTarget_databases", nCores=10)
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  
  # Build Network & Score
  genesKept <- geneFiltering(exprMat, scenicOptions)
  exprMat_filtered <- exprMat[genesKept, ]
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- log2(exprMat_filtered + 1)
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  exprMat_log <- log2(exprMat + 1)
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
  scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
  
  setwd(WORK_DIR) # Restore global working directory after completion
  message(sprintf("SCENIC for %s completed successfully.", group_name))
}

# NOTE: SCENIC is extremely computationally intensive.
# It is recommended to run the following loop in a background HPC job.
# for (grp in c("Tumor", "Adjacent", "Healthy")) {
#   sub_obj <- subset(Mono_Macro, group == grp)
#   run_isolated_scenic(sub_obj, grp)
# }

message("=== Master Pipeline Complete! All trajectory, intersections, and modules are generated. ===")