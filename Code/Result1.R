# ==============================================================================
# Title: Pan-cancer Immunosenescence Analysis - GTEx Physiological Aging
# Journal: Gigascience
# Description: This script processes bulk RNA-seq data from multiple normal 
#              tissues (GTEx), calculates ssGSEA scores for immune cell subsets, 
#              and generates radar plots and stacked bar plots to illustrate 
#              myeloid-biased immune aging trajectories.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(GSVA)
library(GSEABase)
library(reshape2)
library(fmsb)
library(ggplot2)

# Set global working directories
WORK_DIR <- "~/JH/Pancancer"
GTEX_DIR <- "~/JH/GTEx"

# ------------------------------------------------------------------------------
# 2. Gene Set Preparation (Generate GMT file)
# ------------------------------------------------------------------------------
setwd(WORK_DIR)

# Define marker genes based on Breast 078 dataset and literature
Immune <- c("CCL4","CXCR4","CCL5","CD69","SRGN","IL7R","RGS1","TYROBP","CD52","DUSP2")
Endothelium <- c("ACKR1","FABP4","PLVAP","RAMP2","VWF","AQP1","CLDN5","SPARCL1","GNG11","PECAM1")
Epithelium <- c("SCGB2A2","CD24","MUCL1","KRT19","PIP","AZGP1","SCGB1D2","KRT18","CLDN4","AGR2")
Fibroblast <- c("COL1A1","COL1A2","DCN","TAGLN","COL3A1","LUM","ACTA2","SFRP2","MYL9","APOD")

marker_list <- list(
  Immune = Immune,
  Endothelium = Endothelium,
  Epithelium = Epithelium,
  Fibroblast = Fibroblast
)

marker_df <- do.call(rbind, lapply(names(marker_list), function(celltype) {
  data.frame(cluster = celltype, gene = marker_list[[celltype]])
})) %>% as_tibble()

# Create GMT format lines
cell_types <- unique(marker_df$cluster)
gmt_lines <- sapply(cell_types, function(ct) {
  genes <- marker_df$gene[marker_df$cluster == ct]
  paste(ct, "NA", paste(genes, collapse = "\t"), sep = "\t")
})

# Write GMT files
writeLines(gmt_lines, "celltype_DEG（Breast筛选Immune和大类）.gmt")

# Load the gene set for downstream GSVA
geneSet <- getGmt("celltype_DEG（Breast筛选Immune和大类）.gmt", geneIdType=SymbolIdentifier())


# ------------------------------------------------------------------------------
# 3. Define Pipeline for GTEx Data Processing & Radar Plot Generation
# ------------------------------------------------------------------------------
# Normalize function for ssGSEA scores
normalize_score <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Core function to process each tissue and draw radar plot
process_and_plot_tissue <- function(tissue, gtex_dir, gene_set) {
  message(paste0("Processing tissue: ", tissue, "..."))
  
  # 3.1 Read data
  file_path <- file.path(gtex_dir, paste0(tissue, "_reads.txt"))
  if(!file.exists(file_path)) {
    warning(paste0("File not found: ", file_path))
    return(NULL)
  }
  
  raw_data <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  # Extract expression matrix (remove first 3 rows which contain metadata)
  expr_data <- raw_data[-c(1:3), ]
  mat <- matrix(as.numeric(as.matrix(expr_data)), nrow = nrow(expr_data), 
                dimnames = list(rownames(expr_data), colnames(expr_data)))
  
  # 3.2 GSVA Analysis (Updated for GSVA >= 1.52)
  gsvapar <- gsvaParam(mat, gene_set, kcdf = 'Gaussian', absRanking = TRUE)
  ssgseaScore <- gsva(gsvapar)
  
  # Normalize scores
  rt_norm <- normalize_score(ssgseaScore)
  ssgseaOut <- t(rt_norm) %>% as.data.frame()
  ssgseaOut$ID <- rownames(ssgseaOut)
  
  # 3.3 Extract metadata (Age is assumed to be in the 2nd row based on original code)
  meta_rt <- raw_data[2, , drop = FALSE] 
  meta_rt <- t(meta_rt) %>% as.data.frame()
  colnames(meta_rt)[1] <- "Sample_age"
  meta_rt$ID <- rownames(meta_rt)
  
  # Merge data
  merged_df <- merge(meta_rt, ssgseaOut, by = "ID")
  rownames(merged_df) <- merged_df$ID
  merged_df$Sample <- merged_df$ID
  merged_df$ID <- NULL
  
  # 3.4 Reshape data
  df_melt <- melt(merged_df, id.vars = c("Sample", "Sample_age"),
                  variable.name = "CellType", value.name = "Score")
  
  # Define Age Groups (Young: 20-39, Old: 40-79+)
  df_melt$Group <- ifelse(df_melt$Sample_age %in% c("20-29", "30-39"), "Young", "Old")
  df_melt$Group <- factor(df_melt$Group, levels = c("Young", "Old"))
  df_melt$Score <- as.numeric(df_melt$Score)
  
  # Calculate mean scores per group
  df_radar <- df_melt %>%
    dplyr::group_by(Group, CellType) %>%
    dplyr::summarise(MeanScore = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = CellType, values_from = MeanScore)
  
  # Format for fmsb::radarchart (Max and Min rows required)
  if(nrow(df_radar) < 2) return(NULL) # Skip if only one age group exists
  
  max_row <- apply(df_radar[, -1], 2, max)
  min_row <- apply(df_radar[, -1], 2, min)
  radar_data <- rbind(max_row, min_row, df_radar[, -1])
  rownames(radar_data) <- c("Max", "Min", as.character(df_radar$Group))
  
  # 3.5 Generate Radar Plot
  colors <- c("#bdd7e7", "#F7BDA0")
  pty_vec <- rep(16, nrow(df_radar))
  pcol_vec <- colors[1:nrow(df_radar)]
  pfcol_vec <- alpha(colors[1:nrow(df_radar)], 0.5)
  plwd_vec <- rep(3, nrow(df_radar))
  
  pdf_file <- file.path(gtex_dir, paste0(tissue, "_Radar.pdf"))
  pdf(pdf_file, width = 7, height = 7)
  
  radarchart(
    radar_data,
    pty = pty_vec, axistype = 1, pcol = pcol_vec, pfcol = pfcol_vec, 
    plwd = plwd_vec, plty = 1, cglcol = "grey60", cglty = 1, cglwd = 1, 
    axislabcol = "grey60", vlcex = 0.8, vlabels = colnames(radar_data), 
    caxislabels = seq(0, 1, 0.2), calcex = 0.8
  )
  
  legend("topright", legend = rownames(radar_data)[3:nrow(radar_data)],
         col = pcol_vec, lty = 1, lwd = plwd_vec, pch = pty_vec, 
         bty = "n", cex = 0.9)
  
  dev.off()
}

# ------------------------------------------------------------------------------
# 4. Execute Pipeline Across All Tissues (Fig. 1E Radar Plots)
# ------------------------------------------------------------------------------
# List of all GTEx tissues based on your script
tissues_to_process <- c("Blood-Vessel", "Brain", "Colon", "Cervix-Uteri", 
                        "Breast", "Blood", "Adipose-Tissue", "Bladder", 
                        "Adrenal-Gland", "Fallopian-Tube", "Esophagus", 
                        "Kidney", "Heart", "Liver", "Lung", "Muscle", 
                        "Ovary", "Nerve", "Pancreas", "Pituitary", 
                        "Salivary-Gland", "Prostate", "Small-Intestine", 
                        "Spleen", "Stomach", "Skin", "Testis", "Uterus", "Vagina", "Thyroid")

# Run the loop
lapply(tissues_to_process, process_and_plot_tissue, gtex_dir = GTEX_DIR, gene_set = geneSet)


# ------------------------------------------------------------------------------
# 5. Cross-Tissue Frequency Analysis (Fig. 1B & 1C Stacked Bar Plots)
# ------------------------------------------------------------------------------
# Data matrix for Lymphocyte vs Myeloid enrichment frequency (Fig 1C)
data_lineage <- data.frame(
  row.names = c("Breast", "Prostate", "Lung", "Liver", "Colon", "Stomach", 
                "Esophagus", "Pancreas", "Cervix-Uteri", "Uterus", "Ovary", 
                "Skin", "Thyroid", "Bladder", "Kidney", "Brain", "Small-Intestine"),
  Myeloid = c(0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1),
  Lymphocyte = c(0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0)
)

# Convert to long format and filter enriched events
df_long_lineage <- melt(as.matrix(data_lineage))
colnames(df_long_lineage) <- c("Tissue", "CellType", "Value")
df_long_lineage <- subset(df_long_lineage, Value == 1)

# Color palette definition for 17 tissues
tissue_colors <- c(
  "Breast" = "#E9998E", "Prostate" = "#ECA79B", "Lung" = "#F0B4A8", 
  "Liver" = "#F3C3B5", "Colon" = "#F7D1C3", "Stomach" = "#EAAA60",
  "Esophagus" = "#F3C284", "Pancreas" = "#FDDAAA", "Cervix-Uteri" = "#A6D96A", 
  "Uterus" = "#C1E495", "Ovary" = "#DCF1C1", "Skin" = "#7da6c6", 
  "Thyroid" = "#A0C3DA", "Bladder" = "#C3E1EE", "Kidney" = "#B7B2D0", 
  "Brain" = "#CDC9DF", "Small-Intestine" = "#E3E0EE"
)

# Generate Stacked Bar Plot
fig1c <- ggplot(df_long_lineage, aes(x = CellType, fill = Tissue)) +
  geom_bar(position = "stack", color = "white", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = tissue_colors, guide = guide_legend(nrow = 3)) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 10),
    labels = function(x) ifelse(x %% 1 == 0, x, "") 
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Tissue Composition of Each Immune Lineage",
    x = "Immune Lineage",
    y = "Number of Enriched Tissues",
    fill = "Tissue"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(fig1c)

# Save the plot
ggsave("Fig_1C_Lineage_Frequency.pdf", plot = fig1c, width = 8, height = 6)
# ==============================================================================
# Title: Pan-cancer Immunosenescence Single-Cell Analysis Pipeline
# Journal: Gigascience
# Description: A unified pipeline for quality control, Harmony batch correction, 
#              and SAG (Senescence-Associated Genes) scoring across multiple 
#              tumor scRNA-seq datasets to demonstrate myeloid senescence bias.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggridges)
library(forcats)
library(modeest) # For density mode calculation

# ------------------------------------------------------------------------------
# 2. Define the Curated Senescence Gene Set (769 Genes)
# ------------------------------------------------------------------------------
# Replace this with your actual SAG_union2 character vector
# Example: SAG_union2 <- read.csv("SAG_769_genes.csv")$GeneSymbol
SAG_union2 <- c("Gene1", "Gene2", "Gene3") # Placeholder

# ------------------------------------------------------------------------------
# 3. Unified Single-Cell Processing Function
# ------------------------------------------------------------------------------
# Matches the text: "processed under a unified quality control and batch correction pipeline"
process_cancer_scRNA <- function(data_dir, cancer_name) {
  message(paste0("\n>>> Step 1: Processing Dataset - ", cancer_name, " ..."))
  
  # Read multi-sample 10X matrices
  sample_list <- basename(list.files(data_dir, recursive = FALSE))
  Object_list <- list()
  
  for (sample in sample_list) {
    filedir <- file.path(data_dir, sample)
    scrna_data <- Read10X(filedir)
    
    # Check if Read10X returns a list (multi-modal) and extract Gene Expression
    if (class(scrna_data) == "list") {
      scrna_data <- scrna_data$`Gene Expression`
    }
    
    obj <- CreateSeuratObject(counts = scrna_data, min.cells = 3, min.features = 200)
    obj$sample <- sample
    Object_list[[sample]] <- obj
  }
  
  # Merge all samples
  sc <- merge(Object_list[[1]], y = Object_list[-1])
  
  # Quality Control
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt < 20)
  
  # Normalization & Dimensionality Reduction
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc)
  sc <- ScaleData(sc)
  sc <- RunPCA(sc, verbose = FALSE)
  
  # Batch Correction using Harmony (Critical for pan-cancer integration)
  sc <- RunHarmony(sc, group.by.vars = "sample")
  
  # Clustering and UMAP using Harmony embeddings
  sc <- RunUMAP(sc, reduction = "harmony", dims = 1:30)
  sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30)
  sc <- FindClusters(sc, resolution = 0.5)
  
  return(sc)
}

# ------------------------------------------------------------------------------
# 4. SAG Scoring & Ridge Plot Function
# ------------------------------------------------------------------------------
# Matches the text: "Ridge plots show SAG scores for myeloid and lymphoid cells"
plot_senescence_ridge <- function(seurat_obj, sag_genes, cancer_name) {
  message(paste0(">>> Step 2: Scoring and Plotting - ", cancer_name, " ..."))
  
  # Add Module Score for Senescence Set
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- AddModuleScore(
    seurat_obj, 
    features = list(sag_genes), 
    ctrl = 100, 
    name = "SAG_score"
  )
  # Correcting column name injected by Seurat (appends '1')
  colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- 'SAG_union2' 
  
  # Extract specific metadata columns for plotting
  metadata_df <- as.data.frame(seurat_obj@meta.data[, c("cell_type", "SAG_union2")])
  metadata_df$cell_type <- as.factor(metadata_df$cell_type)
  
  # Function: Order cell types by the highest density peak (Mode) of the SAG score
  get_mode <- function(paramSampleVec) {
    mcmcDensity <- density(paramSampleVec, na.rm = TRUE)
    mo <- mcmcDensity$x[which.max(mcmcDensity$y)]
    return(mo)
  }
  
  metadata_df <- metadata_df %>%
    mutate(cell_type = fct_reorder(cell_type, SAG_union2, .fun = get_mode, .desc = FALSE))
  
  # Generate Ridge Plot
  p_ridge <- ggplot(metadata_df, aes(x = SAG_union2, y = cell_type, fill = cell_type, color = cell_type)) +
    geom_density_ridges(
      jittered_points = TRUE, scale = 2, rel_min_height = 0.01,
      point_shape = "|", point_size = 3, size = 0.25,
      position = position_points_jitter(height = 0)
    ) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0, 0), name = 'SAG Score (Senescence)') +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    labs(title = paste0(cancer_name, " - Immune Senescence (Myeloid vs Lymphoid)"), y = NULL)
  
  print(p_ridge)
  
  # Save the plot
  ggsave(filename = paste0("Fig1F_", cancer_name, "_SAG_RidgePlot.pdf"), 
         plot = p_ridge, width = 7, height = 5)
  
  return(seurat_obj)
}

# ------------------------------------------------------------------------------
# 5. Execute Pipeline Across All Cancers (Loop implementation)
# ------------------------------------------------------------------------------
# Dictionary containing all data directories
cancer_datasets <- list(
  "Breast"   = "~/JH/GSE176078_RAW",
  "ccRCC"    = "~/JH/GSE207493_RAW",
  "CRC"      = "~/JH/GSE166555_RAW",   # Please update paths
  "Gastric"  = "~/JH/GSE167297_RAW",
  "HNSCC"    = "~/JH/GSE139324_RAW",
  "NSCLC"    = "~/JH/GSE117570_RAW",
  "Ova"      = "~/JH/GSE184880_RAW",
  "PDAC"     = "~/JH/GSE205049_RAW",
  "Prostate" = "~/JH/GSE181294_RAW",
  "RCC"      = "~/JH/GSE159115_RAW"
)

# Run the master loop
for (cancer in names(cancer_datasets)) {
  
  data_path <- cancer_datasets[[cancer]]
  if (!dir.exists(data_path)) next # Skip if path is invalid
  
  # Step 1: Execute Unified Pipeline
  sc_obj <- process_cancer_scRNA(data_dir = data_path, cancer_name = cancer)
  
  # Step 2: Subset only Immune cells for Senescence analysis
  
  if ("cell_type" %in% colnames(sc_obj@meta.data)) {
    sc_immune <- subset(sc_obj, subset = cell_type %in% c("Myeloid", "Lymphocyte", "B cell", "T cell", "Macrophage", "Monocyte"))
    
    # Step 3: Add Module Score and generate Ridge Plot (Fig 1F)
    sc_immune <- plot_senescence_ridge(seurat_obj = sc_immune, sag_genes = SAG_union2, cancer_name = cancer)
    
    # Save the processed immune object
    saveRDS(sc_immune, file = paste0(cancer, "_immune_processed.rds"))
  } else {
    warning(paste0("'cell_type' missing for ", cancer, ". Please annotate first."))
  }
}