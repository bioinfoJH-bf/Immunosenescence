# ==============================================================================
# Machine Learning Pipeline for Prognostic Signature Development (Mime1)
# ==============================================================================

# --- 1. Install Core Dependencies (Bioconductor) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva', 
               'ComplexHeatmap', 'survival', 'survminer', 'timeROC')

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE)
  }
}

# --- 2. Install ML Frameworks and Specific Algorithm Packages ---
# Note: randomForestSRC version 3.3.1 is recommended for var.select compatibility
if (!requireNamespace("CoxBoost", quietly = TRUE)) devtools::install_github("binderh/CoxBoost")
if (!requireNamespace("fastAdaboost", quietly = TRUE)) devtools::install_github("souravc83/fastAdaboost")
if (!requireNamespace("Mime1", quietly = TRUE)) devtools::install_github("l-magnificence/Mime")

# Load libraries
library(Mime1)
library(randomForestSRC)
library(survival)
library(survminer)
library(timeROC)
library(gbm)
library(dplyr)

# --- 3. Data Preprocessing ---
# 1. Import expression matrix (Assumes 'pan_mRNA' object is already in environment)
# 2. Import clinical info and match samples
survival <- read.table('~/JH/TCGA GDC/Pancancer数据/Survival_SupplementalTable_S1_20171025_xena_sp', 
                       sep='\t', header = TRUE)

# Standardize sample names (replace hyphens with dots to match R's data frame column naming)
survival$sample <- gsub("-", ".", survival$sample)
survival <- survival[survival$sample %in% rownames(pan_mRNA), ]

# Align and merge expression data with clinical data
pan_mRNA <- pan_mRNA[match(survival$sample, rownames(pan_mRNA)), ]
sur_info <- survival[, c('sample', 'OS.time', 'OS')]
expr_clinical <- merge(sur_info, pan_mRNA, by.x = "sample", by.y = "row.names")
colnames(expr_clinical)[1] <- 'ID'

# Data Cleaning: Remove missing values and samples with non-positive survival time
expr_clinical <- expr_clinical[complete.cases(expr_clinical) & expr_clinical$OS.time > 0, ]

# Split into Training and Validation sets (70:30 ratio)
set.seed(5201314)
train_index <- sample(1:nrow(expr_clinical), size = 0.7 * nrow(expr_clinical))
list_train_vali_Data <- list(
  Dataset1 = expr_clinical[train_index, ],
  Dataset2 = expr_clinical[-train_index, ]
)

# --- 4. Model Training Preparation ---
# Load candidate gene list (e.g., results from Univariate Cox regression)
rt <- read.table("uniCox.txt", header=TRUE, sep="\t", check.names=FALSE, row.names=1)
genelist <- rownames(rt)

# Fix compatibility issues for randomForestSRC
var.select <- function(object, ...) { UseMethod("var.select") }
var.select.rfsrc <- randomForestSRC:::var.select.rfsrc

# Run Mime1 Auto-ML Pipeline (including up to 117 combinations)
res <- Mime1::ML.Dev.Prog.Sig(
  train_data = list_train_vali_Data$Dataset1,
  list_train_vali_Data = list_train_vali_Data,
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = 0.05,
  candidate_genes = genelist,
  mode = 'all', 
  nodesize = 5,
  seed = 123 
)

# Save the resulting model object
save(res, file='res.Rdata')

# --- 5. Visualization & Evaluation ---

# 1. Comprehensive C-index Visualization
pdf('ML_Cindex_all.pdf', width = 6, height = 12)
cindex_dis_all(res, 
               validate_set = names(list_train_vali_Data)[-1], 
               order = names(list_train_vali_Data), 
               width = 0.35)
dev.off()

# 2. Calculate AUC for all models (Time-consuming operation)
all.auc.1y <- cal_AUC_ml_res(res, train_data = list_train_vali_Data$Dataset1, 
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 1)
all.auc.3y <- cal_AUC_ml_res(res, train_data = list_train_vali_Data$Dataset1, 
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 3)
all.auc.5y <- cal_AUC_ml_res(res, train_data = list_train_vali_Data$Dataset1, 
                             inputmatrix.list = list_train_vali_Data, mode = 'all', AUC_time = 5)

# 3. Survival Curve Comparison for a specific model (e.g., GBM)
survplot_list <- lapply(names(list_train_vali_Data), function(ds) {
  rs_sur(res, model_name = "GBM", dataset = ds, median.line = "hv", cutoff = 0.5, conf.int = TRUE)
})
aplot::plot_list(gglist = survplot_list, ncol = 2)

# --- 6. External Validation Function ---
validate_external_model <- function(ext_data, model_res, combination_name, cohort_label) {
  
  # 1. Extract sub-model and relevant genes
  model_obj <- model_res$ml.res[[combination_name]]
  
  # Extract genes based on different algorithmic prefixes
  if(grepl("StepCox", combination_name)) {
    model_genes <- attr(model_res$ml.res[['StepCox[forward]']]$formula, "term.labels")
  } else if(grepl("RSF", combination_name)) {
    model_genes <- model_res$ml.res[['RSF']]$forest$xvar.names
  } else {
    model_genes <- model_obj$fit$var.names
  }
  
  # 2. Align genes (Fill missing genes with 0)
  missing_genes <- setdiff(model_genes, colnames(ext_data))
  for (g in missing_genes) ext_data[[g]] <- 0
  
  # 3. Predict Risk Score
  test_df <- ext_data[, c("OS.time", "OS", model_genes)]
  test_df$risk <- predict(model_obj$fit, newdata = test_df, n.trees = model_obj$best)
  
  # 4. Calculate C-index
  cindex <- survival::survConcordance(Surv(OS.time, OS) ~ risk, data = test_df)$concordance
  message(paste(cohort_label, "-", combination_name, "C-index:", round(cindex, 4)))
  
  # 5. Kaplan-Meier Curve Generation
  test_df$group <- factor(ifelse(test_df$risk > median(test_df$risk), "High", "Low"), 
                          levels = c("Low", "High"))
  fit <- survfit(Surv(OS.time, OS) ~ group, data = test_df)
  
  p_km <- ggsurvplot(fit, data = test_df, pval = TRUE, conf.int = TRUE, 
                     risk.table = TRUE, title = paste(cohort_label, "Validation"),
                     palette = c("#868686", "#B24745"), xlab = "Days")
  
  # 6. Time-Dependent ROC Calculation
  roc_res <- timeROC(T = test_df$OS.time, delta = test_df$OS, marker = test_df$risk,
                     cause = 1, times = c(365, 1095, 1825), weighting = "marginal")
  
  return(list(plot_km = p_km, roc = roc_res, cindex = cindex))
}

# --- Examples: Running External Validation ---

# Example 1: GSE74777
# gse74_res <- validate_external_model(GSE74777, res, "StepCox[forward] + GBM", "GSE74777")
# print(gse74_res$plot_km)

# Example 2: GSE72094
# GSE72094_clean <- read.csv("GSE72094_exp+survival.csv") %>% 
#   rename(ID = 1, OS = 2, OS.time = 3) # Adjust column indices based on your actual file
# gse72_res <- validate_external_model(GSE72094_clean, res, "StepCox[forward] + GBM", "GSE72094")