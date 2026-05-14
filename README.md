# Pan-cancer Immunosenescence Analysis

![Language](https://img.shields.io/badge/Language-R-blue.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

## Project Overview
This repository contains the core analytical code and pipelines for the study of **pan-cancer myeloid immunosenescence**. The study explores how the tumor microenvironment (TME) hijacks physiological myeloid senescence programs to drive immunosuppression and develops a robust prognostic stratification model (MSAG.SIG) using machine learning algorithms.

The analytical framework integrates bulk RNA-seq (GTEx/TCGA) and multi-cohort single-cell RNA sequencing (scRNA-seq) data to construct high-resolution myeloid atlases, cellular communication networks, and clinical prognostic signatures.

## Repository Structure
The analysis scripts are organized sequentially corresponding to the results sections of our manuscript. All core scripts are located in the `Code/` directory:

* `Code/`
  * [`Result1.R`](Code/Result1.R) - **Physiological Aging & Single-Cell Profiling:** Processes bulk RNA-seq from normal tissues (GTEx) to calculate ssGSEA scores illustrating myeloid-biased immune aging. Executes a unified pipeline for scRNA-seq quality control, Harmony batch correction, and Senescence-Associated Gene (SAG) scoring across multiple tumor datasets.
  * [`Result2.R`](Code/Result2.R) - **TME Senescence Interactions:** Identifies malignant epithelial cells via `inferCNV`, decodes senescence-associated intercellular communication using `CellChat`, and profiles ligand-receptor expression in distinct myeloid subsets.
  * [`Result3.R`](Code/Result3.R) - **High-Resolution Myeloid Atlas & In Silico Knockdown:** Constructs a pan-cancer Monocyte/Macrophage atlas (UMAP, GSVA functional scoring). Performs *in silico* virtual knockdown of specific targets (e.g., CD74) using `scTenifoldKnk` followed by downstream Gene Ontology (GO) enrichment.
  * [`Result4.R`](Code/Result4.R) - **Cellular Trajectories & Transcriptional Networks:** Infers pseudotime developmental trajectories using `Monocle3`. Computes gene co-expression modules and reconstructs single-cell gene regulatory networks (Regulons) via `SCENIC`.
  * [`Result5.R`](Code/Result5.R) - **MSAG Subtyping & Clinical Evaluation:** Translates scRNA-seq signatures to bulk TCGA cohorts. Conducts unsupervised consensus clustering based on the Myeloid Senescence-Associated Gene signature (MSAG), evaluating clinical outcomes (Kaplan-Meier survival) and immune infiltration landscapes (`CIBERSORT` & `ESTIMATE`).
  * [`Result6.R`](Code/Result6.R) - **Machine Learning Prognostic Modeling:** Develops and validates a clinical prognostic signature utilizing the `Mime1` auto-ML framework. Integrates up to 117 machine learning algorithm combinations (e.g., StepCox, RSF, GBM), calculating C-index and Time-ROC across independent validation cohorts.

## Prerequisites & Dependencies
The scripts are written entirely in **R**. The following major R and Bioconductor packages are required to execute the pipeline:

**Single-Cell & Spatial Analysis:**
* `Seurat` (v4)
* `harmony`
* `CellChat`
* `infercnv`
* `monocle3`
* `SCENIC`
* `scTenifoldKnk`

**Bulk RNA-seq & Clinical/Machine Learning:**
* `GSVA` / `GSEABase`
* `survival` / `survminer` / `timeROC`
* `ConsensusClusterPlus`
* `Mime1`
* `randomForestSRC`
* `CoxBoost` / `fastAdaboost` / `gbm`

**Visualization:**
* `ggplot2`, `ggpubr`, `ComplexHeatmap`, `fmsb`, `ggvenn`

## Usage
1. Clone this repository to your local machine:
   ```bash
   git clone [https://github.com/YourUsername/YourRepository.git](https://github.com/YourUsername/YourRepository.git)
