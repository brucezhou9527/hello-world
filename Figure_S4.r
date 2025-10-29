##############################################
# Figure S4 - Myeloid cell analysis
# Author: Bruce Zhou
# Description:
#   This script generates Supplementary Figure S4 for the manuscript,
#   including feature visualization of myeloid subsets and preparation
#   for NicheNet interaction analysis.
##############################################

# ======================
# 1. Environment setup
# ======================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GseaVis)
  library(GSEABase)
  library(enrichplot)
  library(GSVA)
  library(CellChat)
  library(nichenetr)
  library(tidyverse)
  library(pheatmap)
  library(scRNAtoolVis)
  library(ggsci)
})

# Set working directory and parallel settings
setwd("workpath")
plan("multicore", workers = 20)
options(future.globals.maxSize = 1000000 * 1024^2)
set.seed(1234)

cat("✅ Environment ready. Starting Figure S4 processing...\n")

# ======================
# 2. Data loading
# ======================
sce <- readRDS("myeloid.rds")
cat("📦 Data loaded:", ncol(sce), "cells ×", nrow(sce), "genes\n")

# ======================
# 3. Figure S4A: Feature plot
# ======================
# Define key marker genes for myeloid cell subsets
gene_feature <- c(
  "CD14", "FCGR3A", "FCGR3B", "CXCR2", "CSF3R",
  "CD1C", "FCER1A", "CLEC9A", "LAMP3",
  "THBS1", "C1QA", "APOE", "GPX3", "VCAN", "SPP1", "TREM2"
)

cat("🎯 Plotting", length(gene_feature), "marker genes...\n")

# Generate high-resolution feature panel
tiff("figure_S4A_featureplot.tiff",
     units = "in", width = 10, height = 10,
     res = 300, compression = "lzw")

FeatureCornerAxes(
  object = sce,
  reduction = "tsne",
  groupFacet = NULL,       # disable faceting
  relLength = 0.4,
  relDist = 0.2,
  pSize = 0.01,
  aspect.ratio = NULL,
  nLayout = 4,
  cornerTextSize = 2,
  cornerVariable = 1,
  minExp = 0, maxExp = 4,
  features = gene_feature
)

dev.off()
cat("✅ Figure S4A saved: figure_S4A_featureplot.tiff\n")

# ======================
# 4. Figure S4B: NicheNet setup
# ======================
# (Full analysis is performed in `figure_4.R`)

