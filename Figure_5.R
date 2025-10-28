##############################################################
# Title: EGR1+ B cell subtype analysis and visualization
# Description: Reproducible pipeline for B cell subtype scoring,
#              trajectory analysis, and visualization (Figure 5)
# Author: Bruce Zhou
# Date: 2025-10-29
##############################################################

rm(list = ls())
gc()

##############################################################
# 1. Load required packages
##############################################################
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggsci)
library(ggpubr)
library(latex2exp)
library(rstatix)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSEABase)
library(enrichplot)
library(GSVA)
library(msigdbr)
library(AUCell)
plan("multicore", workers = 20)

##############################################################
# 2. Load data
##############################################################
setwd("your_project_directory")  # <-- 请修改为项目路径

sce <- readRDS("rds/B_more.rds")
seu <- readRDS("rds/GSE16026_B.rds")

##############################################################
# 3. Label transfer: Map subtype annotations from reference to query
##############################################################
B_anchors <- FindTransferAnchors(
  reference = sce,
  query = seu,
  dims = 1:30,
  reference.reduction = "harmony"
)

predictions <- TransferData(
  anchorset = B_anchors,
  refdata = sce$B_type,
  dims = 1:30
)

seu <- AddMetaData(seu, metadata = predictions)
seu$B_type <- seu$predicted.id

# Visualization check
DimPlot(seu, group.by = c("B_type", "celltype"), reduction = "tsne")

##############################################################
# 4. Factor order and DEG-based signature scoring
##############################################################
seu$B_type <- factor(seu$B_type,
  levels = c("Naive_B", "CD27+_Bm", "CD27-_Bm", "Plasma",
             "EGR1+_B", "FCRL4+_Bm", "pre_GC", "GC_B")
)
Idents(seu) <- "B_type"

# Find DEGs in EGR1+ B cells
deg <- FindMarkers(
  sce,
  ident.1 = "EGR1+_B",
  min.pct = 0.5,
  logfc.threshold = 0.5,
  only.pos = TRUE
)
egr1_bm_gene <- rownames(top_n(deg, n = 100, wt = avg_log2FC))

##############################################################
# 5. AUCell scoring for EGR1+ B signature
##############################################################
geneSets <- GeneSet(egr1_bm_gene, setName = "EGR1_B_signature")

cells_rankings <- AUCell_buildRankings(
  seu@assays$RNA@data,
  nCores = 20,
  plotStats = TRUE,
  splitByBlocks = TRUE
)

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
seu$EGR1_B_score <- as.numeric(getAUC(cells_AUC)["EGR1_B_signature", ])

##############################################################
# 6. Visualization: violin plot across pathologic stages
##############################################################
seu <- subset(seu, tumor_type %in% c("ESCC"))
seu$pathologic_stage <- factor(seu$pathologic_stage, levels = c("I", "II", "III"))

data <- seu@meta.data[, c("pathologic_stage", "EGR1_B_score")]

# Violin + boxplot + significance
p1 <- ggplot(data, aes(pathologic_stage, EGR1_B_score)) +
  geom_violin(aes(fill = pathologic_stage), color = "white") +
  geom_boxplot(fill = "#a6a7ac", color = "#a6a7ac", width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("#A6CEE3FF", "grey", "#E93E3F")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("Pathologic stage") +
  ylab("EGR1+ B score") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  geom_signif(
    comparisons = list(c("I", "II"), c("II", "III"), c("I", "III")),
    map_signif_level = TRUE,
    tip_length = 0,
    y_position = c(0.7, 0.75, 0.8),
    size = 0.5,
    textsize = 5,
    color = "black",
    test = "t.test"
  )

# Save high-resolution figure
tiff("figure5H_EGR1_B_score.tiff", units = "in", width = 6, height = 5.5, res = 600, compression = "lzw")
print(p1)
dev.off()

##############################################################
# 7. Trajectory analysis (Monocle)
##############################################################
# load local package (if not installed)
devtools::load_all("D:/R/sc_env/package/monocle/")

HSMM <- readRDS("rds/gse16_monocle_egr1_B.rds")

# Plot trajectory
a1 <- plot_cell_trajectory(
  HSMM,
  markers_linear = TRUE,
  cell_size = 0.7,
  show_branch_points = FALSE,
  color_by = "pathologic_stage"
) +
  scale_color_manual(values = c("#A6CEE3FF", "grey", "#E93E3F"))

tiff("figure5F_trajectory.tiff", units = "in", width = 6, height = 5, res = 600, compression = "lzw")
print(a1)
dev.off()

##############################################################
# 8. Pseudotime density visualization
##############################################################
df <- pData(HSMM)
color_panel <- c("I" = "#A6CEE3FF", "II" = "#B2DF8AFF", "III" = "#E93E3F")

a5 <- ggplot(df, aes(Pseudotime, color = pathologic_stage, fill = pathologic_stage)) +
  geom_density(bw = 0.85, size = 0.8, alpha = 0.2) +
  theme_classic2() +
  scale_fill_manual(values = color_panel) +
  scale_color_manual(values = color_panel)

tiff("figure5G_pseudotime_density.tiff", units = "in", width = 6, height = 4.5, res = 600, compression = "lzw")
print(a5)
dev.off()
