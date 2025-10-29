################################################################################
# Figure S2 | B cell trajectory and differentiation potential analysis
# Author: Bruce Zhou
# Description: CytoTRACE2 and Monocle pseudotime analysis for B cells
# Date: 2025-10-29
################################################################################

#=============================#
#   Environment Preparation   #
#=============================#
rm(list = ls())
gc()

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

# Set working directory
setwd("workpath")

plan("multicore", workers = 20)

################################################################################
# Figure S2A: Slingshot trajectory (reference in Figure2.r)
################################################################################
# -- See "Figure2.r" for detailed Slingshot analysis workflow --

################################################################################
# Figure S2B: CytoTRACE2 analysis of B cell differentiation potential
################################################################################
library(CytoTRACE2)
library(export)
library(reticulate)
library(SeuratData)

# Load Seurat object
sce <- readRDS("B_rds")

# Run CytoTRACE2
cytotrace2_result_sce <- cytotrace2(
  sce,
  is_seurat = TRUE,
  slot_type = "counts",
  species = "human",
  ncores = 20,
  seed = 1234
)

# Prepare phenotype annotation
annotation <- data.frame(phenotype = sce@meta.data$B_type) %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  tibble::column_to_rownames("cell_id")

# Generate CytoTRACE2 plots
plots <- plotData(
  cytotrace2_result = cytotrace2_result_sce,
  annotation = annotation,
  is_seurat = TRUE
)

# Extract individual plots
p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$CytoTRACE2_Boxplot_byPheno

# Preview combined layout
library(patchwork)
(p1 + p2 + p3 + p4) + plot_layout(ncol = 2)

#=== Export Figure S2B ===#
tiff("figure_S2B_cytotrace.tiff", units = "in", width = 10, height = 8,
     res = 300, compression = "lzw")
p4
dev.off()


################################################################################
# Figure S2C-D: Monocle pseudotime trajectory analysis
################################################################################
# Load Monocle (from local source if needed)
devtools::load_all("rpackages/monocle/")

# Load Seurat data
sce <- readRDS("B_rds")
seu <- subset(sce, B_type %in% c("Naive_B", "CD27+_Bm", "EGR1+_B"))

# Convert Seurat → Monocle CellDataSet
data <- as(as.matrix(seu@assays$RNA@counts), "sparseMatrix")
gene_ann <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
sample_ann <- seu@meta.data

fd <- new("AnnotatedDataFrame", data = gene_ann)
pd <- new("AnnotatedDataFrame", data = sample_ann)

HSMM <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.01,
  expressionFamily = negbinomial.size()
)

# Estimate size factors & dispersions
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

# Filter low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 1)
expressed_genes <- rownames(subset(fData(HSMM), num_cells_expressed >= 5))

# Identify ordering genes
diff_test_res <- differentialGeneTest(
  HSMM[expressed_genes, ],
  fullModelFormulaStr = "~ seurat_clusters",
  cores = 16, verbose = TRUE
)
ordering_genes <- rownames(subset(diff_test_res, qval < 0.01))

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

# Dimensionality reduction & trajectory ordering
HSMM <- reduceDimension(
  HSMM,
  max_components = 3,
  num_dim = 10,
  method = "DDRTree"
)
HSMM <- orderCells(HSMM)

#============================#
#    Visualization Section   #
#============================#
colour <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080",
            "#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE",
            "#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0",
            "#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22",
            "#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347",
            "#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

mycol <- c("#0073C2FF", "#EFC000FF", "#7AA6DCFF", "#E64B35FF", "#00A087FF", "black")

# Plot trajectories
a1 <- plot_cell_trajectory(HSMM, cell_size = 0.5, color_by = "B_type") +
  scale_color_manual(values = colour)
a2 <- plot_cell_trajectory(HSMM, cell_size = 0.5, color_by = "State") +
  scale_color_manual(values = colour)
a3 <- plot_cell_trajectory(HSMM, cell_size = 0.3, show_branch_points = FALSE, color_by = "Pseudotime")
a4 <- plot_cell_trajectory(HSMM, cell_size = 0.3, show_branch_points = FALSE, color_by = "B_main") +
  scale_color_manual(values = mycol)

p1 <- a1 + a3

#=== Export Figure S2C ===#
tiff("figure_S2C_trajectory_B.tiff", units = "in", width = 5, height = 6,
     res = 300, compression = "lzw")
p1
dev.off()

# Save monocle object
saveRDS(HSMM, file = "monocle_B.rds")

#============================#
#    Pseudotime Heatmap      #
#============================#
diff_test_res <- differentialGeneTest(
  HSMM[expressed_genes, ],
  fullModelFormulaStr = "~sm.ns(Pseudotime)",
  cores = 8
)

sig_gene_names <- rownames(subset(diff_test_res, qval < 0.05))

# Plot pseudotime heatmap
pseudoplot <- plot_pseudotime_heatmap(
  HSMM[sig_gene_names, ],
  num_clusters = 4,
  cores = 8,
  show_rownames = FALSE,
  return_heatmap = TRUE
)

#=== Export Figure S2D ===#
tiff("figure_S2D_pseudotime_heatmap.tiff", units = "in", width = 8, height = 6,
     res = 300, compression = "lzw")
pseudoplot
dev.off()
