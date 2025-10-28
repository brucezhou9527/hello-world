###############################################################
# Figure 3: SCENIC & DoRothEA Transcriptional Activity Analysis
# Author: Bruce Zhou
# Date: 2025-10
# Output:
#   Figure3a1: SCENIC Binary Heatmap
#   Figure3a2: DoRothEA Top TF Heatmap
#   Figure3b: EGR1 FeaturePlot

###############################################################

###########################
# 0. Environment Setup
###########################
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(SCENIC)
library(AUCell)
library(harmony)
library(decoupleR)
library(dorothea)
library(viridis)
library(scCustomize)
library(OmnipathR)
plan("multicore", workers = 20)
options(future.globals.maxSize = 1000000 * 1024^2)

# Working directory
setwd("work_path")
sce <- readRDS("B_more.rds")
scenicOptions <- readRDS("int/scenicOptions.rds")

dir.create("figure_3A", showWarnings = FALSE)

###############################################################
# 1. SCENIC Analysis
###############################################################
## --- Expression and metadata preparation ---
exprMat <- as.matrix(sce@assays$RNA@data)
cellInfo <- data.frame(sce@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="B_type")] <- "celltype"
cellInfo <- cellInfo[, c("sample","cluster","celltype")]

## --- Load SCENIC results ---
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

## --- Regulon visualization ---
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

## Regulon activity by cell type
regulonActivity_byCellType <- sapply(
  split(rownames(cellInfo), cellInfo$celltype),
  function(cells) rowMeans(getAUC(regulonAUC)[, cells])
)
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))

## --- Heatmap (Continuous AUC) ---
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:70,], name="Regulon activity")

## --- Binary heatmap (for Figure 3a1) ---
minPerc <- 0.2
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[rownames(cellInfo) %in% colnames(binaryRegulonActivity), , drop=FALSE]

regulonActivity_byCellType_Binarized <- sapply(
  split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype),
  function(cells) rowMeans(binaryRegulonActivity[, cells, drop=FALSE])
)

binaryActPerc_subset <- regulonActivity_byCellType_Binarized[
  rowSums(regulonActivity_byCellType_Binarized > minPerc) > 0, 
]

# === Figure 3a1: SCENIC binary heatmap ===
tiff("figure_3A/Figure3a1_SCENIC_binary.tiff", units="in", width=6, height=8, res=600, compression="lzw")
pheatmap(
  binaryActPerc_subset[1:40,],
  name="Regulon activity (%)",
  cluster_rows=TRUE,
  cluster_cols=FALSE,
  show_colnames=TRUE,
  scale="row",
  color=colorRampPalette(c("#0044BB","white","#FF7744"))(100)
)
dev.off()

###############################################################
# 2. DoRothEA Analysis
###############################################################
## --- Load or prepare network ---
net <- readRDS("D:/R/sc_env/package/dorothea/net.rds")

## --- ULM-based TF activity inference ---
mat <- as.matrix(sce@assays$RNA@data)
plan("multisession", workers=16)
acts <- run_ulm(mat, net, minsize=4)

sce[["tfsulm"]] <- acts %>%
  pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
  column_to_rownames("source") %>%
  Seurat::CreateAssayObject()

DefaultAssay(sce) <- "tfsulm"
sce <- ScaleData(sce)
sce@assays$tfsulm@data <- sce@assays$tfsulm@scale.data

## --- Aggregate by cluster ---
n_tfs <- 100
df <- t(as.matrix(sce@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(sce)) %>%
  pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
  group_by(cluster, source) %>%
  summarise(mean=mean(score), .groups="drop")

tfs <- df %>%
  group_by(source) %>%
  summarise(std=sd(mean), .groups="drop") %>%
  arrange(desc(abs(std))) %>%
  head(n_tfs) %>%
  pull(source)

top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols="cluster", names_from="source", values_from="mean") %>%
  column_to_rownames("cluster") %>%
  as.matrix()

# === Figure 3a2: DoRothEA heatmap ===
my_breaks <- c(seq(-3, 0, length.out=51), seq(0.05, 3, length.out=49))
pheatmap(
  t(top_acts_mat),
  border_color=NA,
  color=colorRampPalette(c("navy","white","firebrick3"))(100),
  cluster_rows=TRUE,
  cluster_cols=FALSE,
  breaks=my_breaks
)
ggsave("figure_3A/Figure3a2_DoRothEA_top100.pdf", height=10, width=4)

###############################################################
# 3. FeaturePlot Visualization (Figure 3b)
###############################################################
DefaultAssay(sce) <- "tfsulm"
genes <- c("EGR1","EGR1_11g","EGR1_extended_204g")

plots <- lapply(genes, function(gene){
  FeaturePlot_scCustom(
    seurat_object=sce,
    reduction='tsne',
    colors_use=colorRampPalette(c("#3288BD","white","#D53E4F"))(50),
    features=gene
  ) + NoAxes()
})

p <- wrap_plots(plots, ncol=1)
ggsave("figure_3A/Figure3b_EGR1_FeaturePlot.pdf", plot=p, width=5, height=12, dpi=600)

