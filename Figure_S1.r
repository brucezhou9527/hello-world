##############################################################
# Title: Figure S1 – Single-cell landscape of CD45+ ESCC samples
# Description: Reproducible pipeline to generate Figure S1A–F
# Author: Bruce Zhou
# Date: 2025-10-29

##############################################################
# 1. Load required packages
##############################################################
library(Seurat)
library(dplyr)
library(future)
library(future.apply)
library(msigdbr)
library(clusterProfiler)
library(showtext)
library(ggpubr)
library(RColorBrewer)
library(ggsci)
library(tidyr)
library(tibble)
library(grid)
library(scRNAtoolVis)
library(dittoSeq)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(pheatmap)
library(monocle)
library(Matrix)
library(scran)

plan("multicore", workers = 20)

##############################################################
# 2. Set working directory & define color palettes
##############################################################
setwd("your_project_directory")  # <-- 请修改为你的项目路径
output_dir <- "figures/FigureS1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

my_col <- c(
  "#DC0000FF", "#4DBBD5FF", "#00A087FF", "#F39B7FFF",
  "#3C5488FF", "#8491B4FF", "#6A6599FF", "#FFA319FF"
)

##############################################################
# 3. Load dataset
##############################################################
sce <- readRDS("rds/ESCC_CD45.rds")

##############################################################
# === Figure S1A: UMAP by major cell type, split by tumor type ===
##############################################################
tiff(file.path(output_dir, "FigureS1A_umap_celltype_group.tiff"),
     units = "in", width = 25, height = 7, res = 600, compression = "lzw")

DimPlot(
  sce,
  group.by = "celltype_main",
  split.by = "tumor_type",
  reduction = "umap",
  pt.size = 0.001,
  label = TRUE, repel = TRUE,
  cols = my_col
)
dev.off()

##############################################################
# === Figure S1B: Feature expression patterns (TSNE) ===
##############################################################
T_cell <- c("CD3D", "CD8A", "CD4")
B_cell <- c("CD79A", "MS4A1", "CD19")
plasma <- c("JCHAIN", "MZB1")
mono <- c("CD14", "MS4A7")
NK <- c("NKG7", "GZMA")
DC <- c("FCER1A", "CST3")
Mast_cell <- c("CPA3")
cycling <- c("MKI67")

gene_feature <- c(T_cell, cycling, B_cell, plasma, mono, NK, DC, Mast_cell)

tiff(file.path(output_dir, "FigureS1B_featureplot.tiff"),
     units = "in", width = 8, height = 8, res = 600, compression = "lzw")

FeatureCornerAxes(
  object = sce,
  reduction = "tsne",
  groupFacet = NULL,
  relLength = 0.4,
  relDist = 0.2,
  pSize = 0.01,
  nLayout = 4,
  cornerTextSize = 2,
  cornerVariable = 1,
  minExp = 0, maxExp = 4,
  features = gene_feature
)
dev.off()

##############################################################
# === Figure S1C: DotPlot of canonical markers by cell type ===
##############################################################
list_genes <- list(
  T_cell = c("PTPRC","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","CCR7","LEF1","GZMH","GZMK"),
  B_cell = c("CD79A","CD79B","MS4A1","CD22","CD19"),
  Plasma = c("XBP1","JCHAIN","MZB1"),
  Mono = c("CD14","S100A8","S100A9","FCGR3A","MS4A7"),
  NK = c("NKG7","GZMB","GZMA","CST7"),
  DC = c("FCER1A","CST3","CLEC10A"),
  pDC = c("TCF4","IL3RA","PLD4"),
  Mast_cell = c("KIT","CPA3","TPSAB1","TPSB2")
)

p1 <- DotPlot(
  sce,
  features = list_genes,
  cols = c("grey", "red"),
  group.by = "celltype_main"
) + RotatedAxis() +
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin = margin(b = 3, unit = "mm")),
    axis.line = element_blank()
  ) + labs(x = "", y = "")

# Remove facet borders except top
q <- ggplotGrob(p1)
lg <- linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(0, 0) + 0.2, "npc"),
                gp = gpar(col = "black", lwd = 4))
for (k in grep("strip-t", q$layout$name)) {
  q$grobs[[k]]$grobs[[1]]$children[[1]] <- lg
}

tiff(file.path(output_dir, "FigureS1C_dotplot_marker.tiff"),
     units = "in", width = 14, height = 4, res = 600, compression = "lzw")
grid.draw(q)
dev.off()

##############################################################
# === Figure S1D: Stacked bar plot – cell composition per sample ===
##############################################################
stackData <- data.frame(sample = sce$orig.ident, Celltype = sce$celltype_main)
df1 <- stackData %>%
  group_by(sample) %>%
  table() %>%
  as.data.frame.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Celltype") %>%
  gather(key = "sample", value = "count", -Celltype)

df1$Celltype <- factor(df1$Celltype, levels = c(
  "T cells", "proliferating T", "B cells", "Plasma cells", "proliferating B",
  "NK cells", "Mono/Macro", "DC", "pDC", "Mast cells", "Others"
))

p3 <- ggplot(df1, aes(fill = Celltype, y = count, x = sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = my_col) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) +
  xlab("") + ylab("Proportion")

tiff(file.path(output_dir, "FigureS1D_cell_composition.tiff"),
     units = "in", width = 6, height = 6, res = 600, compression = "lzw")
print(p3)
dev.off()

##############################################################
# === Figure S1E: Heatmap of top markers (dittoSeq) ===
##############################################################
sce <- ScaleData(sce, verbose = TRUE, block.size = 5000)
b.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- b.markers %>%
  group_by(cluster) %>%
  top_n(n = 23, wt = avg_log2FC)

# remove housekeeping genes
hspGenes <- grep("^HSP", top_markers$gene, value = TRUE)
mtGenes <- grep("^MT-", top_markers$gene, value = TRUE)
rplGenes <- grep("^RPL", top_markers$gene, value = TRUE)
rpsGenes <- grep("^RPS", top_markers$gene, value = TRUE)
gene_blocked <- c(hspGenes, mtGenes, rplGenes, rpsGenes)

top_filtered <- top_markers[!top_markers$gene %in% gene_blocked, ]
write.csv(top_filtered, "rds/b_top20.csv", row.names = FALSE)

top30 <- read.csv("rds/b_top20.csv")

mycolors <- c("#4F9F3B","#FDAC4F","#A6CEE3","#D1AAB7",
              "#E93E3F","#A99099","#B15928","#72B29C")
getPalette <- colorRampPalette(brewer.pal(12, "Paired"))
mycolors_1 <- sample(getPalette(4))

p1 <- dittoHeatmap(
  sce,
  top30$gene,
  annot.by = c("B_type", "tumor_type"),
  annot.colors = c(mycolors, mycolors_1),
  heatmap.colors = colorRampPalette(
    c("#1f294e", "#5390b5", "#eaebea", "#d56e5e", "#57121d")
  )(20)
)

tiff(file.path(output_dir, "FigureS1E_marker_heatmap.tiff"),
     units = "in", width = 8, height = 6, res = 600, compression = "lzw")
print(p1)
dev.off()

##############################################################
# === Figure S1F: DotPlot of gene set–based B cell signatures ===
##############################################################
# Load signature gene lists (Excel sheets 1–5)
sig_file <- "signature_gene_b.xlsx"
sig_names <- c("SASP", "Breg", "TLS", "CSR", "Costimulation")
sig_idx <- list()

for (i in 1:5) {
  genes <- read_excel(sig_file, sheet = i, col_names = FALSE)[[1]][-1]
  idx <- match(genes, rownames(sce))
  sig_idx[[i]] <- idx[!is.na(idx)]
}

# Calculate mean expression per signature
expr_mat <- sapply(sig_idx, function(idx) colMeans(as.matrix(sce@assays$RNA@data[idx, ]), na.rm = TRUE))
rownames(expr_mat) <- colnames(sce)
expr_mat <- t(expr_mat)
rownames(expr_mat) <- sig_names

# z-score normalization
expr_scaled <- t(scale(t(expr_mat)))
sce[["RNA"]]@data <- rbind(sce[["RNA"]]@data, expr_scaled)

# DotPlot of B cell signature modules
markers.to.plot <- rev(sig_names)
p1 <- DotPlot(
  sce,
  features = markers.to.plot,
  group.by = "B_type",
  cols = c("grey", "red"),
  dot.scale = 8
) + RotatedAxis() + ggtitle("B cell functional modules")

tiff(file.path(output_dir, "FigureS1F_Bcell_signatures.tiff"),
     units = "in", width = 6, height = 5, res = 600, compression = "lzw")
print(p1)
dev.off()
