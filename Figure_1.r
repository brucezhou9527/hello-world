# =========================================
# CRM_2025_Figure1
# Author: Bruce Zhou
# R4.2.1
# Purpose: Generate Figure 1 (1B-1F) for Cell Reports Medicine submission
# Notes: High-resolution TIFF/PDF output (300 dpi), modular workflow
# =========================================

# -------------------------------
# Load required libraries
# -------------------------------
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(ggpubr)
library(latex2exp)
library(rstatix)
library(RColorBrewer)
library(tibble)
library(tidyr)
library(stringr)
library(reshape2)
library(GSEABase)
library(GSVA)
library(clusterProfiler)
library(survival)
library(survminer)
library(IOBR)
library(scRNAtoolVis)
library(pheatmap)
library(data.table)

# -------------------------------
# Set working directory
# -------------------------------
setwd("work_path")
getwd()

# ========================================
# 1. Figure 1B: IBD cohort GSVA + KM + CIBERSORT
# ========================================

# 1.1 Read IBD counts
counts_ibd <- read.csv("ibd.csv", header = TRUE) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  drop_na()
rownames(counts_ibd) <- counts_ibd$SYMBOL
counts_ibd <- counts_ibd[, -1]

# 1.2 GSVA analysis
gobp_genes <- read.gmt("GOBP_CELLULAR_SENESCENCE.gmt")
geneSets <- split(gobp_genes$gene, gobp_genes$term)

gsva_results <- gsva(
  expr = as.matrix(counts_ibd),
  gset.idx.list = geneSets,
  method = "ssgsea",
  kcdf = "Gaussian",
  mx.diff = TRUE,
  parallel.sz = 24,
  abs.ranking = FALSE
)
gsva_score <- as.data.frame(t(gsva_results))

# 1.3 Combine with clinical data
clinical_info <- read.csv("filename_anon.csv", header = TRUE)
ibd_data <- cbind(clinical_info, gsva_score)

# 1.4 Kaplan-Meier analysis
km_data <- ibd_data %>% filter(CANCER_TYPE == "mUC")
km_data$status <- ifelse(km_data$GOBP_CELLULAR_SENESCENCE > median(km_data$GOBP_CELLULAR_SENESCENCE, na.rm = TRUE),
                         "High", "Low")

fit <- survfit(Surv(OS_MONTHS, OS_CNSR_REV) ~ status, data = km_data)
km_plot <- ggsurvplot(
  fit,
  surv.scale = "percent",
  break.time.by = 6,
  size = 1.5,
  xlim = c(0, 24),
  risk.table = TRUE,
  tables.height = 0.3,
  tables.y.text.col = FALSE,
  linetype = c(1, 1),
  legend.labs = c("High", "Low"),
  font.legend = c(8, "plain", "black"),
  palette = c("#E93E3F", "#A6CEE3FF"),
  ylab = "Overall survival probability",
  xlab = "Months after ICB",
  pval = TRUE,
  conf.int = FALSE
)

tiff("figure_1B_KM.tiff", units = "in", width = 6, height = 7, res = 300, compression = "lzw")
print(km_plot)
dev.off()

# 1.5 CIBERSORT immune deconvolution
counts_ibd_filtered <- counts_ibd[, colnames(counts_ibd) %in% rownames(km_data)]
cibersort_res <- deconvo_tme(eset = counts_ibd_filtered, method = "cibersort", arrays = TRUE, perm = 100)
rownames(cibersort_res) <- cibersort_res$ID
km_data <- cbind(km_data, cibersort_res)

# 1.6 Boxplot of immune composition
data_barplot <- cibersort_res
colnames(data_barplot)[1:2] <- c("Group", "Sample")
TME_New <- melt(data_barplot)
colnames(TME_New) <- c("Group", "Sample", "Celltype", "Composition")
TME_New$Celltype <- sub("^(.*)_CIBERSORT$", "\\1", TME_New$Celltype)

plot_order <- TME_New %>%
  filter(Group == "High") %>%
  group_by(Celltype) %>%
  summarise(m = median(Composition)) %>%
  arrange(desc(m)) %>%
  pull(Celltype)
TME_New$Celltype <- factor(TME_New$Celltype, levels = plot_order)

mytheme <- theme(
  plot.title = element_text(size = 12, color = "black", hjust = 0.5),
  axis.title = element_text(size = 12, color = "black"),
  axis.text = element_text(size = 8, color = "black"),
  panel.grid = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1),
  legend.position = "top",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12)
)

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition)) +
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#E93E3F", "#A6CEE3FF")) +
  labs(y = "Cell composition", x = NULL, title = "TME Cell composition") +
  theme_classic() + mytheme +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = TRUE)

tiff("figure_1B_TME_boxplot.tiff", units = "in", width = 8, height = 5, res = 300, compression = "lzw")
print(box_TME)
dev.off()

# ========================================
# 2. Figure 1C-F: scRNA-seq B cell analysis
# ========================================

# 2.1 Load scRNA-seq data
sce <- readRDS("../../step4_B/B_more.rds")
sce$B_type <- factor(sce$B_type, levels = c("Naive_B","CD27+_Bm","CD27-_Bm","Plasma",
                                            "EGR1+_B","FCRL4+_Bm","pre_GC","GC_B"))
Idents(sce) <- "B_type"

mycolors <- c("#4F9F3B","#FDAC4F","#A6CEE3","#D1AAB7","#E93E3F","#A99099","#B15928","#72B29C")

# 2.2 tSNE plot (Figure 1C)
df_tsne <- sce@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cluster = sce$B_type)

cluster_centers <- df_tsne %>%
  group_by(cluster) %>%
  summarise(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2))

# Optional manual adjustment to avoid overlap
cluster_centers$label_x <- cluster_centers$tSNE_1 + c(14, -21, 8, -5, 15, -5, 2, -8)
cluster_centers$label_y <- cluster_centers$tSNE_2 + c(16, 5, -9, 18, -9, -2, 17, 15)

tsne_plot <- ggplot(df_tsne, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
  geom_point(size = 0.3) +
  stat_ellipse(aes(fill = cluster), geom = "polygon", linetype = 2, alpha = 0.1, show.legend = FALSE, level = 0.9) +
  geom_text(data = cluster_centers, aes(x = label_x, y = label_y, label = cluster), fontface = "bold") +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1)) +
  labs(title = "tSNE of B cell subtypes")

tiff("figure_1C_TSNE_B.tiff", units = "in", width = 7, height = 6, res = 300, compression = "lzw")
print(tsne_plot)
dev.off()

# 2.3 Stacked bar plot (Figure 1D)
stackData <- data.frame(sample = sce$treated_patholage, Celltype = sce$B_type)
df_bar <- stackData %>%
  group_by(sample, Celltype) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(Celltype = factor(Celltype, levels = c("Naive_B","CD27+_Bm","CD27-_Bm","Plasma",
                                                "pre_GC","GC_B","FCRL4+_Bm","EGR1+_B")),
         sample = factor(sample, levels = c("sd","mpr","pcr")))

mycolors_bar <- c('#4F9F3B','#FDAC4F','#A6CEE3','#D1AAB7','#B15928','#72B29C','#A99099','#E93E3F')

stacked_plot <- ggplot(df_bar, aes(x = sample, y = count, fill = Celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = mycolors_bar) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9)) +
  ylab("Percentage") +
  xlab("")

tiff("figure_1D_stacked_B_cells.tiff", units = "in", width = 4, height = 5, res = 300, compression = "lzw")
print(stacked_plot)
dev.off()

# 2.4 Ro/e plots (Figure 1D)
sce$patient_id <- str_extract(sce$orig.ident, "s\\d+")
metainfo <- data.frame(sce@meta.data[, c("patient_id", "B_type", "tumor_type")])
metainfo$UniqueCell_ID <- rownames(metainfo)
metainfo$pathology <- factor(metainfo$pathology, levels = c("Pre", "Post_PCR", "Post_MPR", "Post_SD"))

source("distribution_Roe.R")
# By tumor_type
distribution_Roe(meta_data = metainfo, celltype_column = "B_type",
                 condition_column = "tumor_type", add_label = "sign",
                 tile_color = "grey", tile_fill = "D")
# By pathology
distribution_Roe(meta_data = metainfo, celltype_column = "B_type",
                 condition_column = "pathology", add_label = "number",
                 condition_level = c("Pre", "Post_PCR", "Post_MPR", "Post_SD"), tile_color = NA)

ggsave("figure_1D_RoE.pdf", width = 6.5, height = 8.5, units = "in", dpi = 300)

# 2.5 Heatmap (Figure 1E)
DefaultAssay(sce) <- "RNA"
cluster <- sce$B_type
my_order <- levels(cluster)
IDX <- lapply(my_order, function(x) which(cluster == x))

geneList_selected <- c('GENE1','GENE2','GENE3') # <-- replace with actual genes
idx <- match(geneList_selected, rownames(sce))

expression <- sapply(IDX, function(cells) {
  rowMeans(as.matrix(sce@assays$RNA@data[idx, cells, drop=FALSE]), na.rm = TRUE)
})
rownames(expression) <- geneList_selected
colnames(expression) <- my_order

heat_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
p_heatmap <- pheatmap(expression, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
                      show_colnames = TRUE, show_rownames = TRUE, color = heat_colors)

tiff('figure_1E_heatmap_B_cells.tiff', units="in", width=3.5, height=9, res=300, compression='lzw')
print(p_heatmap)
dev.off()

# 2.6 FeaturePlot (Figure 1F)
T_cell_genes <- c("CD3D", "CD8A", "CD4")
B_cell_genes <- c("CD79A", "MS4A1", "CD19")
plasma_genes <- c("JCHAIN", "MZB1")
mono_genes <- c("CD14", "MS4A7")
NK_genes <- c("NKG7", "GZMA")
DC_genes <- c("FCER1A", "CST3")
Mast_cell_genes <- c("CPA3")
cycling_genes <- c('MKI67')

gene_feature <- c(T_cell_genes, B_cell_genes, plasma_genes, mono_genes, NK_genes, DC_genes, Mast_cell_genes, cycling_genes)

gc()
tiff('figure_1F_featureplot_B_cells.tiff', units="in", width=8, height=8, res=300, compression='lzw')
FeatureCornerAxes(object = sce, reduction = 'tsne', groupFacet = NULL,
                  relLength = 0.4, relDist = 0.2, pSize = 0.01, aspect.ratio = NULL,
                  nLayout = 4, cornerTextSize = 2, cornerVariable = 1,
                  minExp = 0, maxExp = 4, features = gene_feature)
dev.off()
