################################################################################
# Figure 4 Analysis Pipeline
# Study: ESCC & nICT scRNA-seq
# Author: Bruce Zhou
# Last Update: 2025-10-29
# Description:
#   This script generates Figure 4 for the manuscript, including:
#     - Myeloid cell subtype visualization
#     - DEG & GSEA analyses
#     - ssGSEA enrichment (M1/M2/T-cell regulation)
#     - CellChat ligand–receptor network inference
#     - NicheNet ligand–target analysis
# Output:
#   High-resolution figures (.tiff/.pdf) 
################################################################################

rm(list = ls())
gc()

############################################
# Load libraries
############################################

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

############################################
# Environment setup
############################################
base.dir <- "D:/R/sc_env/escc&ncit/"
setwd(base.dir)
dir.create("results/Figure4", showWarnings = FALSE)
setwd("results/Figure4")

################################################################################
# Figure 4A – Myeloid cell subtype visualization
################################################################################
sce <- readRDS(file.path(base.dir, "step5_Macro/Myeloid.rds"))
mylevel <- c('CD1C_DC','CLEC9A_DC','LAMP3_DC','pDC','Neu',
             'CD14_Mono','FCGR3A_Mono',
             'VCAN_Macro','APOE_Macro','GPX3_Macro','TREM2_TAM')
sce$myeloid_type <- factor(sce$myeloid_type, levels = mylevel)

palette <- colorRampPalette(brewer.pal(12, "Paired"))(length(mylevel))
p1 <- DimPlot(sce, group.by = "myeloid_type", label = TRUE, reduction = "tsne") +
  scale_colour_manual(values = palette) +
  theme_bw()

ggsave("Figure4A_myeloid_tsne.pdf", p1, width = 10.5, height = 9, dpi = 300)

################################################################################
# Figure 4B – DEG + GSVA Analysis (R vs NR)
################################################################################
Idents(sce) <- "tumor_type"
deg <- FindMarkers(sce, ident.1 = "ESCC_NR", ident.2 = "ESCC_R",
                   min.pct = 0.01, logfc.threshold = 0.01)

genelist <- deg$avg_log2FC
names(genelist) <- toupper(rownames(deg))
genelist <- sort(genelist, decreasing = TRUE)

# Import curated gene sets
gmt_files <- c("COATES_MACROPHAGE_M1_VS_M2_DN.v2023.2.Hs.gmt",
               "GOBP_NEGATIVE_REGULATION_OF_T_CELL_CYTOKINE_PRODUCTION.v2023.2.Hs.gmt",
               "GOBP_NEGATIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY.v2023.2.Hs.gmt")
gsets <- do.call(rbind, lapply(gmt_files, read.gmt))

egmt <- GSEA(genelist, TERM2GENE = gsets,
             minGSSize = 0.1, maxGSSize = 100000,
             pvalueCutoff = 0.3, pAdjustMethod = "BH")

p2 <- gseaNb(egmt, geneSetID = unique(gsets$term)[1:3],
             subPlot = 2, termWidth = 35,
             legend.position = c(0.74, 0.8),
             addPval = TRUE)
ggsave("Figure4B_GSVA.pdf", p2, width = 7, height = 7, dpi = 300)

################################################################################
# Figure 4C – ssGSEA (M1/M2 and T-cell suppression)
################################################################################
exp <- as.matrix(sce@assays$RNA@counts)
geneSets <- split(gsets$gene, gsets$term)
GSEA_hall <- gsva(expr = exp, gset.idx.list = geneSets, method = "ssgsea",
                  kcdf = "Poisson", mx.diff = TRUE, abs.ranking = FALSE, parallel.sz = 24)
GSEA_score <- t(GSEA_hall)

sce <- AddMetaData(sce, GSEA_score)
sce$M2_score <- sce$COATES_MACROPHAGE_M1_VS_M2_DN
sce$NEGATIVE_REGULATION_OF_T <- sce$GOBP_NEGATIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY

features_plot <- c("M2_score", "NEGATIVE_REGULATION_OF_T", "IL10", "SPP1", "TREM2", "TGFB1", "CD163")
p3 <- DotPlot(sce, features = features_plot, group.by = "myeloid_type", scale = TRUE) +
  scale_color_distiller(palette = "RdYlBu") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
ggsave("Figure4C_M2_signature.pdf", p3, width = 6, height = 6, dpi = 300)

################################################################################
# Figure 4D–E – CellChat analysis
################################################################################
cellchat_dir <- "results/Figure4/cellchat"
dir.create(cellchat_dir, showWarnings = FALSE)
setwd(cellchat_dir)

sce <- readRDS(file.path(base.dir, "step5_Macro/Myeloid.rds"))
sce$celltype_chat <- sce$celltype

# Example: comparing ESCC_NR vs ESCC_R
R_chat <- readRDS("../ESCC_NCIT_R_cellchat2_TM.rds")
NR_chat <- readRDS("../ESCC_NCIT_NR_cellchat2_TM.rds")

object.list <- list(NCIT_R = R_chat, NCIT_NR = NR_chat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

tiff("Figure4E_CellChat_circle.tiff", width = 12, height = 6, units = "in", res = 300, compression = "lzw")
par(mfrow = c(1, 2), xpd = TRUE)
for (i in seq_along(object.list)) {
  count <- object.list[[i]]@net$count
  netVisual_circle(count, weight.scale = TRUE, label.edge = TRUE,
                   edge.weight.max = 20, edge.width.max = 12,
                   title.name = paste0("Interactions - ", names(object.list)[i]))
}
dev.off()

################################################################################
# Figure 4F – NicheNet ligand-target inference
################################################################################
nichenet_dir <- "results/Figure4/nichenet"
dir.create(nichenet_dir, showWarnings = FALSE)
setwd(nichenet_dir)

ligand_target_matrix <- readRDS("D:/R/sc_env/package/nichenete/ligand_target_matrix.rds")
lr_network <- readRDS("D:/R/sc_env/package/nichenete/lr_network.rds")
weighted_networks <- readRDS("D:/R/sc_env/package/nichenete/weighted_networks.rds")

sce <- readRDS(file.path(base.dir, "EGR1_B_TAM.rds"))
sce$celltype_chat <- sce$celltype

nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj = sce,
  receiver = c('EGR1+_B','TREM2_TAM','VCAN_Macro','APOE_Macro'),
  condition_colname = "tumor_type",
  condition_oi = "ESCC_NR",
  condition_reference = "ESCC_R",
  sender = c("EGR1+_B"),
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  weighted_networks = weighted_networks
)

# Ligand-target heatmap
tiff("Figure4F_NicheNet_targets.tiff", width = 12, height = 6, units = "in", res = 300, compression = "lzw")
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke", high = "#E93E3F")
dev.off()
ggsave("Figure4F_NicheNet_targets.pdf", nichenet_output$ligand_target_heatmap, width = 12, height = 6, dpi = 300)

################################################################################
# End of script
################################################################################
