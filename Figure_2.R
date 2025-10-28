# =========================================
# Figure 2: B cell trajectory & EGR1+ analysis
# Author: Bruce Zhou
# Purpose: Generate Figure 2 panels (2A-2F)
# High-resolution output (300 dpi)
# =========================================

# -------------------------------
# Load required libraries
# -------------------------------
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
library(slingshot)
library(SingleCellExperiment)
library(tradeSeq)
library(scales)
library(patchwork)
library(cowplot)
library(ggrepel)
library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(limma)
library(DESeq2)
library(showtext)
library(ggthemes)
library(ggprism)

# -------------------------------
# Setup parallelization
# -------------------------------
plan("multicore", workers = 20)
options(future.globals.maxSize = 1000000 * 1024^2)

# -------------------------------
# Set working directory
# -------------------------------
setwd("your_working_directory")

# -------------------------------
# Load B cell Seurat object
# -------------------------------
sce <- readRDS('../../step4_B/B_more.rds')
Idents(sce) <- 'B_type'

# -------------------------------
# Figure 2A: Slingshot trajectory
# -------------------------------
seuratObject <- as.SingleCellExperiment(sce)
seuratObject <- slingshot(seuratObject, 
                          clusterLabels = 'B_type', 
                          reducedDim = 'TSNE',
                          start.clus= 'Naive_B',
                          stretch=0, 
                          approx_points = 100)

# Compute pseudotime0 as max across all lineages
seuratObject$pseudotime1 <- ifelse(is.na(seuratObject$slingPseudotime_1), 0, seuratObject$slingPseudotime_1)
seuratObject$pseudotime2 <- ifelse(is.na(seuratObject$slingPseudotime_2), 0, seuratObject$slingPseudotime_2)
seuratObject$pseudotime3 <- ifelse(is.na(seuratObject$slingPseudotime_3), 0, seuratObject$slingPseudotime_3)
seuratObject$pseudotime0 <- pmax(seuratObject$pseudotime1, seuratObject$pseudotime2, seuratObject$pseudotime3)

# t-SNE plot colored by pseudotime
tsne_coords <- as.data.frame(reducedDims(seuratObject)$TSNE)
colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
tsne_coords$pseudotime <- seuratObject$pseudotime0
colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)

ggplot(tsne_coords, aes(x = tSNE_1, y = tSNE_2, color = pseudotime)) +
  geom_point(size = 0.3) +
  scale_color_gradientn(colors = colors) +
  labs(color = "Pseudotime", x = "tSNE_1", y = "tSNE_2") +
  theme_minimal() +
  theme(legend.position = "right")
ggsave('figure_2A_slingshot.pdf', width=6, height=6.5, dpi=300)

# -------------------------------
# Figure 2B: Pseudotime expression trends
# -------------------------------
# Prepare data for plotting
psts <- as.data.frame(slingPseudotime(seuratObject))
psts$cells <- rownames(psts)
psts$conditions <- seuratObject$tumor_type

gene_list <- c('GENE1','GENE2','GENE3') # <-- replace with actual genes
data <- t(as.data.frame(sce@assays$RNA@data[gene_list,]))
AUC_SASP <- sce@meta.data[,"AUC_SASP", drop=FALSE]
data_psts <- cbind(psts, data, AUC_SASP)
data_long <- pivot_longer(data_psts, cols = all_of(gene_list), names_to = "gene", values_to = "value")

# Define lineage mapping
data_long$lineages_1 <- case_when(
  data_long$lineages %in% 'Lineage1' ~ 'CS_Lineage',
  data_long$lineages %in% 'Lineage2' ~ 'Plasma_Lineage',
  data_long$lineages %in% 'Lineage3' ~ 'GC_Lineage',
  TRUE ~ NA_character_
)
data_long <- data_long %>% filter(conditions %in% c('ESCC_R','ESCC_NR'), !is.na(pseudotime))

# Plot smoothed pseudotime expression
P2 <- ggplot(data_long, aes(x = pseudotime, y = value, color = lineages_1)) +
  geom_smooth(aes(fill = lineages_1), se = FALSE) +
  facet_wrap(~gene, scales = "free_y") +
  scale_color_manual(values = c("#E93E3F","#A6CEE3FF","grey")) +
  xlab('Pseudotime') + ylab('Relative expression') +
  theme_minimal() +
  theme(legend.position = "right", strip.text = element_text(size = 12))

tiff('figure_2B.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
P2
dev.off()

# -------------------------------
# Figure 2D: MA plot for EGR1+ B vs Other B
# -------------------------------
sce$B_type <- ifelse(!sce$B_main %in% 'EGR1+_B', 'Other_B','EGR1_B')
av <- AggregateExpression(sce, group.by = c("orig.ident", "B_type"), assays = "RNA", slot = "counts", return.seurat = FALSE)[[1]]
colData <- data.frame(samples = colnames(av)) %>%
  mutate(condition = ifelse(grepl('Other', samples), 'Other_B', 'EGR1_Bm')) %>%
  column_to_rownames('samples')
dds <- DESeqDataSetFromMatrix(countData = av, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "EGR1_Bm", "Other_B"))
DEG_deseq2 <- na.omit(as.data.frame(res))
DEG_deseq2 <- DEG_deseq2 %>%
  mutate(Type = case_when(
    padj > 0.05 ~ "stable",
    abs(log2FoldChange) < 1 ~ "stable",
    log2FoldChange >= 1 ~ "up",
    TRUE ~ "down"
  )) %>% rownames_to_column("Symbol")

# Highlight genes
genes_to_label <-c('GENE1','GENE2','GENE3') # <-- replace with actual genes
DEG_deseq2$label <- ifelse(DEG_deseq2$Symbol %in% genes_to_label, DEG_deseq2$Symbol, "")

ggplot(DEG_deseq2, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = Type), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("up" = "#E93E3F", "down" = "#A6CEE3FF", "stable" = "gray")) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_text_repel(aes(label = label), size=3.5, max.overlaps = Inf, segment.size=0.3) +
  scale_x_log10() +
  labs(title = "MA Plot: EGR1⁺ B cells vs Other B cells", x = "Mean expression (log10 scale)", y = "Log2 Fold Change") +
  theme_bw()
ggsave(file="MA_plot.pdf", width=8, height=6, dpi=300)

# -------------------------------
# Figure 2E: GSEA for EGR1+ B cells
# -------------------------------
deg <- FindMarkers(sce, ident.1 = "EGR1+_B", min.pct = 0.01, logfc.threshold = 0.01)
genelist <- setNames(sort(deg$avg_log2FC, decreasing = TRUE), toupper(rownames(deg)))
g4 <- read.gmt('REACTOME_CELLULAR_SENESCENCE.v2022.1.Hs.gmt')
g6 <- read.gmt('REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP.v2022.1.Hs.gmt')
g7 <- read.gmt('SAUL_SEN_MAYO.v2025.1.Hs.gmt')
c2 <- rbind(g4,g6,g7)
egmt <- GSEA(genelist, TERM2GENE = c2, minGSSize = 10, maxGSSize = 10000, pvalueCutoff = 1, pAdjustMethod = 'BH')
p2 <- gseaNb(egmt, geneSetID = c("REACTOME_CELLULAR_SENESCENCE",
                                 "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                                 "SAUL_SEN_MAYO"), subPlot = 2, termWidth = 35,
             legend.position = c(0.74,0.785), addPval = T, pvalX = 0.01, pvalY = 0.17)
ggsave("figure_2F_EGR1_GSEA.pdf", p2, width=6.5, height=6, units="in", dpi=300)

# -------------------------------
# Figure 2F: Single-cell GSVA & limma barplot
# -------------------------------
m_df <- msigdbr(species = "human", category = "H")
geneSets <- split(m_df$gene_symbol, m_df$gs_name)
exp <- as.matrix(sce@assays$RNA@counts)
meta <- sce@meta.data[, c("orig.ident", "tumor_type", "B_type", "EGR1_state")]

GSVA_hall <- gsva(expr = exp, method = 'ssgsea', gset.idx.list = geneSets, mx.diff = TRUE, kcdf = "Poisson", parallel.sz = 20)
meta$B_type <- ifelse(meta$B_type == 'EGR1+_B', 'EGR1_B', 'Others')
group <- factor(meta$B_type, levels = c('Others','EGR1_B'))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
compare <- makeContrasts(EGR1_B - Others, levels = design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
dat_plot <- data.frame(id = str_replace(rownames(Diff), "HALLMARK_", ""), t = Diff$t)
dat_plot$threshold <- factor(ifelse(dat_plot$t >= 1, "Up", ifelse(dat_plot$t <= -1, "Down", "NoSignifi")),
                             levels = c("Up","Down","NoSignifi"))
dat_plot <- dat_plot %>% arrange(t)
p <- ggplot(dat_plot, aes(x=id, y=t, fill=threshold)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c('Up'= "#E93E3F",'NoSignifi'='#cccccc','Down'='#A6CEE3FF')) +
  geom_hline(yintercept = c(-2,2), color='white', lty='dashed') +
  xlab('') + ylab('t value of ssGSEA score, Others VS EGR1_B') +
  theme_prism(border=T) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
tiff('figure2F_gsea_bar_egr1b.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
p
dev.off()
