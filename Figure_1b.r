# CRM_2025_Figure1B
# Author: Bruce Zhou
# Purpose: Generate Figure 1B for Cell Reports Medicine submission
# Notes: High-resolution TIFF output (300 dpi), modular workflow

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggsci)
library(ggpubr)
library(latex2exp)
library(rstatix)
library(RColorBrewer)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(clusterProfiler)
library(survival)
library(survminer)
library(IOBR)
library(data.table)
library(reshape2)

# Set working directory
setwd("work_path")

###########################
# 1. Read IBD cohort data 
counts_ibd <- read.csv("ibd.csv", header = TRUE)
counts_ibd <- counts_ibd[!duplicated(counts_ibd$SYMBOL), ] # remove duplicates
counts_ibd <- counts_ibd[complete.cases(counts_ibd), ]
rownames(counts_ibd) <- counts_ibd$SYMBOL
counts_ibd <- counts_ibd[, -1]  # remove SYMBOL column if needed

###########################
# 2. GSVA Analysis        
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

# Convert GSVA result to dataframe
gsva_score <- as.data.frame(t(gsva_results))

###########################
# 3. Combine clinical data 
clinical_info <- read.csv("filename_anon.csv", header = TRUE)
ibd_data <- cbind(clinical_info, gsva_score)

###########################
# 4. Kaplan-Meier analysis 
km_data <- ibd_data[ibd_data$CANCER_TYPE %in% c("mUC"), ]
km_data$status <- ifelse(
  km_data$GOBP_CELLULAR_SENESCENCE > median(km_data$GOBP_CELLULAR_SENESCENCE, na.rm = TRUE),
  "High", "Low"
)

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

# Save KM plot as high-res TIFF
tiff("ibd_cs_os.tiff", units = "in", width = 6, height = 7, res = 300, compression = "lzw")
print(km_plot)
dev.off()

###########################
# 5. CIBERSORT analysis   
counts_ibd_filtered <- counts_ibd[, colnames(counts_ibd) %in% rownames(km_data)]
cibersort_res <- deconvo_tme(eset = counts_ibd_filtered, method = "cibersort", arrays = TRUE, perm = 100)
rownames(cibersort_res) <- cibersort_res$ID
km_data <- cbind(km_data, cibersort_res)

###########################
# 6. Boxplot of immune composition 
data_barplot <- cibersort_res
colnames(data_barplot)[1:2] <- c("Group", "Sample")
TME_New <- melt(data_barplot)
colnames(TME_New) <- c("Group", "Sample", "Celltype", "Composition")
TME_New$Celltype <- sub("^(.*)_CIBERSORT$", "\\1", TME_New$Celltype)

# Order by median of High group
plot_order <- TME_New %>%
  filter(Group == "High") %>%
  group_by(Celltype) %>%
  summarise(m = median(Composition)) %>%
  arrange(desc(m)) %>%
  pull(Celltype)

TME_New$Celltype <- factor(TME_New$Celltype, levels = plot_order)

# Log-transform composition for plotting
TME_New$Composition_1 <- log10(1 + TME_New$Composition)

# Boxplot theme
mytheme <- theme(
  plot.title = element_text(size = 12, color = "black", hjust = 0.5),
  axis.title = element_text(size = 12, color = "black"),
  axis.text = element_text(size = 8, color = "black"),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1),
  panel.grid = element_blank(),
  legend.position = "top",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 12)
)

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition)) +
  labs(y = "Cell composition", x = NULL, title = "TME Cell composition") +
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#E93E3F", "#A6CEE3FF")) +
  theme_classic() + mytheme +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = TRUE)

# Save boxplot as high-res TIFF
tiff("filename.tiff", units = "in", width = 8, height = 5, res = 300, compression = "lzw")
print(box_TME)
dev.off()
