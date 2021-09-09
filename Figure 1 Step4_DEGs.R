###############################################################################################
# Pre-processing for data used to generate Figure 1 from Sayed, Kodama, Fan, et al. 2021 
# Total of 46 human AD R47H vs CV samples
# This script is: STEP 4 of 5 - identifying DEGs for different comparisons

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration")
R47H_all_integrated <- readRDS("R47H_all_integrated_Annotation.rds")
setwd("/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration/DEGs")
R47H_all_integrated$TREM2.Sex <- paste(R47H_all_integrated$TREM2, R47H_all_integrated$Sex, sep = "_")

Cluster_EN <- subset(R47H_all_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(R47H_all_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(R47H_all_integrated, idents = "microglia")
Cluster_AST <- subset(R47H_all_integrated, idents = "astrocytes")
Cluster_OL <- subset(R47H_all_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(R47H_all_integrated, idents = "OPCs")
Cluster_EC <- subset(R47H_all_integrated, idents = "endothelial cells")

rm(R47H_all_integrated)

Idents(Cluster_EN) <- "TREM2.Sex"
Idents(Cluster_IN) <- "TREM2.Sex"
Idents(Cluster_MG) <- "TREM2.Sex"
Idents(Cluster_AST) <- "TREM2.Sex"
Idents(Cluster_OL) <- "TREM2.Sex"
Idents(Cluster_OPC) <- "TREM2.Sex"
Idents(Cluster_EC) <- "TREM2.Sex"

EN_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_EN, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_M_R47H_vs_CV_DEGs, "EN_M_R47H_vs_CV_DEGs.csv")
IN_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_IN, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_M_R47H_vs_CV_DEGs, "IN_M_R47H_vs_CV_DEGs.csv")
MG_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_MG, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_M_R47H_vs_CV_DEGs, "MG_M_R47H_vs_CV_DEGs.csv")
AST_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_AST, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_M_R47H_vs_CV_DEGs, "AST_M_R47H_vs_CV_DEGs.csv")
OL_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_OL, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_M_R47H_vs_CV_DEGs, "OL_M_R47H_vs_CV_DEGs.csv")
OPC_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_M_R47H_vs_CV_DEGs, "OPC_M_R47H_vs_CV_DEGs.csv")
EC_M_R47H_vs_CV_DEGs <- FindMarkers(Cluster_EC, ident.1 = "R47H_M", ident.2 = "WT_M", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EC_M_R47H_vs_CV_DEGs, "EC_M_R47H_vs_CV_DEGs.csv")

EN_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_EN, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EN_F_R47H_vs_CV_DEGs, "EN_F_R47H_vs_CV_DEGs.csv")
IN_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_IN, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(IN_F_R47H_vs_CV_DEGs, "IN_F_R47H_vs_CV_DEGs.csv")
MG_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_MG, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(MG_F_R47H_vs_CV_DEGs, "MG_F_R47H_vs_CV_DEGs.csv")
AST_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_AST, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(AST_F_R47H_vs_CV_DEGs, "AST_F_R47H_vs_CV_DEGs.csv")
OL_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_OL, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OL_F_R47H_vs_CV_DEGs, "OL_F_R47H_vs_CV_DEGs.csv")
OPC_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_OPC, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(OPC_F_R47H_vs_CV_DEGs, "OPC_F_R47H_vs_CV_DEGs.csv")
EC_F_R47H_vs_CV_DEGs <- FindMarkers(Cluster_EC, ident.1 = "R47H_F", ident.2 = "WT_F", logfc.threshold = 0.1, min.pct = 0.25, only.pos = F,test.use = "MAST")
write.csv(EC_F_R47H_vs_CV_DEGs, "EC_F_R47H_vs_CV_DEGs.csv")











