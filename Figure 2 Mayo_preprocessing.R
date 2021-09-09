###############################################################################################
# Pre-processing for data used to generate Figure 2 from Sayed, Kodama, Fan, et al. 2021 
# Total of 46 human AD R47H vs CV samples
# This script is to integrate the microglia cells only from the Mayo samples

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====
#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Human_R47H/data_analysis/DF_2ndRound")
Human_R47H_1 <- readRDS(file = "Human_R47H_1_singlets_PCA.rds")
Human_R47H_2 <- readRDS(file = "Human_R47H_2_singlets_PCA.rds")
Human_R47H_3 <- readRDS(file = "Human_R47H_3_singlets_PCA.rds")
Human_R47H_4 <- readRDS(file = "Human_R47H_4_singlets_PCA.rds")
Human_R47H_5 <- readRDS(file = "Human_R47H_5_singlets_PCA.rds")
Human_R47H_6 <- readRDS(file = "Human_R47H_6_singlets_PCA.rds")
Human_R47H_7 <- readRDS(file = "Human_R47H_7_singlets_PCA.rds")
Human_R47H_8 <- readRDS(file = "Human_R47H_8_singlets_PCA.rds")
Human_R47H_9 <- readRDS(file = "Human_R47H_9_singlets_PCA.rds")
Human_R47H_10 <- readRDS(file = "Human_R47H_10_singlets_PCA.rds")
Human_R47H_11 <- readRDS(file = "Human_R47H_11_singlets_PCA.rds")
Human_R47H_13 <- readRDS(file = "Human_R47H_13_singlets_PCA.rds")
Human_R47H_14 <- readRDS(file = "Human_R47H_14_singlets_PCA.rds")
Human_R47H_15 <- readRDS(file = "Human_R47H_15_singlets_PCA.rds")
Human_R47H_16 <- readRDS(file = "Human_R47H_16_singlets_PCA.rds")
Human_R47H_18 <- readRDS(file = "Human_R47H_18_singlets_PCA.rds")
Human_R47H_19 <- readRDS(file = "Human_R47H_19_singlets_PCA.rds")
Human_R47H_20 <- readRDS(file = "Human_R47H_20_singlets_PCA.rds")
Human_R47H_21 <- readRDS(file = "Human_R47H_21_singlets_PCA.rds")
Human_R47H_22 <- readRDS(file = "Human_R47H_22_singlets_PCA.rds")
Human_R47H_23 <- readRDS(file = "Human_R47H_23_singlets_PCA.rds")
Human_R47H_24 <- readRDS(file = "Human_R47H_24_singlets_PCA.rds")
Human_R47H_25 <- readRDS(file = "Human_R47H_25_singlets_PCA.rds")
Human_R47H_26 <- readRDS(file = "Human_R47H_26_singlets_PCA.rds")
Human_R47H_27 <- readRDS(file = "Human_R47H_27_singlets_PCA.rds")
Human_R47H_29 <- readRDS(file = "Human_R47H_29_singlets_PCA.rds")
Human_R47H_30 <- readRDS(file = "Human_R47H_30_singlets_PCA.rds")
Human_R47H_31 <- readRDS(file = "Human_R47H_31_singlets_PCA.rds")
Human_R47H_33 <- readRDS(file = "Human_R47H_33_singlets_PCA.rds")
Human_R47H_35 <- readRDS(file = "Human_R47H_35_singlets_PCA.rds")
Human_R47H_37 <- readRDS(file = "Human_R47H_37_singlets_PCA.rds")
Human_R47H_38 <- readRDS(file = "Human_R47H_38_singlets_PCA.rds")
Human_R47H_39 <- readRDS(file = "Human_R47H_39_singlets_PCA.rds")
Human_R47H_40 <- readRDS(file = "Human_R47H_40_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Human_R47H/data_analysis/integration_new")

E2E4 <- c(Human_R47H_1, Human_R47H_2)
anchors_E2E4 <- FindIntegrationAnchors(object.list = E2E4, dims = 1:30)
E2E4_integrated <- IntegrateData(anchorset = anchors_E2E4, dims = 1:30)
rm(Human_R47H_1, Human_R47H_2)

R47H_E3E3_F <- c(Human_R47H_3, Human_R47H_7)
anchors_R47H_E3E3_F <- FindIntegrationAnchors(object.list = R47H_E3E3_F, dims = 1:30)
R47H_E3E3_F_integrated <- IntegrateData(anchorset = anchors_R47H_E3E3_F, dims = 1:30)
rm(Human_R47H_3, Human_R47H_7, R47H_E3E3_F)

R47H_E3E3_M <- c(Human_R47H_5, Human_R47H_9)
anchors_R47H_E3E3_M <- FindIntegrationAnchors(object.list = R47H_E3E3_M, dims = 1:30)
R47H_E3E3_M_integrated <- IntegrateData(anchorset = anchors_R47H_E3E3_M, dims = 1:30)
rm(Human_R47H_5, Human_R47H_9, R47H_E3E3_M)

WT_E3E3_F <- c(Human_R47H_4, Human_R47H_8)
anchors_WT_E3E3_F <- FindIntegrationAnchors(object.list = WT_E3E3_F, dims = 1:30)
WT_E3E3_F_integrated <- IntegrateData(anchorset = anchors_WT_E3E3_F, dims = 1:30)
rm(Human_R47H_4, Human_R47H_8, WT_E3E3_F)

WT_E3E3_M <- c(Human_R47H_6, Human_R47H_10)
anchors_WT_E3E3_M <- FindIntegrationAnchors(object.list = WT_E3E3_M, dims = 1:30)
WT_E3E3_M_integrated <- IntegrateData(anchorset = anchors_WT_E3E3_M, dims = 1:30)
rm(Human_R47H_6, Human_R47H_10, WT_E3E3_M)

R47H_E3E3 <- c(R47H_E3E3_F_integrated, R47H_E3E3_M_integrated)
anchors_R47H_E3E3 <- FindIntegrationAnchors(object.list = R47H_E3E3, dims = 1:30)
R47H_E3E3_integrated <- IntegrateData(anchorset = anchors_R47H_E3E3, dims = 1:30)
rm(R47H_E3E3_F_integrated, R47H_E3E3_M_integrated, R47H_E3E3)

WT_E3E3 <- c(WT_E3E3_F_integrated, WT_E3E3_M_integrated)
anchors_WT_E3E3 <- FindIntegrationAnchors(object.list = WT_E3E3, dims = 1:30)
WT_E3E3_integrated <- IntegrateData(anchorset = anchors_WT_E3E3, dims = 1:30)
rm(WT_E3E3_F_integrated, WT_E3E3_M_integrated, WT_E3E3)

E3E3 <- c(R47H_E3E3_integrated, WT_E3E3_integrated)
anchors_E3E3 <- FindIntegrationAnchors(object.list = E3E3, dims = 1:30)
E3E3_integrated <- IntegrateData(anchorset = anchors_E3E3, dims = 1:30)
rm(R47H_E3E3_integrated, WT_E3E3_integrated, E3E3)

R47H_E3E4_F <- c(Human_R47H_11, Human_R47H_13, Human_R47H_15, Human_R47H_21,Human_R47H_23,Human_R47H_27,Human_R47H_31)
anchors_R47H_E3E4_F <- FindIntegrationAnchors(object.list = R47H_E3E4_F, dims = 1:30)
R47H_E3E4_F_integrated <- IntegrateData(anchorset = anchors_R47H_E3E4_F, dims = 1:30)
rm(Human_R47H_11, Human_R47H_13, Human_R47H_15, Human_R47H_21,Human_R47H_23,Human_R47H_27,Human_R47H_31, R47H_E3E4_F)

R47H_E3E4_M <- c(Human_R47H_19, Human_R47H_25, Human_R47H_29)
anchors_R47H_E3E4_M <- FindIntegrationAnchors(object.list = R47H_E3E4_M, dims = 1:30)
R47H_E3E4_M_integrated <- IntegrateData(anchorset = anchors_R47H_E3E4_M, dims = 1:30)
rm(Human_R47H_19, Human_R47H_25, Human_R47H_29, R47H_E3E4_M)

WT_E3E4_F <- c(Human_R47H_14, Human_R47H_16, Human_R47H_18, Human_R47H_22,Human_R47H_24)
anchors_WT_E3E4_F <- FindIntegrationAnchors(object.list = WT_E3E4_F, dims = 1:30)
WT_E3E4_F_integrated <- IntegrateData(anchorset = anchors_WT_E3E4_F, dims = 1:30)
rm(Human_R47H_14, Human_R47H_16, Human_R47H_18, Human_R47H_22,Human_R47H_24, WT_E3E4_F)

WT_E3E4_M <- c(Human_R47H_20, Human_R47H_26, Human_R47H_30)
anchors_WT_E3E4_M <- FindIntegrationAnchors(object.list = WT_E3E4_M, dims = 1:30)
WT_E3E4_M_integrated <- IntegrateData(anchorset = anchors_WT_E3E4_M, dims = 1:30)
rm(Human_R47H_20, Human_R47H_26, Human_R47H_30, WT_E3E4_M)

R47H_E3E4 <- c(R47H_E3E4_F_integrated, R47H_E3E4_M_integrated)
anchors_R47H_E3E4 <- FindIntegrationAnchors(object.list = R47H_E3E4, dims = 1:30)
R47H_E3E4_integrated <- IntegrateData(anchorset = anchors_R47H_E3E4, dims = 1:30)
rm(R47H_E3E4_F_integrated, R47H_E3E4_M_integrated, R47H_E3E4)

WT_E3E4 <- c(WT_E3E4_F_integrated, WT_E3E4_M_integrated)
anchors_WT_E3E4 <- FindIntegrationAnchors(object.list = WT_E3E4, dims = 1:30)
WT_E3E4_integrated <- IntegrateData(anchorset = anchors_WT_E3E4, dims = 1:30)
rm(WT_E3E4_F_integrated, WT_E3E4_M_integrated, WT_E3E4)

E3E4 <- c(R47H_E3E4_integrated, WT_E3E4_integrated)
anchors_E3E4 <- FindIntegrationAnchors(object.list = E3E4, dims = 1:30)
E3E4_integrated <- IntegrateData(anchorset = anchors_E3E4, dims = 1:30)
rm(R47H_E3E4_integrated, WT_E3E4_integrated, E3E4)

R47H_E4E4_F <- c(Human_R47H_33, Human_R47H_39)
anchors_R47H_E4E4_F <- FindIntegrationAnchors(object.list = R47H_E4E4_F, dims = 1:30)
R47H_E4E4_F_integrated <- IntegrateData(anchorset = anchors_R47H_E4E4_F, dims = 1:30)
rm(Human_R47H_33, Human_R47H_39, R47H_E4E4_F)

R47H_E4E4_M <- c(Human_R47H_35, Human_R47H_37)
anchors_R47H_E4E4_M <- FindIntegrationAnchors(object.list = R47H_E4E4_M, dims = 1:30)
R47H_E4E4_M_integrated <- IntegrateData(anchorset = anchors_R47H_E4E4_M, dims = 1:30)
rm(Human_R47H_35, Human_R47H_37, R47H_E4E4_M)

R47H_E4E4 <- c(R47H_E4E4_F_integrated, R47H_E4E4_M_integrated)
anchors_R47H_E4E4 <- FindIntegrationAnchors(object.list = R47H_E4E4, dims = 1:30)
R47H_E4E4_integrated <- IntegrateData(anchorset = anchors_R47H_E4E4, dims = 1:30)
rm(R47H_E4E4_F_integrated, R47H_E4E4_M_integrated, R47H_E4E4)

WT_E4E4 <- c(Human_R47H_40, Human_R47H_38)
anchors_WT_E4E4 <- FindIntegrationAnchors(object.list = WT_E4E4, dims = 1:30)
WT_E4E4_integrated <- IntegrateData(anchorset = anchors_WT_E4E4, dims = 1:30)
rm(Human_R47H_40, Human_R47H_38, WT_E4E4)

E4E4 <- c(R47H_E4E4_integrated, WT_E4E4_integrated)
anchors_E4E4 <- FindIntegrationAnchors(object.list = E4E4, dims = 1:30)
E4E4_integrated <- IntegrateData(anchorset = anchors_E4E4, dims = 1:30)
rm(R47H_E4E4_integrated, WT_E4E4_integrated, E4E4)

R47H_all <- c(E2E4_integrated, E3E3_integrated, E3E4_integrated, E4E4_integrated)
anchors_R47H_all <- FindIntegrationAnchors(object.list = R47H_all, dims = 1:30)
R47H_all_integrated <- IntegrateData(anchorset = anchors_R47H_all, dims = 1:30)
rm(E2E4_integrated, E3E3_integrated, E3E4_integrated, E4E4_integrated, R47H_all)

pdf("R47H_all_integrated_QC_1.pdf", width=16, height=12)
Idents(R47H_all_integrated) <- "orig.ident"
VlnPlot(object = R47H_all_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size=0, idents=NULL)
dev.off()
pdf("R47H_all_integrated_QC_2.pdf", width=16, height=4)
Idents(R47H_all_integrated) <- "Condition"
VlnPlot(object = R47H_all_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(R47H_all_integrated, file = "R47H_all_integrated.rds")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

R47H_all_integrated <- readRDS("R47H_all_integrated.rds")
DefaultAssay(R47H_all_integrated) <- 'integrated'

R47H_all_integrated <- ScaleData(R47H_all_integrated, verbose = FALSE)
R47H_all_integrated <- RunPCA(R47H_all_integrated, features = VariableFeatures(object = R47H_all_integrated), verbose = FALSE)

R47H_all_integrated <- FindNeighbors(R47H_all_integrated, dims = 1:20)
R47H_all_integrated <- FindClusters(R47H_all_integrated, resolution = 0.1)
R47H_all_integrated <- RunUMAP(R47H_all_integrated, dims = 1: 20)

str(R47H_all_integrated)

DefaultAssay(R47H_all_integrated) <- 'RNA'
R47H_all_integrated <- NormalizeData(R47H_all_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
R47H_all_integrated <- ScaleData(R47H_all_integrated, features = rownames(R47H_all_integrated))

pdf("R47H_all_integrated_umap.pdf", width=5, height=4)
DimPlot(R47H_all_integrated, reduction = 'umap', label = T)
dev.off()
pdf("R47H_all_integrated_umap_split_individual.pdf", width=14, height=12)
DimPlot(R47H_all_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 6)
dev.off()
pdf("R47H_all_integrated_umap_split_Condition.pdf", width=6, height=12)
DimPlot(R47H_all_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

saveRDS(R47H_all_integrated, file = 'R47H_all_integrated_PCA_0.1.rds')


R47H_all_integrated <- readRDS("R47H_all_integrated_PCA_0.1.rds")

DefaultAssay(R47H_all_integrated) <- 'RNA'
pdf("R47H_all_integrated_umap_test.pdf", width=8, height=6)
DimPlot(R47H_all_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

pdf("R47H_all_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("MAP2","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2", "PLP1", "MBP", "MOBP","SCRG1", "OLIG1", "VTN", "MGP", "IGFBP7","CX3CR1", "P2RY12", "CSF1R","CLU", "PLA2G7", "CD3G", "CCL5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = R47H_all_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

R47H_all_markers <- FindAllMarkers(R47H_all_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(R47H_all_markers, "R47H_all_markers.csv")

R47H_all_integrated <- RenameIdents(R47H_all_integrated,
                                    `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="microglia",
                                    `4`="OPCs", `5`="inhibitory neurons", `6`="excitatory neurons", `7`="inhibitory neurons",
                                    `8`="inhibitory neurons", `9`="excitatory neurons", `10`="excitatory neurons", `11`="excitatory neurons",
                                    `12`="excitatory neurons", `13`="endothelial cells", `14`="inhibitory neurons",
                                    `15`="excitatory neurons", `16`="excitatory neurons", `17`="endothelial cells"
)

R47H_all_integrated$celltype.orig.ident <- paste(Idents(R47H_all_integrated), R47H_all_integrated$orig.ident, sep = "_")
R47H_all_integrated$celltype <- Idents(R47H_all_integrated)

Idents(R47H_all_integrated) <- "celltype"

R47H_all_MG <- subset(R47H_all_integrated, idents = "microglia")

saveRDS(R47H_all_integrated, file = "R47H_all_integrated_Annotation.rds")
saveRDS(R47H_all_MG, file = "R47H_all_MG.rds")




