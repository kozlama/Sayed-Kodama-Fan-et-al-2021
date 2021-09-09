###############################################################################################
# Pre-processing for data used to generate Figure 2 from Sayed, Kodama, Fan, et al. 2021 
# Total of 46 human AD R47H vs CV samples
# This script is to integrate the microglia cells only from the UPenn samples

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====
#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Human_UPenn/data_analysis/DF_2ndRound")
Gan_43 <- readRDS(file = "Gan_43_singlets.rds")
Gan_44 <- readRDS(file = "Gan_44_singlets.rds")
Gan_45 <- readRDS(file = "Gan_45_singlets.rds")
Gan_46 <- readRDS(file = "Gan_46_singlets.rds")
Gan_47 <- readRDS(file = "Gan_47_singlets.rds")
Gan_48 <- readRDS(file = "Gan_48_singlets.rds")
Gan_49 <- readRDS(file = "Gan_49_singlets.rds")
Gan_50 <- readRDS(file = "Gan_50_singlets.rds")
Gan_51 <- readRDS(file = "Gan_51_singlets.rds")
Gan_52 <- readRDS(file = "Gan_52_singlets.rds")
Gan_54 <- readRDS(file = "Gan_54_singlets.rds")
Gan_55 <- readRDS(file = "Gan_55_singlets.rds")
Gan_56 <- readRDS(file = "Gan_56_singlets.rds")
Gan_59 <- readRDS(file = "Gan_59_singlets.rds")
Gan_60 <- readRDS(file = "Gan_60_singlets.rds")
Gan_61 <- readRDS(file = "Gan_61_singlets.rds")
Gan_62 <- readRDS(file = "Gan_62_singlets.rds")
Gan_63 <- readRDS(file = "Gan_63_singlets.rds")
Gan_64 <- readRDS(file = "Gan_64_singlets.rds")
Gan_65 <- readRDS(file = "Gan_65_singlets.rds")
Gan_66 <- readRDS(file = "Gan_66_singlets.rds")

setwd("/athena/ganlab/scratch/lif4001/Human_UPenn/data_analysis/integration_singlets")

E3E3_F <- c(Gan_43, Gan_44)
anchors_E3E3_F <- FindIntegrationAnchors(object.list = E3E3_F, dims = 1:30)
E3E3_F_integrated <- IntegrateData(anchorset = anchors_E3E3_F, dims = 1:30)
rm(Gan_43, Gan_44, E3E3_F)

E3E3 <- c(E3E3_F_integrated, Gan_46)
anchors_E3E3 <- FindIntegrationAnchors(object.list = E3E3, dims = 1:30)
E3E3_integrated <- IntegrateData(anchorset = anchors_E3E3, dims = 1:30)
rm(E3E3_F_integrated, Gan_46, E3E3)

E4E4 <- c(Gan_45, Gan_50)
anchors_E4E4 <- FindIntegrationAnchors(object.list = E4E4, dims = 1:30)
E4E4_integrated <- IntegrateData(anchorset = anchors_E4E4, dims = 1:30)
rm(Gan_45, Gan_50, E4E4)

E3E4 <- c(Gan_47, Gan_48, Gan_49)
anchors_E3E4 <- FindIntegrationAnchors(object.list = E3E4, dims = 1:30)
E3E4_integrated <- IntegrateData(anchorset = anchors_E3E4, dims = 1:30)
rm(Gan_47, Gan_48, Gan_49, E3E4)

AD_WT <- c(E3E3_integrated, E4E4_integrated, E3E4_integrated)
anchors_AD_WT <- FindIntegrationAnchors(object.list = AD_WT, dims = 1:30)
AD_WT_integrated <- IntegrateData(anchorset = anchors_AD_WT, dims = 1:30)
rm(E3E3_integrated, E4E4_integrated, E3E4_integrated, AD_WT)

R47H_E3E3 <- c(Gan_52, Gan_54)
anchors_R47H_E3E3 <- FindIntegrationAnchors(object.list = R47H_E3E3, dims = 1:30)
R47H_E3E3_integrated <- IntegrateData(anchorset = anchors_R47H_E3E3, dims = 1:30)
rm(Gan_52, Gan_54, R47H_E3E3)

R47H_E3E4 <- c(Gan_55, Gan_56)
anchors_R47H_E3E4 <- FindIntegrationAnchors(object.list = R47H_E3E4, dims = 1:30)
R47H_E3E4_integrated <- IntegrateData(anchorset = anchors_R47H_E3E4, dims = 1:30)
rm(Gan_55, Gan_56, R47H_E3E4)

AD_R47H <- c(Gan_51, R47H_E3E3_integrated, R47H_E3E4_integrated)
anchors_AD_R47H <- FindIntegrationAnchors(object.list = AD_R47H, dims = 1:30)
AD_R47H_integrated <- IntegrateData(anchorset = anchors_AD_R47H, dims = 1:30)
rm(Gan_51, R47H_E3E3_integrated, R47H_E3E4_integrated, AD_R47H)

AD <- c(AD_WT_integrated, AD_R47H_integrated)
anchors_AD <- FindIntegrationAnchors(object.list = AD, dims = 1:30)
AD_integrated <- IntegrateData(anchorset = anchors_AD, dims = 1:30)
rm(AD_WT_integrated, AD_R47H_integrated, AD)

Non_E2E3_F <- c(Gan_59, Gan_60)
anchors_Non_E2E3_F <- FindIntegrationAnchors(object.list = Non_E2E3_F, dims = 1:30)
Non_E2E3_F_integrated <- IntegrateData(anchorset = anchors_Non_E2E3_F, dims = 1:30)
rm(Gan_59, Gan_60, Non_E2E3_F)

Non_E2E3_M <- c(Gan_62, Gan_63)
anchors_Non_E2E3_M <- FindIntegrationAnchors(object.list = Non_E2E3_M, dims = 1:30)
Non_E2E3_M_integrated <- IntegrateData(anchorset = anchors_Non_E2E3_M, dims = 1:30)
rm(Gan_62, Gan_63, Non_E2E3_M)

Non_E2E3 <- c(Non_E2E3_F_integrated, Non_E2E3_M_integrated)
anchors_Non_E2E3 <- FindIntegrationAnchors(object.list = Non_E2E3, dims = 1:30)
Non_E2E3_integrated <- IntegrateData(anchorset = anchors_Non_E2E3, dims = 1:30)
rm(Non_E2E3_F_integrated, Non_E2E3_M_integrated, Non_E2E3)

Non_E3E3 <- c(Gan_65, Gan_66)
anchors_Non_E3E3 <- FindIntegrationAnchors(object.list = Non_E3E3, dims = 1:30)
Non_E3E3_integrated <- IntegrateData(anchorset = anchors_Non_E3E3, dims = 1:30)
rm(Gan_65, Gan_66, Non_E3E3)

Non <- c(Non_E2E3_integrated, Gan_61, Gan_64, Non_E3E3_integrated)
anchors_Non <- FindIntegrationAnchors(object.list = Non, dims = 1:30)
Non_integrated <- IntegrateData(anchorset = anchors_Non, dims = 1:30)
rm(Non_E2E3_integrated, Gan_61, Gan_64, Non_E3E3_integrated, Non)

UPenn <- c(AD_integrated, Non_integrated)
anchors_UPenn <- FindIntegrationAnchors(object.list = UPenn, dims = 1:30)
UPenn_integrated <- IntegrateData(anchorset = anchors_UPenn, dims = 1:30)
rm(AD_integrated, Non_integrated, UPenn)

pdf("UPenn_integrated_QC_1.pdf", width=16, height=12)
Idents(UPenn_integrated) <- "orig.ident"
VlnPlot(object = UPenn_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size=0, idents=NULL)
dev.off()
pdf("UPenn_integrated_QC_2.pdf", width=16, height=4)
Idents(UPenn_integrated) <- "Condition"
VlnPlot(object = UPenn_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(UPenn_integrated, file = "UPenn_integrated.rds")

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

#UPenn_integrated <- readRDS("UPenn_integrated.rds")
DefaultAssay(UPenn_integrated) <- 'integrated'

#UPenn_integrated <- NormalizeData(UPenn_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
#UPenn_integrated <- FindVariableFeatures(UPenn_integrated, selection.method = "vst", nfeatures = 3000)

UPenn_integrated <- ScaleData(UPenn_integrated, verbose = FALSE)
UPenn_integrated <- RunPCA(UPenn_integrated, features = VariableFeatures(object = UPenn_integrated), verbose = FALSE)

UPenn_integrated <- FindNeighbors(UPenn_integrated, dims = 1:20)
UPenn_integrated <- FindClusters(UPenn_integrated, resolution = 0.1)
UPenn_integrated <- RunUMAP(UPenn_integrated, dims = 1: 20)

str(UPenn_integrated)

DefaultAssay(UPenn_integrated) <- 'RNA'
UPenn_integrated <- NormalizeData(UPenn_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
UPenn_integrated <- ScaleData(UPenn_integrated, features = rownames(UPenn_integrated))

pdf("UPenn_integrated_umap.pdf", width=6, height=4)
DimPlot(UPenn_integrated, reduction = 'umap', label = T)
dev.off()
pdf("UPenn_integrated_umap_split_individual.pdf", width=14, height=12)
DimPlot(UPenn_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 6)
dev.off()
pdf("UPenn_integrated_umap_split_Condition.pdf", width=6, height=12)
DimPlot(UPenn_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()

saveRDS(UPenn_integrated, file = 'UPenn_integrated_PCA_0.1.rds')

DefaultAssay(UPenn_integrated) <- 'RNA'
pdf("UPenn_integrated_umap_test.pdf", width=8, height=6)
DimPlot(UPenn_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

pdf("UPenn_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("MAP2","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2", "PLP1", "MBP", "MOBP","SCRG1", "OLIG1", "VTN", "MGP", "IGFBP7","CX3CR1", "P2RY12", "CSF1R","CLU", "PLA2G7", "CD3G", "CCL5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = UPenn_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

UPenn_integrated <- RenameIdents(UPenn_integrated,
                                 `0` = "oligodendrocytes", `1`="astrocytes", `2`="microglia", `3`="excitatory neurons",
                                 `4`="OPCs", `5`="inhibitory neurons", `6`="inhibitory neurons", `7`="excitatory neurons",
                                 `8`="excitatory neurons", `9`="inhibitory neurons", `10`="oligodendrocytes", `11`="excitatory neurons",
                                 `12`="endothelial cells", `13`="endothelial cells", `14`="inhibitory neurons",
                                 `15`="excitatory neurons"
)

UPenn_integrated$celltype.orig.ident <- paste(Idents(UPenn_integrated), UPenn_integrated$orig.ident, sep = "_")
UPenn_integrated$celltype <- Idents(UPenn_integrated)

Idents(UPenn_integrated) <- "celltype"

UPenn_MG <- subset(UPenn_integrated, idents = "microglia")

saveRDS(UPenn_integrated, file = "UPenn_integrated_Annotation.rds")
saveRDS(UPenn_MG, file = "UPenn_MG.rds")


