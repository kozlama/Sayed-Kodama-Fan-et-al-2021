###############################################################################################
# Pre-processing for data used to generate Figure 6 from Sayed, Kodama, Fan, et al. 2021 
# This script is STEP 2 of 3 - integration of all samples

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/LG72/DF_2ndRound")
CV_NT_Ctrl_1 <- readRDS(file = "CV_NT_Ctrl_1_singlets_PCA.rds")
CV_NT_Ctrl_2 <- readRDS(file = "CV_NT_Ctrl_2_singlets_PCA.rds")
CV_NT_Ctrl_3 <- readRDS(file = "CV_NT_Ctrl_3_singlets_PCA.rds")
CV_NT_Ctrl_4 <- readRDS(file = "CV_NT_Ctrl_4_singlets_PCA.rds")

R47H_NT_Ctrl_1 <- readRDS(file = "R47H_NT_Ctrl_1_singlets_PCA.rds")
R47H_NT_Ctrl_2 <- readRDS(file = "R47H_NT_Ctrl_2_singlets_PCA.rds")
R47H_NT_Ctrl_3 <- readRDS(file = "R47H_NT_Ctrl_3_singlets_PCA.rds")

R47H_NT_MK_1 <- readRDS(file = "R47H_NT_MK_1_singlets_PCA.rds")
R47H_NT_MK_2 <- readRDS(file = "R47H_NT_MK_2_singlets_PCA.rds")
R47H_NT_MK_3 <- readRDS(file = "R47H_NT_MK_3_singlets_PCA.rds")
R47H_NT_MK_4 <- readRDS(file = "R47H_NT_MK_4_singlets_PCA.rds")

R47H_PS19_Ctrl_1 <- readRDS(file = "R47H_PS19_Ctrl_1_singlets_PCA.rds")
R47H_PS19_Ctrl_2 <- readRDS(file = "R47H_PS19_Ctrl_2_singlets_PCA.rds")
R47H_PS19_Ctrl_3 <- readRDS(file = "R47H_PS19_Ctrl_3_singlets_PCA.rds")
R47H_PS19_Ctrl_4 <- readRDS(file = "R47H_PS19_Ctrl_4_singlets_PCA.rds")

R47H_PS19_MK_1 <- readRDS(file = "R47H_PS19_MK_1_singlets_PCA.rds")
R47H_PS19_MK_2 <- readRDS(file = "R47H_PS19_MK_2_singlets_PCA.rds")
R47H_PS19_MK_3 <- readRDS(file = "R47H_PS19_MK_3_singlets_PCA.rds")
R47H_PS19_MK_4 <- readRDS(file = "R47H_PS19_MK_4_singlets_PCA.rds")


setwd("/athena/ganlab/scratch/lif4001/LG72/integration_5genotypes")
CV_NT_Ctrl <- c(CV_NT_Ctrl_1, CV_NT_Ctrl_2,CV_NT_Ctrl_3,CV_NT_Ctrl_4)
anchors_CV_NT_Ctrl <- FindIntegrationAnchors(object.list = CV_NT_Ctrl, dims = 1:30)
CV_NT_Ctrl_integrated <- IntegrateData(anchorset = anchors_CV_NT_Ctrl, dims = 1:30)
rm(CV_NT_Ctrl_1, CV_NT_Ctrl_2,CV_NT_Ctrl_3,CV_NT_Ctrl_4, CV_NT_Ctrl)

R47H_NT_Ctrl <- c(R47H_NT_Ctrl_1, R47H_NT_Ctrl_2,R47H_NT_Ctrl_3)
anchors_R47H_NT_Ctrl <- FindIntegrationAnchors(object.list = R47H_NT_Ctrl, dims = 1:30)
R47H_NT_Ctrl_integrated <- IntegrateData(anchorset = anchors_R47H_NT_Ctrl, dims = 1:30)
rm(R47H_NT_Ctrl_1, R47H_NT_Ctrl_2,R47H_NT_Ctrl_3, R47H_NT_Ctrl)

R47H_NT_MK <- c(R47H_NT_MK_1, R47H_NT_MK_2,R47H_NT_MK_3,R47H_NT_MK_4)
anchors_R47H_NT_MK <- FindIntegrationAnchors(object.list = R47H_NT_MK, dims = 1:30)
R47H_NT_MK_integrated <- IntegrateData(anchorset = anchors_R47H_NT_MK, dims = 1:30)
rm(R47H_NT_MK_1, R47H_NT_MK_2,R47H_NT_MK_3,R47H_NT_MK_4, R47H_NT_MK)

R47H_PS19_Ctrl <- c(R47H_PS19_Ctrl_1, R47H_PS19_Ctrl_2,R47H_PS19_Ctrl_3,R47H_PS19_Ctrl_4)
anchors_R47H_PS19_Ctrl <- FindIntegrationAnchors(object.list = R47H_PS19_Ctrl, dims = 1:30)
R47H_PS19_Ctrl_integrated <- IntegrateData(anchorset = anchors_R47H_PS19_Ctrl, dims = 1:30)
rm(R47H_PS19_Ctrl_1, R47H_PS19_Ctrl_2,R47H_PS19_Ctrl_3,R47H_PS19_Ctrl_4, R47H_PS19_Ctrl)

R47H_PS19_MK <- c(R47H_PS19_MK_1, R47H_PS19_MK_2,R47H_PS19_MK_3,R47H_PS19_MK_4)
anchors_R47H_PS19_MK <- FindIntegrationAnchors(object.list = R47H_PS19_MK, dims = 1:30)
R47H_PS19_MK_integrated <- IntegrateData(anchorset = anchors_R47H_PS19_MK, dims = 1:30)
rm(R47H_PS19_MK_1, R47H_PS19_MK_2,R47H_PS19_MK_3,R47H_PS19_MK_4, R47H_PS19_MK)

LG72 <- c(CV_NT_Ctrl_integrated, R47H_NT_Ctrl_integrated, R47H_NT_MK_integrated,R47H_PS19_Ctrl_integrated, R47H_PS19_MK_integrated)
anchors_LG72 <- FindIntegrationAnchors(object.list = LG72, dims = 1:30)
LG72_integrated <- IntegrateData(anchorset = anchors_LG72, dims = 1:30)
rm(CV_NT_Ctrl_integrated, R47H_NT_Ctrl_integrated, R47H_NT_MK_integrated,R47H_PS19_Ctrl_integrated, R47H_PS19_MK_integrated, LG72)


pdf("LG72_QC.pdf", width=12, height=4)
Idents(LG72_integrated) <- "orig.ident"
VlnPlot(object = LG72_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(LG72_integrated) <- "Condition"
VlnPlot(object = LG72_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

#saveRDS(LG72_integrated, file = "LG72_integrated.rds")

DefaultAssay(LG72_integrated) <- 'integrated'

# LG72_integrated <- NormalizeData(LG72_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG72_integrated <- FindVariableFeatures(LG72_integrated, selection.method = "vst", nfeatures = 3000)

LG72_integrated <- ScaleData(LG72_integrated, verbose = FALSE)
LG72_integrated <- RunPCA(LG72_integrated, features = VariableFeatures(object = LG72_integrated), verbose = FALSE)

LG72_integrated <- FindNeighbors(LG72_integrated, dims = 1:15)
LG72_integrated <- FindClusters(LG72_integrated, resolution = 0.1)
LG72_integrated <- RunUMAP(LG72_integrated, dims = 1: 15)

DefaultAssay(LG72_integrated) <- 'RNA'
LG72_integrated <- NormalizeData(LG72_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG72_integrated <- ScaleData(LG72_integrated, features = rownames(LG72_integrated))

pdf("LG72_integrated_umap.pdf", width=5, height=4)
DimPlot(LG72_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG72_integrated_umap_split_individual.pdf", width=16, height=9)
DimPlot(LG72_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("LG72_integrated_umap_split_Condition.pdf", width=15, height=3)
DimPlot(LG72_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 5)
dev.off()

saveRDS(LG72_integrated, file = 'LG72_integrated_PCA_0.1.rds')

DefaultAssay(LG72_integrated) <- 'RNA'
pdf("LG72_umap_test.pdf", width=4, height=3)
DimPlot(LG72_integrated, reduction = 'umap', label = T)
dev.off()

DefaultAssay(LG72_integrated) <- 'RNA'

LG72_markers <- FindAllMarkers(LG72_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(LG72_markers, "LG72_markers.csv")

LG72_markers <- read.csv(file = "LG72_markers.csv", header=T,row.names =1)
top5 <- LG72_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5$gene <- as.character(top5$gene)
pdf("LG72_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(LG72_integrated, features = top5$gene) + NoLegend()
dev.off()

DefaultAssay(LG72_integrated) <- 'RNA'
pdf("LG72_umap_test.pdf", width=8, height=6)
DimPlot(LG72_integrated, reduction = 'umap', label = T)
dev.off()
#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG72_annotation_combine.pdf", width=10, height=5)
DotPlot(object = LG72_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

