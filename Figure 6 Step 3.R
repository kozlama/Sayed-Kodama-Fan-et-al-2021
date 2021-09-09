###############################################################################################
# Pre-processing for data used to generate Figure 6 from Sayed, Kodama, Fan, et al. 2021 
# This script is STEP 3 of 3:
# annotation of cell types based on markers
# subsetting microglia 
# sub-clustering of microglia 
# generation of figure panels

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/LG72/integration_5genotypes")
LG72_integrated <- readRDS("LG72_integrated_PCA_0.1.rds")
Idents(LG72_integrated) <- "seurat_clusters"
LG72_integrated <- RenameIdents(LG72_integrated,
                                 `0` = "oligodendrocytes", `1`="excitatory neurons", `2`="astrocytes", `3`="excitatory neurons",
                                 `4`="excitatory neurons", `5`="microglia", `6`="excitatory neurons", `7`="excitatory neurons",
                                 `8`="excitatory neurons", `9`="inhibitory neurons", `10`="OPCs", `11`="excitatory neurons",
                                 `12`="inhibitory neurons", `13`="endothelial cells", `14`="inhibitory neurons", `15`="excitatory neurons",
                                 `16`="excitatory neurons", `17`="astro-ependymal cells"
)

setwd("/athena/ganlab/scratch/lif4001/LG72/integration_5genotypes/annotation")
### Fig.10A
pdf("LG72_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(LG72_integrated, reduction = 'umap', label = F)
dev.off()

LG72_integrated$celltype.orig.ident <- paste(Idents(LG72_integrated), LG72_integrated$orig.ident, sep = "_")
LG72_integrated$celltype <- Idents(LG72_integrated)

Idents(LG72_integrated) <- "Condition"
LG72_integrated <- RenameIdents(LG72_integrated,
                                `CV_NT_Ctrl` = "mTrem2+/+ Veh", `R47H_NT_Ctrl` = "R47H/+ Veh", `R47H_NT_MK` = "R47H/+ MK",
                                `R47H_PS19_Ctrl` = "P301S R47H/+ Veh", `R47H_PS19_MK` = "P301S R47H/+ MK")
LG72_integrated$Condition <- Idents(LG72_integrated)

Idents(LG72_integrated) <- "Sample_Name"
LG72_integrated <- RenameIdents(LG72_integrated,
                                `CV_NT_Ctrl_1` = "mTrem2+/+ Veh_1", `CV_NT_Ctrl_2` = "mTrem2+/+ Veh_2", `CV_NT_Ctrl_3` = "mTrem2+/+ Veh_3", `CV_NT_Ctrl_4` = "mTrem2+/+ Veh_4", 
                                `R47H_NT_Ctrl_1` = "R47H/+ Veh_1", `R47H_NT_Ctrl_2` = "R47H/+ Veh_2", `R47H_NT_Ctrl_3` = "R47H/+ Veh_3", 
                                `R47H_NT_MK_1` = "R47H/+ MK_1",`R47H_NT_MK_2` = "R47H/+ MK_2",`R47H_NT_MK_3` = "R47H/+ MK_3",`R47H_NT_MK_4` = "R47H/+ MK_4",
                                `R47H_PS19_Ctrl_1` = "P301S R47H/+ Veh_1",`R47H_PS19_Ctrl_2` = "P301S R47H/+ Veh_2",`R47H_PS19_Ctrl_3` = "P301S R47H/+ Veh_3",`R47H_PS19_Ctrl_4` = "P301S R47H/+ Veh_4",
                                `R47H_PS19_MK_1` = "P301S R47H/+ MK_1",`R47H_PS19_MK_2` = "P301S R47H/+ MK_2",`R47H_PS19_MK_3` = "P301S R47H/+ MK_3",`R47H_PS19_MK_4` = "P301S R47H/+ MK_4")
LG72_integrated$Sample_Name <- Idents(LG72_integrated)


data <- LG72_integrated
# Fig.10C: calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Condition,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/LG72/integration_5genotypes/annotation",
       width=4,height=4,units="in")

data <- LG72_integrated
# Fig.10D: calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$Sample_Name,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/LG72/integration_5genotypes/annotation",
       width=8,height=4,units="in")

### Fig.10 E, F, G
pdf("LG72_QC_sample.pdf", width=16, height=4)
Idents(LG72_integrated) <- "Sample_Name"
VlnPlot(object = LG72_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

### Fig.10 H, I
pdf("LG72_FeatureScatter_gene_UMI.pdf", width=6, height=4)
FeatureScatter(object = LG72_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample_Name", pt.size=0.5)
dev.off()
pdf("LG72_FeatureScatter_mt_UMI_1.pdf", width=6, height=6)
FeatureScatter(object = LG72_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample_Name", pt.size=0.5)
dev.off()

data <- LG72_integrated
#Fig.10B: markers for annotation
pdf("annotation_1.pdf", width=9, height=3)
DotPlot(data, features = c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Gad1", "Gad2","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1","Clic6","Ttr")) + RotatedAxis()
dev.off()

MG <- subset(LG72_integrated, idents="5")
saveRDS(MG, file = "LG72_MG_subset.rds")

setwd("~/Desktop/data_analysis/LG72/integration_5genotypes/MG_only/final")

DefaultAssay(MG) <- 'integrated'

MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:12)
MG <- FindClusters(MG, resolution = 0.1)
MG <- RunUMAP(MG, dims = 1: 12)

DimPlot(MG, reduction = 'umap', label = T)

MG <- subset(MG, idents = c("0","1","2","3"))

n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)

### Fig. 6G
DimPlot(MG, reduction = 'umap', label = T)
### Fig. 6H
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 5)

#7/26 Table S10
DefaultAssay(MG) <- 'RNA'
LG72_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0, test.use = "MAST",min.pct = 0.25)
write.csv(LG72_MG_markers, "LG72_MG_markers_5genotypes_final.csv")

LG72_MG4_markers <- FindMarkers(MG, logfc.threshold = 0, test.use = "MAST", ident.1 = "4",min.pct = 0.25)
write.csv(LG72_MG4_markers, "LG72_MG4_final_DEG_for_IPA.csv")

















