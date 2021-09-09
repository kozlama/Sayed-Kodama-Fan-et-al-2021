###############################################################################################
# Pre-processing for data used to generate Figure 1 from Sayed, Kodama, Fan, et al. 2021 
# Total of 46 human AD R47H vs CV samples
# This script is: STEP 5 of 5 - Generation of Figure panels as below

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
R47H_all_integrated$TREM2.Sex <- paste(R47H_all_integrated$TREM2, R47H_all_integrated$Sex, sep = "_")

# Define an order of cluster identities
my_levels <- c("WT_F","R47H_F","WT_M","R47H_M")
my_levels_celltype <- c("astrocytes","excitatory neurons","inhibitory neurons","microglia","oligodendrocytes","OPCs","endothelia cells")
# Relevel object@ident
R47H_all_integrated$TREM2.Sex <- factor(x = R47H_all_integrated$TREM2.Sex, levels = my_levels)
R47H_all_integrated$celltype <- factor(x = R47H_all_integrated$celltype, levels = my_levels_celltype)


setwd("/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration/Figures")
DefaultAssay(R47H_all_integrated) <- 'RNA'
### Fig.1B
pdf("R47H_all_integrated_umap_test_1.pdf", width=6.5, height=4)
DimPlot(R47H_all_integrated, reduction = 'umap', label = F, cols = c("#00B6EB","#F8766D","#C49A00","#00C094","#A58AFF","#006838","#FB61D7"))
dev.off()
### Fig.S1 E, F, G
pdf("R47H_all_integrated_QC.pdf", width=12, height=12)
Idents(R47H_all_integrated) <- "orig.ident"
VlnPlot(object = R47H_all_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size=0, idents=NULL)
dev.off()
### Fig.S1 H, I
pdf("R47H_all_integrated_FeatureScatter.pdf", width=10, height=5)
FeatureScatter(object = R47H_all_integrated, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", pt.size=0.1)
FeatureScatter(object = R47H_all_integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size=0.1)
dev.off()

# Fig. 1C: calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(R47H_all_integrated$TREM2.Sex,R47H_all_integrated$celltype))
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

ggsave("genotype_celltype_distribution_1.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration/Figures",
       width=8,height=8,units="in")

### Fig.S1J
pdf("R47H_all_integrated_Annotation.pdf", width=10, height=6)
FeaturePlot(R47H_all_integrated, features = c("FLT1","CLDN5","EBF1","GAD1","GAD2","PDGFRA","VCAN","CD74","C3","CSF1R","SLC17A7","CAMK2A","NRGN", "AQP4", "GFAP",
                                              "PLP1","MBP","MOBP"))
dev.off()


# Fig. S1K: calculate ratio of each cell type in each sample
a<-as.data.frame(table(R47H_all_integrated$orig.ident,R47H_all_integrated$celltype))
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

ggsave("sample_celltype_distribution_1.pdf",plot=last_plot(),path="/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration/Figures",
       width=10,height=4,units="in")




