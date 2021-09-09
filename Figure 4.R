#################################################################################
# Generation of Figure 4 from Sayed, Kodama, Fan, et al. 2021 
# Mouse single-cell microglia SMART-seq analysis post-QC (metrics defined in Methods)

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Lay Kodama
#################################################################################

library(Seurat)
library(plotly)
library(dplyr)
library(ggplot2)

#load in post QC data
setwd()
load("DataObj_QC_F_All.Robj", envir = parent.frame(), verbose = FALSE)

#Input metadata / genotype information
setwd()
pdata <- read.csv("Metadata_F_corrected.csv", header=T)
pdata.2<-read.csv("metaData_Faten.csv",header=T)
pdata$FACS_Tmem119 <- pdata.2$FACS_Tmem119[match(pdata$Sample.name, pdata.2$Sample_Name)]
pdata$FACS_CD45 <- pdata.2$FACS_CD45[match(pdata$Sample.name, pdata.2$Sample_Name)]

Column_names_pbmc<-as.data.frame(rownames(DataObj_QC_F_All@meta.data))
colnames(Column_names_pbmc)<-"Sample.name"
Column_names_pbmc$genotype<-pdata$Genotype.1[match(Column_names_pbmc$Sample.name,pdata$Sample.name)]
DataObj_QC_F_All <- AddMetaData(object = DataObj_QC_F_All , metadata = Column_names_pbmc$genotype, col.name = "Genotype_new")

#normalize
data <- NormalizeData(DataObj_QC_F_All, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
data <- ScaleData(object = data)
#perform and visualize PCA
data <- RunPCA(object = data, features = rownames(x = data), verbose = FALSE)

#selecting dimensions to run dimensional reduction
ElbowPlot(data, ndims = 50)
data <- RunTSNE(object = DataObj_QC_F_All, dims = 1:8, do.fast=TRUE, perplexity=100, max.iter=10000)
DimPlot(object = data, reduction = 'tsne')

#clustering and plotting tsNE  (Figure 4A,B)
data <- FindNeighbors(data, dims=1:6)
data <- FindClusters(data, resolution = 0.2)
DimPlot(object = data, reduction = 'tsne')

DimPlot(object = data, reduction = 'tsne',split.by = "Genotype_new")

#finding marker genes for clusters
setwd()
markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers,"cluster_markers_2.csv")
exp<-AverageExpression(object = data,return.seurat = TRUE)
avg<-GetAssayData(data)
write.csv(avg,"Expression.csv")

#proportion of clusters (Figure 4C)
prop.table(table(Idents(data)))
table <- as.matrix(table(Idents(data), data$Genotype_new))
sum=colSums(table)
table<-rbind(table[1,]/sum, table[2,]/sum)
barplot(as.matrix(table))

#DE markers
DE_Analysis_cluster<- FindMarkers(data,0)
write.csv(DE_Analysis_cluster, "DE_mg_cluster0.csv")

DE_Analysis_cluster<- FindMarkers(data,1)
write.csv(DE_Analysis_cluster, "DE_mg_cluster1.csv")

#volcano plots (Figure 4D)
Cluster1 = read.csv("DE_mg_cluster1.csv", header = TRUE)
with(Cluster1, plot(avg_logFC, -log10(p_val_adj), pch=20, cex=2, main="cluster 2", xlab ="log2FC[cluster 2 vs all]", ylab ="-log[FDR]", col="gray"))
with(subset(Cluster1, p_val_adj < 0.05 & avg_logFC >= 0.25), points(avg_logFC, -log10(p_val_adj), pch=20, cex=2, col="red"))
with(subset(Cluster1, p_val_adj < 0.05 & avg_logFC <= -0.25), points(avg_logFC, -log10(p_val_adj), pch=20, cex=2, col="blue"))

genes = c("Apoe",
          "Ms4a7",
          "Lyz2",
          "Axl",
          "Cybb",
          "Apobec1",
          "Lilrb4",
          "Ctsc",
          'Mrc1',
          "Tmem176a",
          "Sparc",
          "Cst3",
          "Selplg",
          "Slc2a5",
          "P2ry12",
          "Sall1",
          "Fscn1",
          "Gpr56",
          "Ctsl",
          "Cd9",
          "Hexb","Clec12a")
label = Cluster1[Cluster1$X %in% genes,]
with(subset(label),points(avg_logFC, -log10(p_val_adj), pch=1, cex=1.5, col="black"))

#feature plots (Figure 4E)
RidgePlot(data,features = c("Ptprc", "Itgam"), ncol = 2)
VlnPlot(data,features = c("Ptprc", "Itgam","Tmem119","P2ry12"), ncol = 2,pt.size=0)
FeaturePlot(data,features=c("Hexb"), cols=c("yellow","dark blue"))
VlnPlot(data,features = c("Cst7"), ncol = 2)

#comparison with DAM vs MG cluster1 (Figure 4F) ======
dam<-read.csv("DAM_signature.csv",header=T) # DAM genes published prerviously
mg<-read.csv("DE_mg_cluster1.csv",header=T)
damstage<-read.csv("Damstages.csv",header=T)

dam<- subset(dam,X.log10.DAM..p.value..Mann.Whitney.>-log10(0.05))
mg<- subset(mg,p_val_adj<0.05)
damstage<-subset(damstage,X.log10.p.value..Mann.Whitney>-log10(0.05))

mg$damlogFC <- dam$Fold.change..DAM.to.homeostatic.microglia.[match(mg$X,dam$Gene.name)]
mg$stagelogFC <- damstage$FC[match(mg$X,damstage$UNIQUD)]

mg.up<-mg[mg$avg_logFC>0,]
dam.up<-dam[dam$Fold.change..DAM.to.homeostatic.microglia.>1,]
overlap.up<-mg.up[mg.up$X%in%dam.up$Gene.name,]

mg.down<-mg[mg$avg_logFC<0,]
dam.down<-dam[dam$Fold.change..DAM.to.homeostatic.microglia.<-1,]
overlap.down<-mg.down[mg.down$X%in%dam.down$Gene.name,]

#scatterplot
mg$color <- "NC"
mg$color[mg$avg_logFC >0 & mg$damlogFC >0]<-"red"
mg$color[mg$avg_logFC <0 & mg$damlogFC <0]<-"blue"
write.csv(mg,"cluster1 overlap with DAM.csv")

ggplot(data=mg, aes(x=avg_logFC, y=damlogFC)) + geom_point(shape=19, alpha=0.5, size=3)  +
  theme_classic() +
  geom_hline(yintercept = c(0), linetype = "dotted") +
  geom_vline(xintercept = c(0), linetype = "dotted") +
  ylab("logFC[DAM vs homeostatic]") + xlab("logFC[Cluster1 vs all]") 
cor(mg$avg_logFC, mg$damlogFC, use = "complete.obs") #0.7908248
cor.test(mg$avg_logFC, mg$damlogFC, use = "complete.obs", method = "pearson")
