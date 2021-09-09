#################################################################################
# Generation of Figure 2 and Figure S2 and S3 from Sayed, Kodama, Fan, et al. 2021 
# Human microglia analysis post-processing of 46 human sample
# For pre-processing steps from raw data, refer to Figure 1 pre-processing codes

# by Lay Kodama
#################################################################################

library(scales)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(reshape2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(Seurat)

# Import pre-processed and integrated RDS file of microglia from 46 human samples (Figure 2A)
setwd()
mg <- readRDS("seurat_integrated_microglia_AD.rds")
colnames(mg@meta.data)[8] <- "TREM2_genotype"

DimPlot(mg,group.by = "orig.ident")

mg <- FindNeighbors(mg, dims = 1:30)
mg <- FindClusters(mg, resolution = 0.4)
DimPlot(mg, label=T, label.size = 5) + NoLegend()

# Find markers of microglia subclusters ======
DefaultAssay(mg) <- "RNA"
markers <- FindAllMarkers(mg)
write.csv(markers, "microglia_markers_with_names.csv")

# rename seurat clusters based on markers from above:
cluster_names <- data.frame(seurat_cluster = 0:11, 
                            cluster_names = c("MG1","MG2","MG3","MG4","MG5","N1","MG6","OG1","MAC1","MG7","MAC2","N2"))

# Dot plot of markers (Figure 2B):
mg@meta.data$cluster_names <- cluster_names$cluster_names[match(mg@meta.data$seurat_clusters, cluster_names$seurat_cluster)]

Idents(mg) <- "cluster_names"
DefaultAssay(mg) <- "RNA"
markers.to.plot <- c("SIGLEC1","MRC1","PTPRC","CD14","P2RY12","CX3CR1","TREM2","SYT1","NRXN1","NRG3","MOG","PLP1")
DotPlot(mg, features = markers.to.plot, col.max = 1, col.min = -1, dot.scale = 5, scale.by = "radius",scale.max = 40) +
  scale_colour_gradient2(low = "darkblue", mid = "white", high = "red", breaks = c(-1,0,1))+
  RotatedAxis()+
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

identities <- levels(mg@active.ident)
my_color_palette <- hue_pal()(length(identities))
DimPlot(mg, label = T, label.size = 5) + NoLegend() +
  scale_color_manual(values = c(my_color_palette[8:9],my_color_palette[1:7],"grey","grey","grey"))

# plot proportion of cells per genotype (Figure 2C) =====
mg@meta.data$Trem2_sex <- as.character(paste(mg@meta.data$TREM2,mg@meta.data$Sex,sep="_"))

a <- as.data.frame(table(mg@meta.data$cluster_names, mg@meta.data$orig.ident))
a$Sex <- mg@meta.data$Sex[match(a$Var2,mg@meta.data$orig.ident)]
a$Trem2 <- mg@meta.data$TREM2[match(a$Var2,mg@meta.data$orig.ident)]
a$Trem2_sex <- mg@meta.data$Trem2_sex[match(a$Var2,mg@meta.data$orig.ident)]
a$gan_number <- mg@meta.data$gan_number[match(a$Var2,mg@meta.data$orig.ident)]

colnames(a)<-c("clusters","sample","cell.no","sex","Trem2","Trem2_sex","ID")

agg<-aggregate(cell.no~sample,a,sum)
a$cluster.total <- agg$cell.no[match(a$sample,agg$sample)]
a$ratio<-a$cell.no/a$cluster.total
a$Trem2 <- factor(a$Trem2, levels=c("WT","R47H"))
mg_only <- subset(a, clusters == "MG1"|
                    clusters  == "MG2"|
                    clusters =="MG3"|
                    clusters =="MG4"|
                    clusters =="MG5"|
                    clusters =="MG6"|
                    clusters =="MG7")
ggplot(mg_only,aes(x=Trem2, y=ratio, fill=clusters))+
  geom_boxplot()+
  geom_point(size = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cluster ratio per genotype")+
  facet_grid(~clusters)+
  NoLegend()+
  scale_fill_manual(values = my_color_palette[1:7])

# cell ratio by Trem2 genotype and sex (Figure S2B)
a$Trem2_sex <- factor(a$Trem2_sex, levels=c("WT_F","R47H_F","WT_M","R47H_M"))
mg_only <- subset(a, clusters == "MG1"|
                    clusters  == "MG2"|
                    clusters =="MG3"|
                    clusters =="MG4"|
                    clusters =="MG5"|
                    clusters =="MG6"|
                    clusters =="MG7")
ggplot(mg_only,aes(x=Trem2_sex, y=ratio, fill=clusters))+
  geom_boxplot()+
  geom_point(size = 1)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cluster ratio per genotype")+
  facet_grid(~clusters)+
  NoLegend()+
  scale_fill_manual(values = my_color_palette[1:7])

# Calculate statistics between R47H vs CV cluster ratios (Figure 2C)====
abundances <- table(mg@meta.data$seurat_clusters, mg@meta.data$orig.ident)
abundances <- unclass(abundances)
abundances <- abundances[c(1:5, 7,9,10,11),]

# Extract meta data from Seurat object
metadata <- mg@meta.data
extra.info <- metadata[match(colnames(abundances), mg@meta.data$orig.ident),]

# Create EdgeR object
y.ab <- DGEList(abundances, samples=extra.info)
y.ab

# Design matrix
design <- model.matrix(~factor(project) + factor(Sex) + Age + factor(APOE) + factor(TREM2_genotype), y.ab$samples)

# Estimate dispersion
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)

# QL dispersion
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

# Calculate GLM for Condition of sample
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
results <- topTags(res)

# GSEA analysis in comparison with microglia markers published in literature (Figure 2D)=====
# file export for GSEA
setwd()

markers <- read.csv("microglia_markers_with_names.csv")

markers$p_val_adj[markers$p_val_adj == 0] <- 1*10^-300
for (i in (0:11)){
  marker.sub <- subset(markers, cluster == i & p_val_adj < 0.05)
  marker.sub$sign <- sign(marker.sub$avg_logFC)
  marker.sub$rank <- marker.sub$sign * 1/marker.sub$p_val_adj
  rank <- data.frame(gene = marker.sub$gene, rank = marker.sub$rank)
  rank$gene <- as.character(rank$gene)
  write.table(rank, file = paste("cluster_",i,".txt",sep=""), sep = "\t", quote = F,
              row.names = FALSE, col.names = TRUE)
}

# plot of GSEA results
# Note: cluster 5 and 7 are not microglia cells based on markers
setwd()
gsea.result <- read.csv("GSEA_results_compiled.csv",header=T)
gsea.result$NOM.p.val[gsea.result$NOM.p.val == 0] <- 0.0000000000001

gsea.result <- gsea.result %>%
  dplyr::mutate(p_sign = sign(NES) * -log10(NOM.p.val))
gsea_cast <- dcast(gsea.result, Cluster_num ~ NAME, value.var = "p_sign")
row.names(gsea_cast) <- gsea_cast[,1]

p_star <- matrix(ncol= ncol(gsea_cast[,-1]), nrow = nrow(gsea_cast[-c(6,8),]))
p_star[gsea_cast[-c(6,8),-1] > -log10(0.05) | gsea_cast[-c(6,8),-1] < log10(0.05)] <- "*"
p_star[gsea_cast[-c(6,8),-1] > -log10(0.01) | gsea_cast[-c(6,8),-1] < log10(0.01)] <- "**"
p_star[gsea_cast[-c(6,8),-1] > -log10(0.001) | gsea_cast[-c(6,8),-1] < log10(0.001)] <- "***"
p_star[gsea_cast[-c(6,8),-1] > -log10(0.0001) | gsea_cast[-c(6,8),-1] < log10(0.0001)] <- "****"
p_star[is.na(p_star)] <- ""

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2))
Heatmap(as.matrix(gsea_cast[-c(6,8),-1]), cluster_rows = F, cluster_columns = F, col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", p_star[i, j]), x, y, gp = gpar(fontsize = 10))
        })

# volcano plot of markers of MG4, focusing on DAM genes (Figure 2E) ======
setwd()
dam <- read.csv("DAM.csv",header=T)
dam$gene_cap <- toupper(dam$gene)

marker_vp <- function(cluster_number){
  cluster4.de <- subset(markers, cluster == cluster_number & p_val_adj < 0.05)
  cluster4.de$color <- "black"
  cluster4.de$color[cluster4.de$avg_logFC > 0.25] <- "red"
  cluster4.de$color[cluster4.de$avg_logFC < -0.25] <- "blue"
  cluster4.de$dam <- "no"
  
  cluster4.de$dam[subset(cluster4.de, avg_logFC > 0)$gene %in% subset(dam,log.pvalue. > -log10(0.05) & foldchange > 0)$gene_cap] <- "yes"
  cluster4.de$dam[subset(cluster4.de, avg_logFC < 0)$gene %in% subset(dam,log.pvalue. > -log10(0.05) & foldchange < 0)$gene_cap] <- "yes"
  
  ggplot(cluster4.de, aes(x = avg_logFC, y = -log(p_val_adj)))+
    geom_point(aes(color= color))+
    scale_color_manual(values = c("black" = "black", "red" = "salmon", "blue" = "dodgerblue"))+
    theme_classic()+
    geom_vline(xintercept = c(0), linetype = "dotted")+
    geom_hline(yintercept = -log(0.05), linetype = "dotted")+
    geom_text_repel(subset(cluster4.de, dam == "yes" & abs(avg_logFC) > 0.4 ), mapping = aes(label = gene))+
    ggtitle(paste("Cluster ",cluster_number, " markers", sep=""))+
    ylab("-Log10(FDR)")+
    xlab("Normalized Log2(FC)")+
    NoLegend()
}

marker_vp(3) # MG4

# Pathway analysis of MG4 markers ==================
# Barplot of MG4 Hallmark GSEA overlap results (Figure 2F)
setwd()
gsea_mg4 <- read.csv("Hallmark_GSEA_MG4_overlap.csv",header=T)
ggplot(gsea_mg4, aes(x=reorder(Gene.Set.Name, -log10(FDR.q.value)), y=-log10(FDR.q.value)))+
  geom_bar(stat = "identity", fill = "chartreuse3")+
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "seashell4")+
  coord_flip()+
  theme_classic()+
  xlab("")+
  ylab("Predicted activation z score")+
  ggtitle("MG4 marker pathways")

# Save bar plot
ggsave("GSEA hallmark MG4 markers.pdf", plot = last_plot(), device = "pdf",
       scale = 1, width = 6, height = 3, units = c("in"))

# IPA upstream analysis for cluster 4 markers (Figure 2G)
setwd()
df <- read.csv("TREM2_IPA.csv",header=T)
df <- df[,c(1,5,6)]
colnames(df) <- c("upstream", "activation", "pvalue")
df <- subset(df, ! is.na(activation))
df <- as.data.frame(df)
ggplot(df, aes(x=reorder(upstream,activation), y=activation, fill=-log(pvalue)))+
  geom_bar(stat = "identity")+
  scale_fill_gradient2(low = "white",
                       high = "red")+
  coord_flip()+
  theme_classic()+
  xlab("")+
  ylab("Predicted activation z score")+
  ggtitle("Trem2 signaling")

# Save bar plot
ggsave("Trem2 signaling MG4 IPA.pdf", plot = last_plot(), device = "pdf",
       scale = 1, width = 4.5, height = 3, units = c("in"))

# heatmap of overlapping genes between MG markers and DAM (Figure S3) =====
dam <- read.csv("DAM.csv",header=T)
dam$gene_cap <- toupper(dam$gene)
dam_sig <- subset(dam, log.pvalue. > -log10(0.05))

mg4_markers <- subset(markers, p_val_adj < 0.05 & cluster == 3)
mg_markers <- subset(markers, p_val_adj < 0.05 &
                       gene %in% dam_sig$gene_cap &
                       gene %in% mg4_markers$gene)
mg_markers <- subset(mg_markers, ! cluster == 5 &
                       ! cluster == 7 &
                       ! cluster == 8 &
                       ! cluster == 10 &
                       ! cluster == 11)
mg_cast <- dcast(mg_markers, gene~cluster, value.var = "avg_logFC")
mg_cast$dam <- dam$foldchange[match(mg_cast$gene, dam$gene_cap)]

row.names(mg_cast) <- mg_cast[,1]
mg_cast <- mg_cast[,-1]

mg_cast[is.na(mg_cast)] <- 0

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-3, 3))
Heatmap(mg_cast, cluster_columns = T, cluster_rows = T, col = col_fun)
