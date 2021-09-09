###############################################################################################
# Pre-processing for data used to generate Figure 1 from Sayed, Kodama, Fan, et al. 2021 
# Total of 46 human AD R47H vs CV samples
# This script is: STEP 3 of 5 - Annotating cell types based on markers

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)


setwd("/athena/ganlab/scratch/lif4001/Human_AD_Mayo_UPenn/data_analysis/integration")
R47H_all_integrated <- readRDS("R47H_all_integrated_PCA_0.1.rds")
DefaultAssay(R47H_all_integrated) <- 'RNA'

pdf("R47H_all_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("SYT1","SNAP25","GRIN1","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2","PLP1", "MBP", "MOBP","AQP4","GFAP", 
           "CD74","CSF1R","C3","PDGFRA","VCAN","EBF1","IGFBP7","FLT1","CLDN5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = R47H_all_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

R47H_all_integrated <- subset(R47H_all_integrated, idents = "18", invert = T)

R47H_all_integrated <- RenameIdents(R47H_all_integrated,
                                 `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="microglia",
                                 `4`="OPCs", `5`="inhibitory neurons", `6`="excitatory neurons", `7`="excitatory neurons",
                                 `8`="inhibitory neurons", `9`="inhibitory neurons", `10`="inhibitory neurons", `11`="excitatory neurons",
                                 `12`="endothelial cells", `13`="oligodendrocytes", `14`="inhibitory neurons",
                                 `15`="excitatory neurons", `16`="excitatory neurons", `17`="excitatory neurons"
)



R47H_all_integrated$celltype.orig.ident <- paste(Idents(R47H_all_integrated), R47H_all_integrated$orig.ident, sep = "_")
R47H_all_integrated$celltype <- Idents(R47H_all_integrated)

Idents(R47H_all_integrated) <- "celltype"
pdf("R47H_all_integrated_umap_annotation.pdf", width=8, height=6)
DimPlot(R47H_all_integrated, reduction = 'umap', label = TRUE)
dev.off()

saveRDS(R47H_all_integrated, file = "R47H_all_integrated_Annotation.rds")