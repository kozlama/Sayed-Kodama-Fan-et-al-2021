###############################################################################################
# Code used to generate Figure 6J from Sayed, Kodama, Fan, et al. 2021 

# by Li Fan 
###############################################################################################

library(ggplot2)
library(Seurat)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)

# Find markers 
marker <- read.csv(file = "LG72_MG4_final_DEG_for_IPA.csv", header=T,row.names =1)

# Volcano plot of marker genes =====

# Identify DEGs for cluster 5 markers
marker$colours <- c("NC")
marker$colours[marker$avg_log2FC >= 0.25 & marker$p_val_adj <= 0.05] <- c("UP")
marker$colours[marker$avg_log2FC <= -0.25 & marker$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_mature <- c("Bank1", "Stab1", "Siglech", "Selplg","P2ry12","Hexb","Cx3cr1","Elmo1","Csf1r")
genes_to_plot_mature <- marker[row.names(marker) %in% genes_select_mature, ]
genes_to_plot_mature$Cluster <- "Mature"

genes_select_immature <- c("H2-K1", "Apbb2","Apoe","Myo1e","Cd9","Stat1","Arhgap24","Ctsz","Parp14","Ddx60","Arid5b")
genes_to_plot_immature <- marker[row.names(marker)  %in% genes_select_immature, ]
genes_to_plot_immature$Cluster <- c("Immature")

genes_to_plot <- rbind(genes_to_plot_mature, genes_to_plot_immature)

# Set color palette
my_color <- c("#2B8CBE", "#D7301F", "skyblue","seashell3", "plum1")
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

ggplot() + 
  geom_point(data=marker, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=1) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(marker$colours))),
                     labels=rev(names(table(marker$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=3) +
  geom_text_repel(data=genes_to_plot,
                  aes(x=avg_log2FC, y=-log10(p_val_adj), label = row.names(genes_to_plot)), 
                  color="black", fontface = 'bold',size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("MG4 vs others")+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15, face = "bold")) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-2, 2, 1), limits=c(-2, 2))+
  NoLegend()

ggsave("Volcano_subcluster_4_vs_others.pdf", plot = last_plot(), device = "pdf", path = "~/Desktop/data_analysis/LG72/integration_5genotypes/MG_only",
       scale = 0.6, width = 9, height = 9, units = c("in"),
       dpi = 600, limitsize = FALSE)
