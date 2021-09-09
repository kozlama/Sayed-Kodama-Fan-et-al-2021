#################################################################################
# Generation of Figure 5 from Sayed, Kodama, Fan, et al. 2021 
# In vitro mouse primary microglia R47H/+ vs Trem2 +/+
# treated with or without tau fibrils and with or without MK-2206 (Akt inhibitor)

# Code adapted from DESeq2 tutorial: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# by Gillian Carling
#################################################################################

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(VennDiagram)
library(ggrepel)
library(data.table)
library(biomaRt)
library(reshape2)
library(stringr)
library(circlize)
library(ComplexHeatmap)

setwd()

# DESeq - input data 
data <- read.csv("counts_new.csv", header=T)
data <- subset(data, !duplicated(data$gene_name))
row.names(data) <- data$gene_name
data <- data[,-1]
write.csv(data, "counts_new.csv")

data <- read.csv("counts_new.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,]
sample <- colnames(data)
genotype <- c(rep("WT",6),rep("R47H",12))
drug <- c(rep("Vehicle",12),rep("MK2206",6))
tau <- c(rep("CTL",3),rep("Tau",3),rep("CTL",3),rep("Tau",3),rep("CTL",3),rep("Tau",3))
genotype_drug_tau <- paste(genotype,drug,tau,sep="_")
meta <- data.frame(sample=sample, genotype=genotype, drug=drug, tau=tau, genotype_drug_tau=genotype_drug_tau)

all(colnames(filtered) %in% meta$sample)

# DESeq Matrix
dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~genotype_drug_tau)
rld <- rlog(dds, blind = T)

plotPCA(rld, intgroup = "genotype_drug_tau", ntop = 500)+
  geom_text_repel(aes(label = genotype_drug_tau))+
  theme_classic()

rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
pheatmap(rld_cor)

design <- ~ genotype_drug_tau
dds <- DESeq(dds)

# DESeq - Tau fibril treatment condition Only ####
data <- read.csv("counts_tau.csv", header=T, row.names=1)
filtered = data[rowSums(data)>15,]
sample <- colnames(data)
genotype <- c(rep("WT",3),rep("R47H",6))
drug <- c(rep("Vehicle",6),rep("MK2206",3))
genotype_drug <- paste(genotype,drug,sep="_")
meta <- data.frame(sample=sample, genotype=genotype, drug=drug, genotype_drug=genotype_drug)

all(colnames(filtered) %in% meta$sample)

dds <- DESeqDataSetFromMatrix(countData = filtered, colData = meta, design = ~genotype_drug)
rld <- rlog(dds, blind = T)

plotPCA(rld, intgroup = "genotype_drug", ntop = 500)+
  geom_text_repel(aes(label = genotype_drug))+
  theme_classic()

rld_matrix <- assay(rld)
rld_cor <- cor(rld_matrix)
pheatmap(rld_cor)

design <- ~ genotype_drug
dds <- DESeq(dds)

# Determining DEGs for different pairwise conditions (Figure 5A) ======
# Pairwise Comparisons - Full Data ####
# R47H vs. WT CTL
contrast_R47HvsWT <- c("genotype_drug_tau","R47H_Vehicle_CTL","WT_Vehicle_CTL")
res_R47HvsWT_unshrunken <- results(dds,contrast=contrast_R47HvsWT,alpha=0.05)
res_R47HvsWT <- lfcShrink(dds,contrast=contrast_R47HvsWT,res=res_R47HvsWT_unshrunken)
write.csv(res_R47HvsWT, "DE_R47HvsWT_CTL.csv")

# Pairwise Comparisons - Tau Only ####
# R47H vs. WT Tau
contrast_R47HvsWT_tau <- c("genotype_drug","R47H_Vehicle","WT_Vehicle")
res_R47HvsWT_tau_unshrunken <- results(dds,contrast=contrast_R47HvsWT_tau,alpha=0.05)
res_R47HvsWT_tau <- lfcShrink(dds,contrast=contrast_R47HvsWT_tau,res=res_R47HvsWT_tau_unshrunken)
write.csv(res_R47HvsWT_tau, "DE_R47HvsWT_Tau.csv")

# MK2206 vs. WT Tau
contrast_MK2206vsWT_tau <- c("genotype_drug","R47H_MK2206","WT_Vehicle")
res_MK2206vsWT_tau_unshrunken <- results(dds,contrast=contrast_MK2206vsWT_tau,alpha=0.05)
res_MK2206vsWT_tau <- lfcShrink(dds,contrast=contrast_MK2206vsWT_tau,res=res_MK2206vsWT_tau_unshrunken)
write.csv(res_MK2206vsWT_tau, "DE_MK2206vsWT_Tau.csv")

# MK2206 vs. R47H Vehicle Tau
contrast_MK2206vsR47H_tau <- c("genotype_drug","R47H_MK2206","R47H_Vehicle")
res_MK2206vsR47H_tau_unshrunken <- results(dds,contrast=contrast_MK2206vsR47H_tau,alpha=0.05)
res_MK2206vsR47H_tau <- lfcShrink(dds,contrast=contrast_MK2206vsR47H_tau,res=res_MK2206vsR47H_tau_unshrunken)
write.csv(res_MK2206vsR47H_tau, "DE_MK2206vsR47H_Tau.csv")

# Reversed by MK2206 ####
MK2206vsVehicle <- read.csv("DE_MK2206vsR47H_Tau.csv",header=T,row.names=1)
R47HvsWT <- read.csv("DE_R47HvsWT_Tau.csv",header=T,row.names=1)

MK2206.up <- subset(MK2206vsVehicle, log2FoldChange > 0.5 & padj < 0.05)
MK2206.dn <- subset(MK2206vsVehicle, log2FoldChange < -0.5 & padj < 0.05)
R47H.up <- subset(R47HvsWT, log2FoldChange > 0.5 & padj < 0.05)
R47H.dn <- subset(R47HvsWT, log2FoldChange < -0.5 & padj < 0.05)

reversed.up <- subset(R47H.up, row.names(R47H.up) %in% row.names(MK2206.dn))
reversed.dn <- subset(R47H.dn, row.names(R47H.dn) %in% row.names(MK2206.up))

write.csv(reversed.up, "Reversed_R47HvsWT_MK2206_Up.csv")
write.csv(reversed.dn, "Reversed_R47HvsWT_MK2206_Down.csv")

# Heatmap (Figure 5D) =====
MK2206.fixed <- read.csv("Reversed_R47HvsWT_MK2206_ALL.csv",header=T,row.names=1)

colors <- rev(brewer.pal(10,"RdBu"))
counts <- counts(dds,normalized = T)
annotations <- as.data.frame(colData(dds)["genotype_drug"])

overlap <- subset(counts, row.names(counts) %in% row.names(MK2206.fixed))
write.csv(overlap, "Reversed_R47HvsWT_MK2206_Counts.csv")

pheatmap(overlap,
         show_rownames = F,
         annotation = annotations,
         color = colors,
         scale = "row")

# Venn Diagram (Figure 5A) ####
R47HvsWT_CTL <- read.csv("DE_R47HvsWT_CTL.csv",header=T,row.names=1)
R47HvsWT_Tau <- read.csv("DE_R47HvsWT_Tau.csv",header=T,row.names=1)

CTL.up <- subset(R47HvsWT_CTL, log2FoldChange > 0.5 & padj < 0.05)
CTL.dn <- subset(R47HvsWT_CTL, log2FoldChange < -0.5 & padj < 0.05)
Tau.up <- subset(R47HvsWT_Tau, log2FoldChange > 0.5 & padj < 0.05)
Tau.dn <- subset(R47HvsWT_Tau, log2FoldChange < -0.5 & padj < 0.05)
write.csv(Tau.up, "R47HvsWT_Tau_Upregulated.csv")
write.csv(Tau.dn, "R47HvsWT_Tau_Downregulated.csv")

overlap.up <- subset(CTL.up, row.names(CTL.up) %in% row.names(Tau.up))
overlap.dn <- subset(CTL.dn, row.names(CTL.dn) %in% row.names(Tau.dn))

unique.CTL.up <- subset(CTL.up, !(row.names(CTL.up) %in% row.names(overlap.up)))
unique.CTL.down <- subset(CTL.dn,!(row.names(CTL.dn) %in% row.names(overlap.dn)))
unique.Tau.up <- subset(Tau.up, !(row.names(Tau.up) %in% row.names(overlap.up)))
unique.Tau.down <- subset(Tau.dn,!(row.names(Tau.dn) %in% row.names(overlap.dn)))

# KEGG Pathway Analysis (Figuure 5E) ####
# R47H vs. WT (Tau)
KEGG <- read.csv("KEGG_R47HvsWT_Tau_ALL.csv",header=T)
colnames(KEGG) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR", "Type")
KEGG$FDR <- as.numeric(KEGG$FDR)
KEGG$logP <- -log10(KEGG$FDR)
KEGG$enrich <- paste(KEGG$genes_in_data, "/", KEGG$genes_in_gsea, sep=" ")
write.csv(KEGG, "KEGG_R47HvsWT_Tau_Edited.csv")
KEGG <- read.csv("KEGG_R47HvsWT_Tau_Edited.csv",header=T)
KEGG$Name <-  gsub("KEGG_", "", KEGG$name)
KEGG$Name <- gsub("*_", " ", KEGG$Name)
KEGG$Name <-  factor(KEGG$Name, levels=rev(KEGG$Name))

ggplot(data=KEGG, aes(x=Name  , y=logP, fill=Type)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, alpha=0.8) +
  coord_flip() + 
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"))+
  theme(aspect.ratio = 1.5)

# MK2206 Reversed DEGs
Rev <- read.csv("Reversed_R47HvsWT_MK2206_sorted.csv",header=T)
colnames(Rev) <- c("name", "genes_in_gsea", "description", "genes_in_data", "k/K", "p-value", "FDR")
Rev$FDR <- as.numeric(Rev$FDR)
Rev$logP <- -log10(Rev$FDR)
Rev$enrich <- paste(Rev$genes_in_data, "/", Rev$genes_in_gsea, sep=" ")
Rev <- Rev[order(-Rev$logP),]
Rev$Name <-  gsub("KEGG_", "", Rev$name)
Rev$Name <- gsub("*_", " ", Rev$Name)
Rev$Name <-  factor(Rev$Name, levels=rev(Rev$Name))

ggplot(data=Rev, aes(x=Name  , y=logP)) +
  theme_classic() +
  ylab("-Log(FDR)") + xlab(NULL) +
  geom_bar(stat="Identity",  width=0.7, fill="green4", alpha=0.8) +
  geom_text(aes(label=enrich), vjust=0.4, 
            hjust=0, size=4, color="black", 
            stat="Identity", y=0.05*max(Rev$logP)) +
  coord_flip() + 
  theme(aspect.ratio = 1.5)

# Upstream Regulators (Figure 5C) ####
USR <- read.csv("USR_R47HvsWT_Tau.csv",header=T)
USR$p.value <- as.numeric(USR$p.value)
USR$logP <- -log10(USR$p.value)
USR <- USR[order(-USR$z.score),]
USR$name <- factor(USR$name, levels = USR$name[order(USR$z.score)])

ggplot(data=USR, aes(x=name, y=z.score)) +
  theme_classic() +
  ylab("Predicted Activation Z-Score") + xlab(NULL) +
  geom_bar(aes(fill=logP), stat="Identity",  width=0.7, alpha=0.8) +
  coord_flip() + 
  scale_fill_distiller(palette = "Reds") +
  theme(aspect.ratio = 1.5)

# MagPix Graphs (Figure 5G) ####
magpix <- read.csv("MagPix_MK2206_Reversed.csv", header=T)
df <- melt(magpix, id = "X")
df$genotype <- c(rep("WT_Vehicle",18),rep("R47H_Vehicle",18),rep("R47H_MK2206",18))
df$genotype <- factor(df$genotype, levels = c("WT_Vehicle", "R47H_Vehicle", "R47H_MK2206"), ordered = TRUE)

# GCSF
GCSF <- subset(df, X=="GCSF")
ggplot(GCSF,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("GCSF")

# CCL2
CCL2 <- subset(df, X=="CCL2")
ggplot(CCL2,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("CCL2")

# CCL3
CCL3 <- subset(df, X=="CCL3")
ggplot(CCL3,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("CCL3")

# CXCL9
CXCL9 <- subset(df, X=="CXCL9")
ggplot(CXCL9,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("CXCL9")

# IL5
IL5 <- subset(df, X=="IL5")
ggplot(IL5,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("IL5")

# IL6
IL6 <- subset(df, X=="IL6")
ggplot(IL6,aes(x=genotype,y=value,fill=genotype))+
  geom_bar(stat="summary")+
  ylab("Protein Concentration (pg/mL)") +
  geom_point()+
  theme_classic()+
  scale_fill_manual(values = c("WT_Vehicle" = "grey", "R47H_Vehicle" = "salmon", "R47H_MK2206" = "chartreuse3"))+
  ggtitle("IL6")


