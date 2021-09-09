#################################################################################
# Generation of Figure 1E-I from Sayed, Kodama, Fan, et al. 2021 
# Human microglia analysis post-processing of 46 human sample
# For pre-processing steps from raw data, refer to Figure 1 pre-processing codes

# by Lay Kodama
#################################################################################

library(ComplexHeatmap)
library(dplyr)
devtools::install_github("mw201608/msigdb")
devtools::install_github("mw201608/GOtest")
library(msigdb)
library(GOtest)
library(ggplot2)
library(ggrepel)

# Plot DEG overlap between cell types (Figure 1E)=======
# Import DEG data (DEGs identified using Step 4 of pre-processing code for Figure 1)
setwd()
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

celltypes <- c("astrocyte","astrocyte","EC","EC","excitatory","excitatory","inhibitory","inhibitory","microglia","microglia",
               "OG","OG","OPC","OPC")
sex <- c(rep(c("female","male"),7))

deg_list <- list()
for (i in 1:length(celltypes)){
  a <- data.frame(myfiles[[i]])
  a$cell_type <- celltypes[i]
  a$sex <- sex[i]
  
  deg_list[[i]] <- a
}

deg_list <- do.call(rbind,deg_list)
deg_list <- subset(deg_list, p_val_adj < 0.05)

write.csv(deg_list, "DEG_full_list.csv")

genes <- data.frame(genes = unique(deg_list$X))
genes$astrocyte <- 0

ast.f <- subset(deg_list, cell_type == "astrocyte" & sex == "female")
ast.m <- subset(deg_list, cell_type == "astrocyte" & sex == "male")
ast.both <- subset(ast.f, X %in% ast.m$X)
genes$astrocyte[genes$genes%in%ast.f$X] <- 1
genes$astrocyte[genes$genes%in%ast.m$X] <- 2
genes$astrocyte[genes$genes%in%ast.both$X] <- 3

genes$EC <- 0
EC.f <- subset(deg_list, cell_type == "EC" & sex == "female")
EC.m <- subset(deg_list, cell_type == "EC" & sex == "male")
EC.both <- subset(EC.f, X %in% EC.m$X)
genes$EC[genes$genes%in%EC.f$X] <- 1
genes$EC[genes$genes%in%EC.m$X] <- 2
genes$EC[genes$genes%in%EC.both$X] <- 3

genes$EN <- 0
EN.f <- subset(deg_list, cell_type == "excitatory" & sex == "female")
EN.m <- subset(deg_list, cell_type == "excitatory" & sex == "male")
EN.both <- subset(EN.f, X %in% EN.m$X)
genes$EN[genes$genes%in%EN.f$X] <- 1
genes$EN[genes$genes%in%EN.m$X] <- 2
genes$EN[genes$genes%in%EN.both$X] <- 3

genes$IN <- 0
IN.f <- subset(deg_list, cell_type == "inhibitory" & sex == "female")
IN.m <- subset(deg_list, cell_type == "inhibitory" & sex == "male")
IN.both <- subset(IN.f, X %in% IN.m$X)
genes$IN[genes$genes%in%IN.f$X] <- 1
genes$IN[genes$genes%in%IN.m$X] <- 2
genes$IN[genes$genes%in%IN.both$X] <- 3

genes$Mg <- 0
Mg.f <- subset(deg_list, cell_type == "microglia" & sex == "female")
Mg.m <- subset(deg_list, cell_type == "microglia" & sex == "male")
Mg.both <- subset(Mg.f, X %in% Mg.m$X)
genes$Mg[genes$genes%in%Mg.f$X] <- 1
genes$Mg[genes$genes%in%Mg.m$X] <- 2
genes$Mg[genes$genes%in%Mg.both$X] <- 3

genes$Oligo <- 0
Oligo.f <- subset(deg_list, cell_type == "OG" & sex == "female")
Oligo.m <- subset(deg_list, cell_type == "OG" & sex == "male")
Oligo.both <- subset(Oligo.f, X %in% Oligo.m$X)
genes$Oligo[genes$genes%in%Oligo.f$X] <- 1
genes$Oligo[genes$genes%in%Oligo.m$X] <- 2
genes$Oligo[genes$genes%in%Oligo.both$X] <- 3

genes$OPC <- 0
OPC.f <- subset(deg_list, cell_type == "OPC" & sex == "female")
OPC.m <- subset(deg_list, cell_type == "OPC" & sex == "male")
OPC.both <- subset(OPC.f, X %in% OPC.m$X)
genes$OPC[genes$genes%in%OPC.f$X] <- 1
genes$OPC[genes$genes%in%OPC.m$X] <- 2
genes$OPC[genes$genes%in%OPC.both$X] <- 3

Heatmap(as.matrix(genes[,-1]),row_labels = F,col=c("white","salmon","dodgerblue","purple"),
        column_order = c("astrocyte","EN","IN","Mg","OPC","Oligo","EC"),
        show_row_dend = F)


# Gene set enrichment analysis using MSIGDB database (Figure 1H,I) ====
universe=curated.genesets(c('HGNC_universe'))$Gene

deg_list$FC <- "up" 
deg_list$FC[deg_list$avg_logFC < 0] <- "down" 
deg_list$category <- paste(deg_list$cell_type, deg_list$sex, deg_list$FC, sep="_")

gene_test <- deg_list[,c(1,ncol(deg_list))]
colnames(gene_test)[1] <- "Gene"
result = msigdb.gsea(x=gene_test, query.population=universe, genesets=c('C5.BP','C5.CC', 'C5.MF', 'C2.CP'), background='query', name.x='DEGs')
head(result)

# create barplots of top 5 pathways per comparison for microglia:
top5<-result %>% 
  arrange(P.adj) %>% 
  group_by(DEGs) %>% slice(1:5)

df <- select(top5, DEGs, MSigDB, Overlap.Size, MSigDB.Size, P.adj)
df$direction <- gsub("^[^.]*_[^.]*_([^.]*)","\\1",df$DEGs)
df$celltype <- gsub("(^[^.]*)_[^.]*_[^.]*","\\1",df$DEGs)
df$sex <- gsub("^[^.]*_([^.]*)_[^.]*","\\1",df$DEGs)

colnames(df) <- c("name", "description", "genes_in_data", "genes_in_gsea", "FDR","direction","celltype","sex")

top <- df %>% group_by(name) %>% top_n(n=1,wt=-log10(FDR))
top$sign <- -1
top$sign[top$direction=='up'] <- 1
top$FDR_sign <- -log10(top$FDR) * top$sign
top$description_celltype <- paste(top$description, top$celltype, sep = "_")

mg_df <- subset(df, celltype == "microglia")

ggplot(subset(mg_df, sex == "female"), 
       aes(x = reorder(description, -FDR), y = -log10(FDR), fill = direction)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(x = "", y = "-Log10(FDR)", title = "R47H vs CV in Females") +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("down"="blue","up"="red"))

ggplot(subset(mg_df, sex == "male"), 
       aes(x = reorder(description, -FDR), y = -log10(FDR), fill = direction)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  labs(x = "", y = "-Log10(FDR)", title = "R47H vs CV in Males") +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("down"="blue","up"="red"))

# volcano plot of DEGs (Figure 1F,) =======
mg_deg_list <- subset(deg_list, cell_type == "microglia" & sex == "female")

setwd()
dam <- read.csv("DAM.csv") # DAM signature genes from published literature
dam_sig <- subset(dam, log.pvalue. > -log10(0.05))
dam_sig$cap <- toupper(dam_sig$gene)
dam_sig$sign <- "down"
dam_sig$sign[dam_sig$foldchange > 0] <- "up"

mg_deg_list_up <- subset(mg_deg_list, FC =="up")
mg_deg_list_up$DAM <- "no"
dam_up <- subset(dam_sig, sign =="up")
mg_deg_list_up$DAM[match(dam_up$cap, mg_deg_list_up$X)] <- "yes"

mg_deg_list_dn <- subset(mg_deg_list, FC =="down")
mg_deg_list_dn$DAM <- "no"
dam_dn <- subset(dam_sig, sign =="down")
mg_deg_list_dn$DAM[match(dam_dn$cap, mg_deg_list_dn$X)] <- "yes"

mg_deg_list <- rbind(mg_deg_list_dn, mg_deg_list_up)

ggplot(mg_deg_list, aes(x = avg_logFC, y = -log10(p_val_adj), color = FC)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("down" = "blue", "up" = "red")) +
  geom_text_repel(subset(mg_deg_list, DAM == "yes"),mapping = aes(label = X)) +
  labs(y = "-Log10(FDR)", x = "Log2(R47H vs CV)", title = "R47H vs CV in Females") +
  theme(legend.position="none")

mg_deg_list <- subset(deg_list, cell_type == "microglia" & sex == "male")

mg_deg_list_up <- subset(mg_deg_list, FC =="up")
mg_deg_list_up$DAM <- "no"
dam_up <- subset(dam_sig, sign =="up")
mg_deg_list_up$DAM[match(dam_up$cap, mg_deg_list_up$X)] <- "yes"

mg_deg_list_dn <- subset(mg_deg_list, FC =="down")
mg_deg_list_dn$DAM <- "no"
dam_dn <- subset(dam_sig, sign =="down")
mg_deg_list_dn$DAM[match(dam_dn$cap, mg_deg_list_dn$X)] <- "yes"

mg_deg_list <- rbind(mg_deg_list_dn, mg_deg_list_up)

ggplot(mg_deg_list, aes(x = avg_logFC, y = -log10(p_val_adj), color = FC)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = c("down" = "blue", "up" = "red")) +
  geom_text_repel(subset(mg_deg_list, DAM == "yes"),mapping = aes(label = X)) +
  labs(y = "-Log10(FDR)", x = "Log2(R47H vs CV)", title = "R47H vs CV in Males") +
  theme(legend.position="none")


