#################################################################################
# Generation of Figure 5C from Sayed, Kodama, Fan, et al. 2021 
# Comparison of IPA pathway enrichment from different modalities (Figure 5C)
# in vitro, human single-nuclei microglia (MG4), mouse single-cell microglia (cluster 2)

# by Lay Kodama
#################################################################################

setwd()
invitro <- read.csv("R47HvsWT_Tau_USRs.csv",header=T) 
deepseq <- read.csv("IPA_deepseq_trem2signaling.csv",header=T)
human <- read.csv("MG4_marker_IPA.csv",header=T) # from new combined human data analysis

upstream <- rbind(data.frame(upstream = invitro$Upstream.Regulator), data.frame(upstream = deepseq$Upstream.Regulator))
upstream <- rbind(upstream, data.frame(upstream = human$Upstream.Regulator))
upstream <- rbind(upstream, data.frame(upstream = mk$Upstream.Regulator))

upstream <- subset(upstream,!duplicated(upstream))

upstream$invitro <- invitro$Activation.z.score[match(upstream$upstream, invitro$Upstream.Regulator)]            
upstream$deepseq <- deepseq$Activation.z.score[match(upstream$upstream, deepseq$Upstream.Regulator)]            
upstream$human <- human$Activation.z.score[match(upstream$upstream, human$Upstream.Regulator)]            

upstream[is.na(upstream)] <- 0

row.names(upstream) <- upstream$upstream
upstream <- upstream[,-1]

heatmap(as.matrix(upstream))
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 3, 6), c("white", "salmon", "red"))
col_fun(seq(0, 6))
Heatmap(as.matrix(upstream),  col = col_fun, column_order = order(as.numeric(gsub("column", "", colnames(upstream)))))

