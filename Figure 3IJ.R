#################################################################################
# Generation of Figure 3I,J from Sayed, Kodama, Fan, et al. 2021 
# R47H vs CV mouse hippocampal bulk RNA seq WGCNA 
# Preprocessing to generation of module trait heatmap
# Workflow adapted from tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

# by Lay Kodama
#################################################################################

library(DESeq2)
library(ggplot2)
library("genefilter")
library(WGCNA)

#=====================================================================================
#
#  STEP 0: Data pre-processing for most variable genes
#
#=====================================================================================

setwd()

#Import counts data for bulk sequencing
Seq = read.csv('counts_R47H_WT_MF.csv', header = TRUE)

Sq = Seq[rowSums(Seq[2:23]) > 15,]
Sq <- Sq[rowSums(Sq==0)<4, ]

#Import sample table
sampleTable <- read.csv("sample_table_R47H_WT_MF.csv", header = TRUE,sep=",")

#Construct data matrix that can be input into DE function:
Data <- Sq[2:23]
dds = DESeqDataSetFromMatrix(countData = Data, colData = sampleTable, design = ~genotype + sex)
#set control as CV Male
dds$control <- relevel(dds$genotype, "WT", ddes$sex, "M")
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
vsd <- vst(dds, blind=FALSE)

vsd_counts <- as.data.frame(assay(vsd))

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 5000) # filter top 5000 most variable genes for WGCNA
mat  <- as.data.frame(assay(vsd)[ topVarGenes, ])

#Export top most variable genes as CSV
genelist = Sq$gene
top_genes = genelist[ topVarGenes ]
mat$genes<-top_genes

write.csv(mat, "top_vsd_R47H_WT_MF.csv")

#=====================================================================================
#
#  STEP 1: Data input for WGCNA
#
#=====================================================================================

options(stringsAsFactors = FALSE);
Data = read.csv("top_vsd_R47H_WT_MF.csv");

# remove extraneous columns and transpose data
Data1 = Data[,-1]
datExpr = as.data.frame(t(Data1));
names(datExpr) = Data$genes;
rownames(datExpr) = names(Data1);

#check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#look for any obvious outliers
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 55, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 55, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#read in trait data
traitData = read.csv("wgcna_sample_table_R47H_WT_MF.csv");
dim(traitData)
names(traitData)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
traitRows = match(Samples, traitData$sample);
datTraits = traitData[traitRows, -1];
rownames(datTraits) = traitData[traitRows, 1];

collectGarbage();

# how clinical traits relate to the samples

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# save the relevant data files
save(datExpr, datTraits, file = "Faten-dataInput_R47H_WT_MF.RData")

#=====================================================================================
#
#  STEP 2: Module identification
#
#=====================================================================================

allowWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Faten-dataInput_R47H_WT_MF.RData");
#The variable lnames contains the names of loaded variables.
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 14;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Faten-stepByStep_R47H_WT_MF.RData")

MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#=====================================================================================
#
#  STEP 3: Module-trait heatmap generation
#
#=====================================================================================

lnames = load(file = "Faten-dataInput_R47H_WT_MF.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Faten-stepByStep_R47H_WT_MF.RData");
lnames

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

# Define variable genotype containing the genotype column of datTrait
genotype = as.data.frame(datTraits$Genotype);
names(genotype) = "genotype"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, genotype, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(genotype), sep="");
names(GSPvalue) = paste("p.GS.", names(genotype), sep="");

#module with high significance for genotype (most red)
module = "lightcyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for genotype",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

names(datExpr)

names(datExpr)[moduleColors=="lightcyan"]

# Create the starting data frame
geneInfo0 = data.frame(genes = names(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, genotype, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.genotype));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")



