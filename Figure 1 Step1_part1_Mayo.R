###############################################################################################
# Pre-processing for data used to generate Figure 1 from Sayed, Kodama, Fan, et al. 2021 
# Total of 55 human AD R47H vs CV samples and non-AD samples
# This script is: STEP 1 of 5 (PART 1 of 2) - 34 human AD R47H vs CV samples from Mayo Clinic

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(ggplot2)
library(DoubletFinder)

#set working directory ====
setwd()

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
Human_R47H_1.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_1")
Human_R47H_1 <- CreateSeuratObject(counts = Human_R47H_1.counts, project = "R47H_E2E4_M", min.cells = 3, min.features = 200)
Human_R47H_1[["Condition"]] = c('R47H_E2E4_M')
Human_R47H_1[["Condition_1"]] = c('R47H_E2E4')
Human_R47H_1[["Condition_2"]] = c('R47H_M')
Human_R47H_1[["Condition_3"]] = c('E2E4_M')
Human_R47H_1[["TREM2"]] = c('R47H')
Human_R47H_1[["Dx"]] = c('AD')
Human_R47H_1[["LBD"]] = c('N/A')
Human_R47H_1[["Braak"]] = c('5')
Human_R47H_1[["Thal"]] = c('5')
Human_R47H_1[["TDP.43"]] = c('0')
Human_R47H_1[["ClinicalDx"]] = c('AD')
Human_R47H_1[["APOE"]] = c('E2E4')
Human_R47H_1[["Age"]] = c('91')
Human_R47H_1[["Sex"]] = c('M')
rm(Human_R47H_1.counts)

#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Human_R47H_1[["percent.mt"]] <- PercentageFeatureSet(object = Human_R47H_1, pattern = "^MT-") #recognize mitochondrial transcripts

Human_R47H_2.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_2")
Human_R47H_2 <- CreateSeuratObject(counts = Human_R47H_2.counts, project = "WT_E2E4_M", min.cells = 3, min.features = 200)
Human_R47H_2[['percent.mt']] <- PercentageFeatureSet(Human_R47H_2, pattern = "^MT-")
Human_R47H_2[["Condition"]] = c('WT_E2E4_M')
Human_R47H_2[["Condition_1"]] = c('WT_E2E4')
Human_R47H_2[["Condition_2"]] = c('WT_M')
Human_R47H_2[["Condition_3"]] = c('E2E4_M')
Human_R47H_2[["TREM2"]] = c('WT')
Human_R47H_2[["Dx"]] = c('AD')
Human_R47H_2[["LBD"]] = c('N/A')
Human_R47H_2[["Braak"]] = c('6')
Human_R47H_2[["Thal"]] = c('5')
Human_R47H_2[["TDP.43"]] = c('N/A')
Human_R47H_2[["ClinicalDx"]] = c('AD')
Human_R47H_2[["APOE"]] = c('E2E4')
Human_R47H_2[["Age"]] = c('83')
Human_R47H_2[["Sex"]] = c('M')
rm(Human_R47H_2.counts)

Human_R47H_3.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_3")
Human_R47H_3 <- CreateSeuratObject(counts = Human_R47H_3.counts, project = "R47H_E3E3_F_1", min.cells = 3, min.features = 200)
Human_R47H_3[['percent.mt']] <- PercentageFeatureSet(Human_R47H_3, pattern = "^MT-")
Human_R47H_3[["Condition"]] = c('R47H_E3E3_F')
Human_R47H_3[["Condition_1"]] = c('R47H_E3E3')
Human_R47H_3[["Condition_2"]] = c('R47H_F')
Human_R47H_3[["Condition_3"]] = c('E3E3_F')
Human_R47H_3[["TREM2"]] = c('R47H')
Human_R47H_3[["Dx"]] = c('AD/DLBD')
Human_R47H_3[["LBD"]] = c('DLBD')
Human_R47H_3[["Braak"]] = c('5')
Human_R47H_3[["Thal"]] = c('5')
Human_R47H_3[["TDP.43"]] = c('0')
Human_R47H_3[["ClinicalDx"]] = c('AD')
Human_R47H_3[["APOE"]] = c('E3E3')
Human_R47H_3[["Age"]] = c('81')
Human_R47H_3[["Sex"]] = c('F')
rm(Human_R47H_3.counts)

Human_R47H_4.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_4")
Human_R47H_4 <- CreateSeuratObject(counts = Human_R47H_4.counts, project = "WT_E3E3_F_1", min.cells = 3, min.features = 200)
Human_R47H_4[['percent.mt']] <- PercentageFeatureSet(Human_R47H_4, pattern = "^MT-")
Human_R47H_4[["Condition"]] = c('WT_E3E3_F')
Human_R47H_4[["Condition_1"]] = c('WT_E3E3')
Human_R47H_4[["Condition_2"]] = c('WT_F')
Human_R47H_4[["Condition_3"]] = c('E3E3_F')
Human_R47H_4[["TREM2"]] = c('WT')
Human_R47H_4[["Dx"]] = c('AD/DLBD')
Human_R47H_4[["LBD"]] = c('DLBD')
Human_R47H_4[["Braak"]] = c('5')
Human_R47H_4[["Thal"]] = c('5')
Human_R47H_4[["TDP.43"]] = c('0')
Human_R47H_4[["ClinicalDx"]] = c('AD')
Human_R47H_4[["APOE"]] = c('E3E3')
Human_R47H_4[["Age"]] = c('85')
Human_R47H_4[["Sex"]] = c('F')
rm(Human_R47H_4.counts)

Human_R47H_5.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_5")
Human_R47H_5 <- CreateSeuratObject(counts = Human_R47H_5.counts, project = "R47H_E3E3_M_1", min.cells = 3, min.features = 200)
Human_R47H_5[['percent.mt']] <- PercentageFeatureSet(Human_R47H_5, pattern = "^MT-")
Human_R47H_5[["Condition"]] = c('R47H_E3E3_M')
Human_R47H_5[["Condition_1"]] = c('R47H_E3E3')
Human_R47H_5[["Condition_2"]] = c('R47H_M')
Human_R47H_5[["Condition_3"]] = c('E3E3_M')
Human_R47H_5[["TREM2"]] = c('R47H')
Human_R47H_5[["Dx"]] = c('AD/CAA')
Human_R47H_5[["LBD"]] = c('N/A')
Human_R47H_5[["Braak"]] = c('5')
Human_R47H_5[["Thal"]] = c('4')
Human_R47H_5[["TDP.43"]] = c('0')
Human_R47H_5[["ClinicalDx"]] = c('FTD')
Human_R47H_5[["APOE"]] = c('E3E3')
Human_R47H_5[["Age"]] = c('65')
Human_R47H_5[["Sex"]] = c('M')
rm(Human_R47H_5.counts)

Human_R47H_6.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_6")
Human_R47H_6 <- CreateSeuratObject(counts = Human_R47H_6.counts, project = "WT_E3E3_M_1", min.cells = 3, min.features = 200)
Human_R47H_6[['percent.mt']] <- PercentageFeatureSet(Human_R47H_6, pattern = "^MT-")
Human_R47H_6[["Condition"]] = c('WT_E3E3_M')
Human_R47H_6[["Condition_1"]] = c('WT_E3E3')
Human_R47H_6[["Condition_2"]] = c('WT_M')
Human_R47H_6[["Condition_3"]] = c('E3E3_M')
Human_R47H_6[["TREM2"]] = c('WT')
Human_R47H_6[["Dx"]] = c('AD/CAA')
Human_R47H_6[["LBD"]] = c('N/A')
Human_R47H_6[["Braak"]] = c('6')
Human_R47H_6[["Thal"]] = c('5')
Human_R47H_6[["TDP.43"]] = c('0')
Human_R47H_6[["ClinicalDx"]] = c('AD_v_FTD')
Human_R47H_6[["APOE"]] = c('E3E3')
Human_R47H_6[["Age"]] = c('65')
Human_R47H_6[["Sex"]] = c('M')
rm(Human_R47H_6.counts)

Human_R47H_7.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_7")
Human_R47H_7 <- CreateSeuratObject(counts = Human_R47H_7.counts, project = "R47H_E3E3_F_2", min.cells = 3, min.features = 200)
Human_R47H_7[['percent.mt']] <- PercentageFeatureSet(Human_R47H_7, pattern = "^MT-")
Human_R47H_7[["Condition"]] = c('R47H_E3E3_F')
Human_R47H_7[["Condition_1"]] = c('R47H_E3E3')
Human_R47H_7[["Condition_2"]] = c('R47H_F')
Human_R47H_7[["Condition_3"]] = c('E3E3_F')
Human_R47H_7[["TREM2"]] = c('R47H')
Human_R47H_7[["Dx"]] = c('AD')
Human_R47H_7[["LBD"]] = c('N/A')
Human_R47H_7[["Braak"]] = c('6')
Human_R47H_7[["Thal"]] = c('5')
Human_R47H_7[["TDP.43"]] = c('0')
Human_R47H_7[["ClinicalDx"]] = c('CBD')
Human_R47H_7[["APOE"]] = c('E3E3')
Human_R47H_7[["Age"]] = c('62')
Human_R47H_7[["Sex"]] = c('F')
rm(Human_R47H_7.counts)

Human_R47H_8.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_8")
Human_R47H_8 <- CreateSeuratObject(counts = Human_R47H_8.counts, project = "WT_E3E3_F_2", min.cells = 3, min.features = 200)
Human_R47H_8[['percent.mt']] <- PercentageFeatureSet(Human_R47H_8, pattern = "^MT-")
Human_R47H_8[["Condition"]] = c('WT_E3E3_F')
Human_R47H_8[["Condition_1"]] = c('WT_E3E3')
Human_R47H_8[["Condition_2"]] = c('WT_F')
Human_R47H_8[["Condition_3"]] = c('E3E3_F')
Human_R47H_8[["TREM2"]] = c('WT')
Human_R47H_8[["Dx"]] = c('AD')
Human_R47H_8[["LBD"]] = c('N/A')
Human_R47H_8[["Braak"]] = c('5')
Human_R47H_8[["Thal"]] = c('5')
Human_R47H_8[["TDP.43"]] = c('0')
Human_R47H_8[["ClinicalDx"]] = c('AD')
Human_R47H_8[["APOE"]] = c('E3E3')
Human_R47H_8[["Age"]] = c('67')
Human_R47H_8[["Sex"]] = c('F')
rm(Human_R47H_8.counts)


Human_R47H_9.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_9")
Human_R47H_9 <- CreateSeuratObject(counts = Human_R47H_9.counts, project = "R47H_E3E3_M_2", min.cells = 3, min.features = 200)
Human_R47H_9[['percent.mt']] <- PercentageFeatureSet(Human_R47H_9, pattern = "^MT-")
Human_R47H_9[["Condition"]] = c('R47H_E3E3_M')
Human_R47H_9[["Condition_1"]] = c('R47H_E3E3')
Human_R47H_9[["Condition_2"]] = c('R47H_M')
Human_R47H_9[["Condition_3"]] = c('E3E3_M')
Human_R47H_9[["TREM2"]] = c('R47H')
Human_R47H_9[["Dx"]] = c('AD/DLBD')
Human_R47H_9[["LBD"]] = c('DLBD')
Human_R47H_9[["Braak"]] = c('5.5')
Human_R47H_9[["Thal"]] = c('5')
Human_R47H_9[["TDP.43"]] = c('0')
Human_R47H_9[["ClinicalDx"]] = c('PDD_v_DLB')
Human_R47H_9[["APOE"]] = c('E3E3')
Human_R47H_9[["Age"]] = c('71')
Human_R47H_9[["Sex"]] = c('M')
rm(Human_R47H_9.counts)

Human_R47H_10.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_10")
Human_R47H_10 <- CreateSeuratObject(counts = Human_R47H_10.counts, project = "WT_E3E3_M_2", min.cells = 3, min.features = 200)
Human_R47H_10[['percent.mt']] <- PercentageFeatureSet(Human_R47H_10, pattern = "^MT-")
Human_R47H_10[["Condition"]] = c('WT_E3E3_M')
Human_R47H_10[["Condition_1"]] = c('WT_E3E3')
Human_R47H_10[["Condition_2"]] = c('WT_M')
Human_R47H_10[["Condition_3"]] = c('E3E3_M')
Human_R47H_10[["TREM2"]] = c('WT')
Human_R47H_10[["Dx"]] = c('AD/DLBD')
Human_R47H_10[["LBD"]] = c('DLBD')
Human_R47H_10[["Braak"]] = c('5.5')
Human_R47H_10[["Thal"]] = c('5')
Human_R47H_10[["TDP.43"]] = c('0')
Human_R47H_10[["ClinicalDx"]] = c('DLB')
Human_R47H_10[["APOE"]] = c('E3E3')
Human_R47H_10[["Age"]] = c('72')
Human_R47H_10[["Sex"]] = c('M')
rm(Human_R47H_10.counts)

Human_R47H_11.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_11")
Human_R47H_11 <- CreateSeuratObject(counts = Human_R47H_11.counts, project = "R47H_E3E4_F_1", min.cells = 3, min.features = 200)
Human_R47H_11[['percent.mt']] <- PercentageFeatureSet(Human_R47H_11, pattern = "^MT-")
Human_R47H_11[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_11[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_11[["Condition_2"]] = c('R47H_F')
Human_R47H_11[["Condition_3"]] = c('E3E4_F')
Human_R47H_11[["TREM2"]] = c('R47H')
Human_R47H_11[["Dx"]] = c('AD/ALB')
Human_R47H_11[["LBD"]] = c('N/A')
Human_R47H_11[["Braak"]] = c('5.5')
Human_R47H_11[["Thal"]] = c('5')
Human_R47H_11[["TDP.43"]] = c('0')
Human_R47H_11[["ClinicalDx"]] = c('depression')
Human_R47H_11[["APOE"]] = c('E3E4')
Human_R47H_11[["Age"]] = c('85')
Human_R47H_11[["Sex"]] = c('F')
rm(Human_R47H_11.counts)

Human_R47H_13.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_13")
Human_R47H_13 <- CreateSeuratObject(counts = Human_R47H_13.counts, project = "R47H_E3E4_F_2", min.cells = 3, min.features = 200)
Human_R47H_13[['percent.mt']] <- PercentageFeatureSet(Human_R47H_13, pattern = "^MT-")
Human_R47H_13[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_13[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_13[["Condition_2"]] = c('R47H_F')
Human_R47H_13[["Condition_3"]] = c('E3E4_F')
Human_R47H_13[["TREM2"]] = c('R47H')
Human_R47H_13[["Dx"]] = c('AD/VaD')
Human_R47H_13[["LBD"]] = c('N/A')
Human_R47H_13[["Braak"]] = c('5')
Human_R47H_13[["Thal"]] = c('5')
Human_R47H_13[["TDP.43"]] = c('1')
Human_R47H_13[["ClinicalDx"]] = c('dementia')
Human_R47H_13[["APOE"]] = c('E3E4')
Human_R47H_13[["Age"]] = c('87')
Human_R47H_13[["Sex"]] = c('F')
rm(Human_R47H_13.counts)

Human_R47H_14.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_14")
Human_R47H_14 <- CreateSeuratObject(counts = Human_R47H_14.counts, project = "WT_E3E4_F_1", min.cells = 3, min.features = 200)
Human_R47H_14[['percent.mt']] <- PercentageFeatureSet(Human_R47H_14, pattern = "^MT-")
Human_R47H_14[["Condition"]] = c('WT_E3E4_F')
Human_R47H_14[["Condition_1"]] = c('WT_E3E4')
Human_R47H_14[["Condition_2"]] = c('WT_F')
Human_R47H_14[["Condition_3"]] = c('E3E4_F')
Human_R47H_14[["TREM2"]] = c('WT')
Human_R47H_14[["Dx"]] = c('AD/VaD')
Human_R47H_14[["LBD"]] = c('N/A')
Human_R47H_14[["Braak"]] = c('5')
Human_R47H_14[["Thal"]] = c('5')
Human_R47H_14[["TDP.43"]] = c('1')
Human_R47H_14[["ClinicalDx"]] = c('AD')
Human_R47H_14[["APOE"]] = c('E3E4')
Human_R47H_14[["Age"]] = c('86')
Human_R47H_14[["Sex"]] = c('F')
rm(Human_R47H_14.counts)

Human_R47H_15.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_15")
Human_R47H_15 <- CreateSeuratObject(counts = Human_R47H_15.counts, project = "R47H_E3E4_F_3", min.cells = 3, min.features = 200)
Human_R47H_15[['percent.mt']] <- PercentageFeatureSet(Human_R47H_15, pattern = "^MT-")
Human_R47H_15[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_15[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_15[["Condition_2"]] = c('R47H_F')
Human_R47H_15[["Condition_3"]] = c('E3E4_F')
Human_R47H_15[["TREM2"]] = c('R47H')
Human_R47H_15[["Dx"]] = c('AD/VaD')
Human_R47H_15[["LBD"]] = c('N/A')
Human_R47H_15[["Braak"]] = c('6')
Human_R47H_15[["Thal"]] = c('5')
Human_R47H_15[["TDP.43"]] = c('1')
Human_R47H_15[["ClinicalDx"]] = c('AD')
Human_R47H_15[["APOE"]] = c('E3E4')
Human_R47H_15[["Age"]] = c('85')
Human_R47H_15[["Sex"]] = c('F')
rm(Human_R47H_15.counts)

Human_R47H_16.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_16")
Human_R47H_16 <- CreateSeuratObject(counts = Human_R47H_16.counts, project = "WT_E3E4_F_2", min.cells = 3, min.features = 200)
Human_R47H_16[['percent.mt']] <- PercentageFeatureSet(Human_R47H_16, pattern = "^MT-")
Human_R47H_16[["Condition"]] = c('WT_E3E4_F')
Human_R47H_16[["Condition_1"]] = c('WT_E3E4')
Human_R47H_16[["Condition_2"]] = c('WT_F')
Human_R47H_16[["Condition_3"]] = c('E3E4_F')
Human_R47H_16[["TREM2"]] = c('WT')
Human_R47H_16[["Dx"]] = c('AD/VaD')
Human_R47H_16[["LBD"]] = c('N/A')
Human_R47H_16[["Braak"]] = c('4')
Human_R47H_16[["Thal"]] = c('3')
Human_R47H_16[["TDP.43"]] = c('0')
Human_R47H_16[["ClinicalDx"]] = c('AD')
Human_R47H_16[["APOE"]] = c('E3E4')
Human_R47H_16[["Age"]] = c('85')
Human_R47H_16[["Sex"]] = c('F')
rm(Human_R47H_16.counts)

Human_R47H_18.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_18")
Human_R47H_18 <- CreateSeuratObject(counts = Human_R47H_18.counts, project = "WT_E3E4_F_3", min.cells = 3, min.features = 200)
Human_R47H_18[['percent.mt']] <- PercentageFeatureSet(Human_R47H_18, pattern = "^MT-")
Human_R47H_18[["Condition"]] = c('WT_E3E4_F')
Human_R47H_18[["Condition_1"]] = c('WT_E3E4')
Human_R47H_18[["Condition_2"]] = c('WT_F')
Human_R47H_18[["Condition_3"]] = c('E3E4_F')
Human_R47H_18[["TREM2"]] = c('WT')
Human_R47H_18[["Dx"]] = c('AD/VaD')
Human_R47H_18[["LBD"]] = c('N/A')
Human_R47H_18[["Braak"]] = c('6')
Human_R47H_18[["Thal"]] = c('5')
Human_R47H_18[["TDP.43"]] = c('0')
Human_R47H_18[["ClinicalDx"]] = c('AD')
Human_R47H_18[["APOE"]] = c('E3E4')
Human_R47H_18[["Age"]] = c('81')
Human_R47H_18[["Sex"]] = c('F')
rm(Human_R47H_18.counts)

Human_R47H_19.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_19")
Human_R47H_19 <- CreateSeuratObject(counts = Human_R47H_19.counts, project = "R47H_E3E4_M_1", min.cells = 3, min.features = 200)
Human_R47H_19[['percent.mt']] <- PercentageFeatureSet(Human_R47H_19, pattern = "^MT-")
Human_R47H_19[["Condition"]] = c('R47H_E3E4_M')
Human_R47H_19[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_19[["Condition_2"]] = c('R47H_M')
Human_R47H_19[["Condition_3"]] = c('E3E4_M')
Human_R47H_19[["TREM2"]] = c('R47H')
Human_R47H_19[["Dx"]] = c('AD')
Human_R47H_19[["LBD"]] = c('N/A')
Human_R47H_19[["Braak"]] = c('5')
Human_R47H_19[["Thal"]] = c('3')
Human_R47H_19[["TDP.43"]] = c('0')
Human_R47H_19[["ClinicalDx"]] = c('AD')
Human_R47H_19[["APOE"]] = c('E3E4')
Human_R47H_19[["Age"]] = c('81')
Human_R47H_19[["Sex"]] = c('M')
rm(Human_R47H_19.counts)

Human_R47H_20.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_20")
Human_R47H_20 <- CreateSeuratObject(counts = Human_R47H_20.counts, project = "WT_E3E4_M_1", min.cells = 3, min.features = 200)
Human_R47H_20[['percent.mt']] <- PercentageFeatureSet(Human_R47H_20, pattern = "^MT-")
Human_R47H_20[["Condition"]] = c('WT_E3E4_M')
Human_R47H_20[["Condition_1"]] = c('WT_E3E4')
Human_R47H_20[["Condition_2"]] = c('WT_M')
Human_R47H_20[["Condition_3"]] = c('E3E4_M')
Human_R47H_20[["TREM2"]] = c('WT')
Human_R47H_20[["Dx"]] = c('AD')
Human_R47H_20[["LBD"]] = c('N/A')
Human_R47H_20[["Braak"]] = c('6')
Human_R47H_20[["Thal"]] = c('5')
Human_R47H_20[["TDP.43"]] = c('0')
Human_R47H_20[["ClinicalDx"]] = c('CBD')
Human_R47H_20[["APOE"]] = c('E3E4')
Human_R47H_20[["Age"]] = c('80')
Human_R47H_20[["Sex"]] = c('M')
rm(Human_R47H_20.counts)

Human_R47H_21.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_21")
Human_R47H_21 <- CreateSeuratObject(counts = Human_R47H_21.counts, project = "R47H_E3E4_F_4", min.cells = 3, min.features = 200)
Human_R47H_21[['percent.mt']] <- PercentageFeatureSet(Human_R47H_21, pattern = "^MT-")
Human_R47H_21[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_21[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_21[["Condition_2"]] = c('R47H_F')
Human_R47H_21[["Condition_3"]] = c('E3E4_F')
Human_R47H_21[["TREM2"]] = c('R47H')
Human_R47H_21[["Dx"]] = c('AD')
Human_R47H_21[["LBD"]] = c('N/A')
Human_R47H_21[["Braak"]] = c('5')
Human_R47H_21[["Thal"]] = c('5')
Human_R47H_21[["TDP.43"]] = c('0')
Human_R47H_21[["ClinicalDx"]] = c('FTD')
Human_R47H_21[["APOE"]] = c('E3E4')
Human_R47H_21[["Age"]] = c('66')
Human_R47H_21[["Sex"]] = c('F')
rm(Human_R47H_21.counts)

Human_R47H_22.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_22")
Human_R47H_22 <- CreateSeuratObject(counts = Human_R47H_22.counts, project = "WT_E3E4_F_4", min.cells = 3, min.features = 200)
Human_R47H_22[['percent.mt']] <- PercentageFeatureSet(Human_R47H_22, pattern = "^MT-")
Human_R47H_22[["Condition"]] = c('WT_E3E4_F')
Human_R47H_22[["Condition_1"]] = c('WT_E3E4')
Human_R47H_22[["Condition_2"]] = c('WT_F')
Human_R47H_22[["Condition_3"]] = c('E3E4_F')
Human_R47H_22[["TREM2"]] = c('WT')
Human_R47H_22[["Dx"]] = c('AD')
Human_R47H_22[["LBD"]] = c('N/A')
Human_R47H_22[["Braak"]] = c('6')
Human_R47H_22[["Thal"]] = c('5')
Human_R47H_22[["TDP.43"]] = c('1')
Human_R47H_22[["ClinicalDx"]] = c('FTD')
Human_R47H_22[["APOE"]] = c('E3E4')
Human_R47H_22[["Age"]] = c('67')
Human_R47H_22[["Sex"]] = c('F')
rm(Human_R47H_22.counts)

Human_R47H_23.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_23")
Human_R47H_23 <- CreateSeuratObject(counts = Human_R47H_23.counts, project = "R47H_E3E4_F_5", min.cells = 3, min.features = 200)
Human_R47H_23[['percent.mt']] <- PercentageFeatureSet(Human_R47H_23, pattern = "^MT-")
Human_R47H_23[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_23[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_23[["Condition_2"]] = c('R47H_F')
Human_R47H_23[["Condition_3"]] = c('E3E4_F')
Human_R47H_23[["TREM2"]] = c('R47H')
Human_R47H_23[["Dx"]] = c('AD')
Human_R47H_23[["LBD"]] = c('N/A')
Human_R47H_23[["Braak"]] = c('6')
Human_R47H_23[["Thal"]] = c('5')
Human_R47H_23[["TDP.43"]] = c('1')
Human_R47H_23[["ClinicalDx"]] = c('AD')
Human_R47H_23[["APOE"]] = c('E3E4')
Human_R47H_23[["Age"]] = c('78')
Human_R47H_23[["Sex"]] = c('F')
rm(Human_R47H_23.counts)

Human_R47H_24.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_24")
Human_R47H_24 <- CreateSeuratObject(counts = Human_R47H_24.counts, project = "WT_E3E4_F_5", min.cells = 3, min.features = 200)
Human_R47H_24[['percent.mt']] <- PercentageFeatureSet(Human_R47H_24, pattern = "^MT-")
Human_R47H_24[["Condition"]] = c('WT_E3E4_F')
Human_R47H_24[["Condition_1"]] = c('WT_E3E4')
Human_R47H_24[["Condition_2"]] = c('WT_F')
Human_R47H_24[["Condition_3"]] = c('E3E4_F')
Human_R47H_24[["TREM2"]] = c('WT')
Human_R47H_24[["Dx"]] = c('AD')
Human_R47H_24[["LBD"]] = c('N/A')
Human_R47H_24[["Braak"]] = c('6')
Human_R47H_24[["Thal"]] = c('5')
Human_R47H_24[["TDP.43"]] = c('0')
Human_R47H_24[["ClinicalDx"]] = c('AD')
Human_R47H_24[["APOE"]] = c('E3E4')
Human_R47H_24[["Age"]] = c('81')
Human_R47H_24[["Sex"]] = c('F')
rm(Human_R47H_24.counts)

Human_R47H_25.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_25")
Human_R47H_25 <- CreateSeuratObject(counts = Human_R47H_25.counts, project = "R47H_E3E4_M_2", min.cells = 3, min.features = 200)
Human_R47H_25[['percent.mt']] <- PercentageFeatureSet(Human_R47H_25, pattern = "^MT-")
Human_R47H_25[["Condition"]] = c('R47H_E3E4_M')
Human_R47H_25[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_25[["Condition_2"]] = c('R47H_M')
Human_R47H_25[["Condition_3"]] = c('E3E4_M')
Human_R47H_25[["TREM2"]] = c('R47H')
Human_R47H_25[["Dx"]] = c('AD')
Human_R47H_25[["LBD"]] = c('N/A')
Human_R47H_25[["Braak"]] = c('4.5')
Human_R47H_25[["Thal"]] = c('5')
Human_R47H_25[["TDP.43"]] = c('0')
Human_R47H_25[["ClinicalDx"]] = c('AD')
Human_R47H_25[["APOE"]] = c('E3E4')
Human_R47H_25[["Age"]] = c('88')
Human_R47H_25[["Sex"]] = c('M')
rm(Human_R47H_25.counts)

Human_R47H_26.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_26")
Human_R47H_26 <- CreateSeuratObject(counts = Human_R47H_26.counts, project = "WT_E3E4_M_2", min.cells = 3, min.features = 200)
Human_R47H_26[['percent.mt']] <- PercentageFeatureSet(Human_R47H_26, pattern = "^MT-")
Human_R47H_26[["Condition"]] = c('WT_E3E4_M')
Human_R47H_26[["Condition_1"]] = c('WT_E3E4')
Human_R47H_26[["Condition_2"]] = c('WT_M')
Human_R47H_26[["Condition_3"]] = c('E3E4_M')
Human_R47H_26[["TREM2"]] = c('WT')
Human_R47H_26[["Dx"]] = c('AD')
Human_R47H_26[["LBD"]] = c('N/A')
Human_R47H_26[["Braak"]] = c('5')
Human_R47H_26[["Thal"]] = c('5')
Human_R47H_26[["TDP.43"]] = c('1')
Human_R47H_26[["ClinicalDx"]] = c('VaD')
Human_R47H_26[["APOE"]] = c('E3E4')
Human_R47H_26[["Age"]] = c('89')
Human_R47H_26[["Sex"]] = c('M')
rm(Human_R47H_26.counts)

Human_R47H_27.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_27")
Human_R47H_27 <- CreateSeuratObject(counts = Human_R47H_27.counts, project = "R47H_E3E4_F_6", min.cells = 3, min.features = 200)
Human_R47H_27[['percent.mt']] <- PercentageFeatureSet(Human_R47H_27, pattern = "^MT-")
Human_R47H_27[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_27[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_27[["Condition_2"]] = c('R47H_F')
Human_R47H_27[["Condition_3"]] = c('E3E4_F')
Human_R47H_27[["TREM2"]] = c('R47H')
Human_R47H_27[["Dx"]] = c('AD')
Human_R47H_27[["LBD"]] = c('N/A')
Human_R47H_27[["Braak"]] = c('5')
Human_R47H_27[["Thal"]] = c('3')
Human_R47H_27[["TDP.43"]] = c('1')
Human_R47H_27[["ClinicalDx"]] = c('AD')
Human_R47H_27[["APOE"]] = c('E3E4')
Human_R47H_27[["Age"]] = c('74')
Human_R47H_27[["Sex"]] = c('F')
rm(Human_R47H_27.counts)

Human_R47H_29.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_29")
Human_R47H_29 <- CreateSeuratObject(counts = Human_R47H_29.counts, project = "R47H_E3E4_M_3", min.cells = 3, min.features = 200)
Human_R47H_29[['percent.mt']] <- PercentageFeatureSet(Human_R47H_29, pattern = "^MT-")
Human_R47H_29[["Condition"]] = c('R47H_E3E4_M')
Human_R47H_29[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_29[["Condition_2"]] = c('R47H_M')
Human_R47H_29[["Condition_3"]] = c('E3E4_M')
Human_R47H_29[["TREM2"]] = c('R47H')
Human_R47H_29[["Dx"]] = c('AD/TLBD')
Human_R47H_29[["LBD"]] = c('TLBD')
Human_R47H_29[["Braak"]] = c('6')
Human_R47H_29[["Thal"]] = c('5')
Human_R47H_29[["TDP.43"]] = c('0')
Human_R47H_29[["ClinicalDx"]] = c('FTD/SD')
Human_R47H_29[["APOE"]] = c('E3E4')
Human_R47H_29[["Age"]] = c('76')
Human_R47H_29[["Sex"]] = c('M')
rm(Human_R47H_29.counts)

Human_R47H_30.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_30")
Human_R47H_30 <- CreateSeuratObject(counts = Human_R47H_30.counts, project = "WT_E3E4_M_3", min.cells = 3, min.features = 200)
Human_R47H_30[['percent.mt']] <- PercentageFeatureSet(Human_R47H_30, pattern = "^MT-")
Human_R47H_30[["Condition"]] = c('WT_E3E4_M')
Human_R47H_30[["Condition_1"]] = c('WT_E3E4')
Human_R47H_30[["Condition_2"]] = c('WT_M')
Human_R47H_30[["Condition_3"]] = c('E3E4_M')
Human_R47H_30[["TREM2"]] = c('WT')
Human_R47H_30[["Dx"]] = c('AD/TLBD')
Human_R47H_30[["LBD"]] = c('TLBD')
Human_R47H_30[["Braak"]] = c('6')
Human_R47H_30[["Thal"]] = c('5')
Human_R47H_30[["TDP.43"]] = c('0')
Human_R47H_30[["ClinicalDx"]] = c('AD')
Human_R47H_30[["APOE"]] = c('E3E4')
Human_R47H_30[["Age"]] = c('79')
Human_R47H_30[["Sex"]] = c('M')
rm(Human_R47H_30.counts)

Human_R47H_31.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_31")
Human_R47H_31 <- CreateSeuratObject(counts = Human_R47H_31.counts, project = "R47H_E3E4_F_7", min.cells = 3, min.features = 200)
Human_R47H_31[['percent.mt']] <- PercentageFeatureSet(Human_R47H_31, pattern = "^MT-")
Human_R47H_31[["Condition"]] = c('R47H_E3E4_F')
Human_R47H_31[["Condition_1"]] = c('R47H_E3E4')
Human_R47H_31[["Condition_2"]] = c('R47H_F')
Human_R47H_31[["Condition_3"]] = c('E3E4_F')
Human_R47H_31[["TREM2"]] = c('R47H')
Human_R47H_31[["Dx"]] = c('AD/TLBD')
Human_R47H_31[["LBD"]] = c('TLBD')
Human_R47H_31[["Braak"]] = c('6')
Human_R47H_31[["Thal"]] = c('5')
Human_R47H_31[["TDP.43"]] = c('1')
Human_R47H_31[["ClinicalDx"]] = c('AD')
Human_R47H_31[["APOE"]] = c('E3E4')
Human_R47H_31[["Age"]] = c('77')
Human_R47H_31[["Sex"]] = c('F')
rm(Human_R47H_31.counts)


Human_R47H_33.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_33")
Human_R47H_33 <- CreateSeuratObject(counts = Human_R47H_33.counts, project = "R47H_E4E4_F_1", min.cells = 3, min.features = 200)
Human_R47H_33[['percent.mt']] <- PercentageFeatureSet(Human_R47H_33, pattern = "^MT-")
Human_R47H_33[["Condition"]] = c('R47H_E4E4_F')
Human_R47H_33[["Condition_1"]] = c('R47H_E4E4')
Human_R47H_33[["Condition_2"]] = c('R47H_F')
Human_R47H_33[["Condition_3"]] = c('E4E4_F')
Human_R47H_33[["TREM2"]] = c('R47H')
Human_R47H_33[["Dx"]] = c('AD/DLBD')
Human_R47H_33[["LBD"]] = c('DLBD')
Human_R47H_33[["Braak"]] = c('6')
Human_R47H_33[["Thal"]] = c('5')
Human_R47H_33[["TDP.43"]] = c('0')
Human_R47H_33[["ClinicalDx"]] = c('AD/VaD/PD')
Human_R47H_33[["APOE"]] = c('E4E4')
Human_R47H_33[["Age"]] = c('82')
Human_R47H_33[["Sex"]] = c('F')
rm(Human_R47H_33.counts)


Human_R47H_35.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_35")
Human_R47H_35 <- CreateSeuratObject(counts = Human_R47H_35.counts, project = "R47H_E4E4_M_1", min.cells = 3, min.features = 200)
Human_R47H_35[['percent.mt']] <- PercentageFeatureSet(Human_R47H_35, pattern = "^MT-")
Human_R47H_35[["Condition"]] = c('R47H_E4E4_M')
Human_R47H_35[["Condition_1"]] = c('R47H_E4E4')
Human_R47H_35[["Condition_2"]] = c('R47H_M')
Human_R47H_35[["Condition_3"]] = c('E4E4_M')
Human_R47H_35[["TREM2"]] = c('R47H')
Human_R47H_35[["Dx"]] = c('AD')
Human_R47H_35[["LBD"]] = c('N/A')
Human_R47H_35[["Braak"]] = c('5')
Human_R47H_35[["Thal"]] = c('5')
Human_R47H_35[["TDP.43"]] = c('0')
Human_R47H_35[["ClinicalDx"]] = c('AD_v_DLB')
Human_R47H_35[["APOE"]] = c('E4E4')
Human_R47H_35[["Age"]] = c('71')
Human_R47H_35[["Sex"]] = c('M')
rm(Human_R47H_35.counts)

Human_R47H_37.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_37")
Human_R47H_37 <- CreateSeuratObject(counts = Human_R47H_37.counts, project = "R47H_E4E4_M_2", min.cells = 3, min.features = 200)
Human_R47H_37[['percent.mt']] <- PercentageFeatureSet(Human_R47H_37, pattern = "^MT-")
Human_R47H_37[["Condition"]] = c('R47H_E4E4_M')
Human_R47H_37[["Condition_1"]] = c('R47H_E4E4')
Human_R47H_37[["Condition_2"]] = c('R47H_M')
Human_R47H_37[["Condition_3"]] = c('E4E4_M')
Human_R47H_37[["TREM2"]] = c('R47H')
Human_R47H_37[["Dx"]] = c('AD/CAA')
Human_R47H_37[["LBD"]] = c('N/A')
Human_R47H_37[["Braak"]] = c('6')
Human_R47H_37[["Thal"]] = c('4')
Human_R47H_37[["TDP.43"]] = c('0')
Human_R47H_37[["ClinicalDx"]] = c('AD')
Human_R47H_37[["APOE"]] = c('E4E4')
Human_R47H_37[["Age"]] = c('79')
Human_R47H_37[["Sex"]] = c('M')
rm(Human_R47H_37.counts)

Human_R47H_38.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_38")
Human_R47H_38 <- CreateSeuratObject(counts = Human_R47H_38.counts, project = "WT_E4E4_M", min.cells = 3, min.features = 200)
Human_R47H_38[['percent.mt']] <- PercentageFeatureSet(Human_R47H_38, pattern = "^MT-")
Human_R47H_38[["Condition"]] = c('WT_E4E4_M')
Human_R47H_38[["Condition_1"]] = c('WT_E4E4')
Human_R47H_38[["Condition_2"]] = c('WT_M')
Human_R47H_38[["Condition_3"]] = c('E4E4_M')
Human_R47H_38[["TREM2"]] = c('WT')
Human_R47H_38[["Dx"]] = c('AD')
Human_R47H_38[["LBD"]] = c('N/A')
Human_R47H_38[["Braak"]] = c('5.5')
Human_R47H_38[["Thal"]] = c('5')
Human_R47H_38[["TDP.43"]] = c('N/A')
Human_R47H_38[["ClinicalDx"]] = c('AD')
Human_R47H_38[["APOE"]] = c('E4E4')
Human_R47H_38[["Age"]] = c('80')
Human_R47H_38[["Sex"]] = c('M')
rm(Human_R47H_38.counts)

Human_R47H_39.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_39")
Human_R47H_39 <- CreateSeuratObject(counts = Human_R47H_39.counts, project = "R47H_E4E4_F_2", min.cells = 3, min.features = 200)
Human_R47H_39[['percent.mt']] <- PercentageFeatureSet(Human_R47H_39, pattern = "^MT-")
Human_R47H_39[["Condition"]] = c('R47H_E4E4_F')
Human_R47H_39[["Condition_1"]] = c('R47H_E4E4')
Human_R47H_39[["Condition_2"]] = c('R47H_F')
Human_R47H_39[["Condition_3"]] = c('E4E4_F')
Human_R47H_39[["TREM2"]] = c('R47H')
Human_R47H_39[["Dx"]] = c('AD/CAA/BLBD')
Human_R47H_39[["LBD"]] = c('BLBD')
Human_R47H_39[["Braak"]] = c('6')
Human_R47H_39[["Thal"]] = c('5')
Human_R47H_39[["TDP.43"]] = c('1')
Human_R47H_39[["ClinicalDx"]] = c('AD')
Human_R47H_39[["APOE"]] = c('E4E4')
Human_R47H_39[["Age"]] = c('82')
Human_R47H_39[["Sex"]] = c('F')
rm(Human_R47H_39.counts)

Human_R47H_40.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_R47H/Raw_matrix/Human_R47H_40")
Human_R47H_40 <- CreateSeuratObject(counts = Human_R47H_40.counts, project = "WT_E4E4_F", min.cells = 3, min.features = 200)
Human_R47H_40[['percent.mt']] <- PercentageFeatureSet(Human_R47H_40, pattern = "^MT-")
Human_R47H_40[["Condition"]] = c('WT_E4E4_F')
Human_R47H_40[["Condition_1"]] = c('WT_E4E4')
Human_R47H_40[["Condition_2"]] = c('WT_F')
Human_R47H_40[["Condition_3"]] = c('E4E4_F')
Human_R47H_40[["TREM2"]] = c('WT')
Human_R47H_40[["Dx"]] = c('AD/CAA/BLBD')
Human_R47H_40[["LBD"]] = c('BLBD')
Human_R47H_40[["Braak"]] = c('6')
Human_R47H_40[["Thal"]] = c('5')
Human_R47H_40[["TDP.43"]] = c('N/A')
Human_R47H_40[["ClinicalDx"]] = c('AD')
Human_R47H_40[["APOE"]] = c('E4E4')
Human_R47H_40[["Age"]] = c('74')
Human_R47H_40[["Sex"]] = c('F')
rm(Human_R47H_40.counts)

setwd("/athena/ganlab/scratch/lif4001/Human_R47H/data_analysis/DF_2ndRound")
all <- Human_R47H_1
pdf("Human_R47H_1_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_1_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_1_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5882) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_229", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_1_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_1_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_229" #visualizing the singlet vs doublet cells
pdf("Human_R47H_1_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_1_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_1_singlets.rds")

singlets<-readRDS("Human_R47H_1_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_1_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_1_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_1_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_1_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_1_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_1_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_1_singlets_PCA.rds")

all <- Human_R47H_2
pdf("Human_R47H_2_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_2_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_2_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8112) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.27, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.27, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.27_495", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_2_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_2_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.27_495" #visualizing the singlet vs doublet cells
pdf("Human_R47H_2_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_2_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_2_singlets.rds")

singlets<-readRDS("Human_R47H_2_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_2_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_2_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_2_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_2_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_2_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_2_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_2_singlets_PCA.rds")


all <- Human_R47H_3
pdf("Human_R47H_3_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_3_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_3_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8200) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_500", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_3_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_3_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_500" #visualizing the singlet vs doublet cells
pdf("Human_R47H_3_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_3_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_3_singlets.rds")

singlets<-readRDS("Human_R47H_3_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_3_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_3_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_3_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_3_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_3_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_3_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_3_singlets_PCA.rds")


all <- Human_R47H_4
pdf("Human_R47H_4_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_4_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_4_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_4_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8002) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_488", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_4_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_4_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_488" #visualizing the singlet vs doublet cells
pdf("Human_R47H_4_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_4_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_4_singlets.rds")

singlets<-readRDS("Human_R47H_4_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_4_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_4_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_4_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_4_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_4_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_4_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_4_singlets_PCA.rds")

all <- Human_R47H_5
pdf("Human_R47H_5_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_5_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_5_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_5_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_5_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_5_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_5_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5579) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_218", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_5_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_5_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_218" #visualizing the singlet vs doublet cells
pdf("Human_R47H_5_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_5_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_5_singlets.rds")

singlets<-readRDS("Human_R47H_5_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_5_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_5_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_5_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_5_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_5_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_5_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_5_singlets_PCA.rds")

all <- Human_R47H_6
pdf("Human_R47H_6_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_6_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_6_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_6_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_6_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_6_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_6_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8652) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_528", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_6_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_6_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_528" #visualizing the singlet vs doublet cells
pdf("Human_R47H_6_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_6_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_6_singlets.rds")

singlets<-readRDS("Human_R47H_6_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_6_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_6_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_6_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_6_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_6_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_6_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_6_singlets_PCA.rds")


all <- Human_R47H_7
pdf("Human_R47H_7_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_7_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_7_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_7_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_7_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_7_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_7_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9367) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_646", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_7_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_7_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.02_646" #visualizing the singlet vs doublet cells
pdf("Human_R47H_7_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_7_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_7_singlets.rds")

singlets<-readRDS("Human_R47H_7_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_7_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_7_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_7_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_7_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_7_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_7_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_7_singlets_PCA.rds")

all <- Human_R47H_8
pdf("Human_R47H_8_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_8_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_8_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_8_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_8_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_8_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_8_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.016*2820) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_45", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_8_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_8_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.01_45" #visualizing the singlet vs doublet cells
pdf("Human_R47H_8_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_8_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_8_singlets.rds")

singlets<-readRDS("Human_R47H_8_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_8_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_8_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_8_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_8_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_8_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_8_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_8_singlets_PCA.rds")


all <- Human_R47H_9
pdf("Human_R47H_9_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_9_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_9_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_9_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_9_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_9_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_9_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7482) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_404", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_9_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_9_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_404" #visualizing the singlet vs doublet cells
pdf("Human_R47H_9_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_9_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_9_singlets.rds")

singlets<-readRDS("Human_R47H_9_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_9_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_9_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_9_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_9_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_9_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_9_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_9_singlets_PCA.rds")


all <- Human_R47H_10
pdf("Human_R47H_10_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_10_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_10_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_10_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_10_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_10_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_10_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6015) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_277", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_10_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_10_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_277" #visualizing the singlet vs doublet cells
pdf("Human_R47H_10_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_10_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_10_singlets.rds")

singlets<-readRDS("Human_R47H_10_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_10_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_10_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_10_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_10_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_10_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_10_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_10_singlets_PCA.rds")


all <- Human_R47H_11
pdf("Human_R47H_11_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_11_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_11_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_11_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_11_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_11_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_11_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6730) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_310", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_11_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_11_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_310" #visualizing the singlet vs doublet cells
pdf("Human_R47H_11_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_11_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_11_singlets.rds")

singlets<-readRDS("Human_R47H_11_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_11_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_11_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_11_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_11_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_11_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_11_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_11_singlets_PCA.rds")


all <- Human_R47H_13
pdf("Human_R47H_13_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_13_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_13_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_13_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_13_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_13_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_13_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5842) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_228", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_13_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_13_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.03_228" #visualizing the singlet vs doublet cells
pdf("Human_R47H_13_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_13_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_13_singlets.rds")

singlets<-readRDS("Human_R47H_13_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_13_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_13_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_13_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_13_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_13_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_13_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_13_singlets_PCA.rds")

all <- Human_R47H_14
pdf("Human_R47H_14_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_14_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_14_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_14_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_14_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_14_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_14_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7130) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_385", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_14_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_14_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_385" #visualizing the singlet vs doublet cells
pdf("Human_R47H_14_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_14_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_14_singlets.rds")

singlets<-readRDS("Human_R47H_14_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_14_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_14_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_14_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_14_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_14_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_14_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_14_singlets_PCA.rds")

all <- Human_R47H_15
pdf("Human_R47H_15_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_15_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_15_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_15_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_15_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_15_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_15_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11081) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_931", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_15_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_15_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_931" #visualizing the singlet vs doublet cells
pdf("Human_R47H_15_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_15_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_15_singlets.rds")

singlets<-readRDS("Human_R47H_15_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_15_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_15_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_15_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_15_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_15_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_15_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_15_singlets_PCA.rds")

all <- Human_R47H_16
pdf("Human_R47H_16_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_16_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_16_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_16_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_16_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_16_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_16_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8320) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_508", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_16_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_16_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_508" #visualizing the singlet vs doublet cells
pdf("Human_R47H_16_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_16_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_16_singlets.rds")

singlets<-readRDS("Human_R47H_16_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_16_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_16_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_16_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_16_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_16_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_16_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_16_singlets_PCA.rds")

all <- Human_R47H_18
pdf("Human_R47H_18_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_18_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_18_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_18_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_18_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_18_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_18_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8305) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_507", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_18_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_18_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_507" #visualizing the singlet vs doublet cells
pdf("Human_R47H_18_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_18_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_18_singlets.rds")

singlets<-readRDS("Human_R47H_18_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_18_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_18_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_18_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_18_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_18_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_18_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_18_singlets_PCA.rds")

all <- Human_R47H_19
pdf("Human_R47H_19_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_19_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_19_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_19_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_19_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_19_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_19_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*4339) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_135", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_19_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_19_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.01_135" #visualizing the singlet vs doublet cells
pdf("Human_R47H_19_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_19_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_19_singlets.rds")

singlets<-readRDS("Human_R47H_19_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_19_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_19_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_19_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_19_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_19_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_19_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_19_singlets_PCA.rds")

all <- Human_R47H_20
pdf("Human_R47H_20_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_20_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_20_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_20_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_20_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_20_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_20_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5026) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_196", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_20_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_20_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.01_196" #visualizing the singlet vs doublet cells
pdf("Human_R47H_20_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_20_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_20_singlets.rds")

singlets<-readRDS("Human_R47H_20_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_20_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_20_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_20_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_20_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_20_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_20_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_20_singlets_PCA.rds")

all <- Human_R47H_21
pdf("Human_R47H_21_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_21_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_21_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_21_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_21_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_21_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_21_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5610) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_219", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_21_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_21_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_219" #visualizing the singlet vs doublet cells
pdf("Human_R47H_21_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_21_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_21_singlets.rds")

singlets<-readRDS("Human_R47H_21_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_21_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_21_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_21_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_21_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_21_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_21_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_21_singlets_PCA.rds")

all <- Human_R47H_22
pdf("Human_R47H_22_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_22_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_22_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_22_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_22_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_22_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_22_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7092) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_383", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_22_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_22_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_383" #visualizing the singlet vs doublet cells
pdf("Human_R47H_22_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_22_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_22_singlets.rds")

singlets<-readRDS("Human_R47H_22_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_22_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_22_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_22_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_22_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_22_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_22_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_22_singlets_PCA.rds")

all <- Human_R47H_23
pdf("Human_R47H_23_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_23_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_23_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_23_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_23_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_23_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_23_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5632) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_220", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_23_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_23_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_220" #visualizing the singlet vs doublet cells
pdf("Human_R47H_23_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_23_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_23_singlets.rds")

singlets<-readRDS("Human_R47H_23_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_23_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_23_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_23_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_23_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_23_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_23_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_23_singlets_PCA.rds")

all <- Human_R47H_24
pdf("Human_R47H_24_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_24_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_24_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_24_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_24_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_24_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_24_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5823) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_227", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_24_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_24_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_227" #visualizing the singlet vs doublet cells
pdf("Human_R47H_24_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_24_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_24_singlets.rds")

singlets<-readRDS("Human_R47H_24_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_24_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_24_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_24_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_24_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_24_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_24_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_24_singlets_PCA.rds")

all <- Human_R47H_25
pdf("Human_R47H_25_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_25_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_25_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_25_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_25_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_25_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_25_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7227) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_390", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_25_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_25_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_390" #visualizing the singlet vs doublet cells
pdf("Human_R47H_25_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_25_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_25_singlets.rds")

singlets<-readRDS("Human_R47H_25_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_25_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_25_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_25_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_25_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_25_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_25_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_25_singlets_PCA.rds")

all <- Human_R47H_26
pdf("Human_R47H_26_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_26_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_26_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_26_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_26_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_26_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_26_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8820) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_538", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_26_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_26_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_538" #visualizing the singlet vs doublet cells
pdf("Human_R47H_26_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_26_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_26_singlets.rds")

singlets<-readRDS("Human_R47H_26_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_26_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_26_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_26_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_26_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_26_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_26_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_26_singlets_PCA.rds")

all <- Human_R47H_27
pdf("Human_R47H_27_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_27_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_27_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_27_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_27_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_27_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_27_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5516) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_215", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_27_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_27_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_215" #visualizing the singlet vs doublet cells
pdf("Human_R47H_27_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_27_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_27_singlets.rds")

singlets<-readRDS("Human_R47H_27_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_27_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_27_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_27_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_27_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_27_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_27_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_27_singlets_PCA.rds")

all <- Human_R47H_29
pdf("Human_R47H_29_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_29_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_29_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_29_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_29_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_29_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_29_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6397) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_294", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_29_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_29_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_294" #visualizing the singlet vs doublet cells
pdf("Human_R47H_29_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_29_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_29_singlets.rds")

singlets<-readRDS("Human_R47H_29_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_29_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_29_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_29_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_29_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_29_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_29_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_29_singlets_PCA.rds")

all <- Human_R47H_30
pdf("Human_R47H_30_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_30_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_30_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_30_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_30_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_30_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_30_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7948) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_429", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_30_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_30_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_429" #visualizing the singlet vs doublet cells
pdf("Human_R47H_30_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_30_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_30_singlets.rds")

singlets<-readRDS("Human_R47H_30_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_30_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_30_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_30_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_30_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_30_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_30_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_30_singlets_PCA.rds")

all <- Human_R47H_31
pdf("Human_R47H_31_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_31_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_31_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_31_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_31_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_31_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_31_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7666) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_414", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_31_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_31_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_414" #visualizing the singlet vs doublet cells
pdf("Human_R47H_31_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_31_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_31_singlets.rds")

singlets<-readRDS("Human_R47H_31_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_31_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_31_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_31_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_31_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_31_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_31_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_31_singlets_PCA.rds")

all <- Human_R47H_33
pdf("Human_R47H_33_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_33_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_33_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_33_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_33_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_33_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_33_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7914) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_427", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_33_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_33_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_427" #visualizing the singlet vs doublet cells
pdf("Human_R47H_33_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_33_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_33_singlets.rds")

singlets<-readRDS("Human_R47H_33_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_33_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_33_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_33_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_33_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_33_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_33_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_33_singlets_PCA.rds")


all <- Human_R47H_35
pdf("Human_R47H_35_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_35_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_35_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_35_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_35_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_35_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_35_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9640) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_665", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_35_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_35_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_665" #visualizing the singlet vs doublet cells
pdf("Human_R47H_35_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_35_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_35_singlets.rds")

singlets<-readRDS("Human_R47H_35_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_35_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_35_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_35_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_35_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_35_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_35_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_35_singlets_PCA.rds")

all <- Human_R47H_37
pdf("Human_R47H_37_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_37_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_37_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_37_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_37_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_37_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_37_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6560) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_302", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_37_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_37_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_302" #visualizing the singlet vs doublet cells
pdf("Human_R47H_37_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_37_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_37_singlets.rds")

singlets<-readRDS("Human_R47H_37_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_37_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_37_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_37_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_37_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_37_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_37_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_37_singlets_PCA.rds")

all <- Human_R47H_38
pdf("Human_R47H_38_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_38_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_38_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_38_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_38_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_38_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_38_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8445) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.28, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.28_515", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_38_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_38_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.28_515" #visualizing the singlet vs doublet cells
pdf("Human_R47H_38_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_38_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_38_singlets.rds")

singlets<-readRDS("Human_R47H_38_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_38_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_38_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_38_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_38_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_38_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_38_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_38_singlets_PCA.rds")

all <- Human_R47H_39
pdf("Human_R47H_39_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_39_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_39_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_39_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_39_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_39_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_39_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7919) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.19, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.19, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.19_428", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_39_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_39_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.19_428" #visualizing the singlet vs doublet cells
pdf("Human_R47H_39_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_39_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_39_singlets.rds")

singlets<-readRDS("Human_R47H_39_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_39_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_39_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_39_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_39_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_39_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_39_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_39_singlets_PCA.rds")

all <- Human_R47H_40
pdf("Human_R47H_40_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Human_R47H_40_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
 
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Human_R47H_40_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Human_R47H_40_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"Human_R47H_40_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("Human_R47H_40_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Human_R47H_40_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7362) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_398", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Human_R47H_40_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Human_R47H_40_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_398" #visualizing the singlet vs doublet cells
pdf("Human_R47H_40_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

saveRDS(all,"Human_R47H_40_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"Human_R47H_40_singlets.rds")

singlets<-readRDS("Human_R47H_40_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("Human_R47H_40_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("Human_R47H_40_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_40_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("Human_R47H_40_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("Human_R47H_40_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Human_R47H_40_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"Human_R47H_40_singlets_PCA.rds")



