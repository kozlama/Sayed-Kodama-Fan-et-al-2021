###############################################################################################
# Pre-processing for data used to generate Figure 1 from Sayed, Kodama, Fan, et al. 2021 
# Total of 55 human AD R47H vs CV samples and non-AD samples
# This script is: STEP 1 of 5 (PART 2 of 2) - 21 human AD R47H vs CV samples from UPenn

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(ggplot2)
library(DoubletFinder)

#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Human_UPenn/data_analysis/DF_2ndRound")

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
Gan_43.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_43/outs/filtered_feature_bc_matrix")
Gan_43 <- CreateSeuratObject(counts = Gan_43.counts, project = "AD_WT_E3E3_F_1", min.cells = 3, min.features = 200)
Gan_43[["Condition"]] = c('AD_WT_E3E3_F')
Gan_43[["Condition_1"]] = c('WT_E3E3')
Gan_43[["Condition_2"]] = c('WT_F')
Gan_43[["Condition_3"]] = c('E3E3_F')
Gan_43[["TREM2"]] = c('WT')
Gan_43[["Dx"]] = c('AD')
Gan_43[["LBD"]] = c('N/A')
Gan_43[["Braak"]] = c('6')
Gan_43[["Thal"]] = c('N/A')
Gan_43[["TDP.43"]] = c('N/A')
Gan_43[["ClinicalDx"]] = c('Dementia with Lewy Bodies')
Gan_43[["APOE"]] = c('E3E3')
Gan_43[["Age"]] = c('83')
Gan_43[["Age_Onset"]] = c('77')
Gan_43[["Sex"]] = c('F')
Gan_43[["PMI"]] = c('4')
Gan_43[["INDDID"]] = c('108953')
rm(Gan_43.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_43[["percent.mt"]] <- PercentageFeatureSet(object = Gan_43, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_44.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_44/outs/filtered_feature_bc_matrix")
Gan_44 <- CreateSeuratObject(counts = Gan_44.counts, project = "AD_WT_E3E3_F_2", min.cells = 3, min.features = 200)
Gan_44[["Condition"]] = c('AD_WT_E3E3_F')
Gan_44[["Condition_1"]] = c('WT_E3E3')
Gan_44[["Condition_2"]] = c('WT_F')
Gan_44[["Condition_3"]] = c('E3E3_F')
Gan_44[["TREM2"]] = c('WT')
Gan_44[["Dx"]] = c('AD')
Gan_44[["LBD"]] = c('N/A')
Gan_44[["Braak"]] = c('6')
Gan_44[["Thal"]] = c('N/A')
Gan_44[["TDP.43"]] = c('N/A')
Gan_44[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_44[["APOE"]] = c('E3E3')
Gan_44[["Age"]] = c('87')
Gan_44[["Age_Onset"]] = c('78')
Gan_44[["Sex"]] = c('F')
Gan_44[["PMI"]] = c('16')
Gan_44[["INDDID"]] = c('105007')
rm(Gan_44.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_44[["percent.mt"]] <- PercentageFeatureSet(object = Gan_44, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_45.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_45/outs/filtered_feature_bc_matrix")
Gan_45 <- CreateSeuratObject(counts = Gan_45.counts, project = "AD_WT_E4E4_F", min.cells = 3, min.features = 200)
Gan_45[["Condition"]] = c('AD_WT_E4E4_F')
Gan_45[["Condition_1"]] = c('WT_E4E4')
Gan_45[["Condition_2"]] = c('WT_F')
Gan_45[["Condition_3"]] = c('E4E4_F')
Gan_45[["TREM2"]] = c('WT')
Gan_45[["Dx"]] = c('AD')
Gan_45[["LBD"]] = c('N/A')
Gan_45[["Braak"]] = c('5')
Gan_45[["Thal"]] = c('N/A')
Gan_45[["TDP.43"]] = c('N/A')
Gan_45[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_45[["APOE"]] = c('E4E4')
Gan_45[["Age"]] = c('72')
Gan_45[["Age_Onset"]] = c('66')
Gan_45[["Sex"]] = c('F')
Gan_45[["PMI"]] = c('5')
Gan_45[["INDDID"]] = c('102475')
rm(Gan_45.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_45[["percent.mt"]] <- PercentageFeatureSet(object = Gan_45, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_46.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_46/outs/filtered_feature_bc_matrix")
Gan_46 <- CreateSeuratObject(counts = Gan_46.counts, project = "AD_WT_E3E3_M", min.cells = 3, min.features = 200)
Gan_46[["Condition"]] = c('AD_WT_E3E3_M')
Gan_46[["Condition_1"]] = c('WT_E3E3')
Gan_46[["Condition_2"]] = c('WT_M')
Gan_46[["Condition_3"]] = c('E3E3_M')
Gan_46[["TREM2"]] = c('WT')
Gan_46[["Dx"]] = c('AD')
Gan_46[["LBD"]] = c('N/A')
Gan_46[["Braak"]] = c('5')
Gan_46[["Thal"]] = c('N/A')
Gan_46[["TDP.43"]] = c('N/A')
Gan_46[["ClinicalDx"]] = c('Dementia of undetermined etiology')
Gan_46[["APOE"]] = c('E3E3')
Gan_46[["Age"]] = c('82')
Gan_46[["Age_Onset"]] = c('74')
Gan_46[["Sex"]] = c('M')
Gan_46[["PMI"]] = c('39')
Gan_46[["INDDID"]] = c('107928')
rm(Gan_46.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_46[["percent.mt"]] <- PercentageFeatureSet(object = Gan_46, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_47.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_47/outs/filtered_feature_bc_matrix")
Gan_47 <- CreateSeuratObject(counts = Gan_47.counts, project = "AD_WT_E3E4_M_1", min.cells = 3, min.features = 200)
Gan_47[["Condition"]] = c('AD_WT_E3E4_M')
Gan_47[["Condition_1"]] = c('WT_E3E4')
Gan_47[["Condition_2"]] = c('WT_M')
Gan_47[["Condition_3"]] = c('E3E4_M')
Gan_47[["TREM2"]] = c('WT')
Gan_47[["Dx"]] = c('AD')
Gan_47[["LBD"]] = c('N/A')
Gan_47[["Braak"]] = c('6')
Gan_47[["Thal"]] = c('N/A')
Gan_47[["TDP.43"]] = c('N/A')
Gan_47[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_47[["APOE"]] = c('E3E4')
Gan_47[["Age"]] = c('64')
Gan_47[["Age_Onset"]] = c('60')
Gan_47[["Sex"]] = c('M')
Gan_47[["PMI"]] = c('6')
Gan_47[["INDDID"]] = c('122417')
rm(Gan_47.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_47[["percent.mt"]] <- PercentageFeatureSet(object = Gan_47, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_48.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_48/outs/filtered_feature_bc_matrix")
Gan_48 <- CreateSeuratObject(counts = Gan_48.counts, project = "AD_WT_E3E4_M_2", min.cells = 3, min.features = 200)
Gan_48[["Condition"]] = c('AD_WT_E3E4_M')
Gan_48[["Condition_1"]] = c('WT_E3E4')
Gan_48[["Condition_2"]] = c('WT_M')
Gan_48[["Condition_3"]] = c('E3E4_M')
Gan_48[["TREM2"]] = c('WT')
Gan_48[["Dx"]] = c('AD')
Gan_48[["LBD"]] = c('N/A')
Gan_48[["Braak"]] = c('6')
Gan_48[["Thal"]] = c('N/A')
Gan_48[["TDP.43"]] = c('N/A')
Gan_48[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_48[["APOE"]] = c('E3E4')
Gan_48[["Age"]] = c('71')
Gan_48[["Age_Onset"]] = c('62')
Gan_48[["Sex"]] = c('M')
Gan_48[["PMI"]] = c('4.5')
Gan_48[["INDDID"]] = c('111593')
rm(Gan_48.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_48[["percent.mt"]] <- PercentageFeatureSet(object = Gan_48, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_49.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_49/outs/filtered_feature_bc_matrix")
Gan_49 <- CreateSeuratObject(counts = Gan_49.counts, project = "AD_WT_E3E4_M_3", min.cells = 3, min.features = 200)
Gan_49[["Condition"]] = c('AD_WT_E3E4_M')
Gan_49[["Condition_1"]] = c('WT_E3E4')
Gan_49[["Condition_2"]] = c('WT_M')
Gan_49[["Condition_3"]] = c('E3E4_M')
Gan_49[["TREM2"]] = c('WT')
Gan_49[["Dx"]] = c('AD')
Gan_49[["LBD"]] = c('N/A')
Gan_49[["Braak"]] = c('6')
Gan_49[["Thal"]] = c('N/A')
Gan_49[["TDP.43"]] = c('N/A')
Gan_49[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_49[["APOE"]] = c('E3E4')
Gan_49[["Age"]] = c('78')
Gan_49[["Age_Onset"]] = c('64')
Gan_49[["Sex"]] = c('M')
Gan_49[["PMI"]] = c('6.5')
Gan_49[["INDDID"]] = c('107833')
rm(Gan_49.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_49[["percent.mt"]] <- PercentageFeatureSet(object = Gan_49, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_50.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_50/outs/filtered_feature_bc_matrix")
Gan_50 <- CreateSeuratObject(counts = Gan_50.counts, project = "AD_WT_E4E4_M", min.cells = 3, min.features = 200)
Gan_50[["Condition"]] = c('AD_WT_E4E4_M')
Gan_50[["Condition_1"]] = c('WT_E4E4')
Gan_50[["Condition_2"]] = c('WT_M')
Gan_50[["Condition_3"]] = c('E4E4_M')
Gan_50[["TREM2"]] = c('WT')
Gan_50[["Dx"]] = c('AD')
Gan_50[["LBD"]] = c('N/A')
Gan_50[["Braak"]] = c('6')
Gan_50[["Thal"]] = c('N/A')
Gan_50[["TDP.43"]] = c('N/A')
Gan_50[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_50[["APOE"]] = c('E4E4')
Gan_50[["Age"]] = c('78')
Gan_50[["Age_Onset"]] = c('70')
Gan_50[["Sex"]] = c('M')
Gan_50[["PMI"]] = c('6')
Gan_50[["INDDID"]] = c('100436')
rm(Gan_50.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_50[["percent.mt"]] <- PercentageFeatureSet(object = Gan_50, pattern = "^MT-") #recognize mitochondrial transcripts


#for loading Cell Ranger counts:
Gan_51.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_51/outs/filtered_feature_bc_matrix")
Gan_51 <- CreateSeuratObject(counts = Gan_51.counts, project = "AD_R47H_E2E4_F", min.cells = 3, min.features = 200)
Gan_51[["Condition"]] = c('AD_R47H_E2E4_F')
Gan_51[["Condition_1"]] = c('R47H_E2E4')
Gan_51[["Condition_2"]] = c('R47H_F')
Gan_51[["Condition_3"]] = c('E2E4_F')
Gan_51[["TREM2"]] = c('R47H')
Gan_51[["Dx"]] = c('AD')
Gan_51[["LBD"]] = c('N/A')
Gan_51[["Braak"]] = c('6')
Gan_51[["Thal"]] = c('N/A')
Gan_51[["TDP.43"]] = c('N/A')
Gan_51[["ClinicalDx"]] = c('PPA_Logopenic')
Gan_51[["APOE"]] = c('E2E4')
Gan_51[["Age"]] = c('71')
Gan_51[["Age_Onset"]] = c('62')
Gan_51[["Sex"]] = c('F')
Gan_51[["PMI"]] = c('4.5')
Gan_51[["INDDID"]] = c('107342')
rm(Gan_51.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_51[["percent.mt"]] <- PercentageFeatureSet(object = Gan_51, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_52.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_52/outs/filtered_feature_bc_matrix")
Gan_52 <- CreateSeuratObject(counts = Gan_52.counts, project = "AD_R47H_E3E3_F", min.cells = 3, min.features = 200)
Gan_52[["Condition"]] = c('AD_R47H_E3E3_F')
Gan_52[["Condition_1"]] = c('R47H_E3E3')
Gan_52[["Condition_2"]] = c('R47H_F')
Gan_52[["Condition_3"]] = c('E3E3_F')
Gan_52[["TREM2"]] = c('R47H')
Gan_52[["Dx"]] = c('AD_ALS')
Gan_52[["LBD"]] = c('N/A')
Gan_52[["Braak"]] = c('5')
Gan_52[["Thal"]] = c('N/A')
Gan_52[["TDP.43"]] = c('N/A')
Gan_52[["ClinicalDx"]] = c('Amyotrophic lateral sclerosis')
Gan_52[["APOE"]] = c('E3E3')
Gan_52[["Age"]] = c('87')
Gan_52[["Age_Onset"]] = c('63')
Gan_52[["Sex"]] = c('F')
Gan_52[["PMI"]] = c('4')
Gan_52[["INDDID"]] = c('107518')
rm(Gan_52.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_52[["percent.mt"]] <- PercentageFeatureSet(object = Gan_52, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_54.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_54/outs/filtered_feature_bc_matrix")
Gan_54 <- CreateSeuratObject(counts = Gan_54.counts, project = "AD_R47H_E3E3_M", min.cells = 3, min.features = 200)
Gan_54[["Condition"]] = c('AD_R47H_E3E3_M')
Gan_54[["Condition_1"]] = c('R47H_E3E3')
Gan_54[["Condition_2"]] = c('R47H_M')
Gan_54[["Condition_3"]] = c('E3E3_M')
Gan_54[["TREM2"]] = c('R47H')
Gan_54[["Dx"]] = c('AD')
Gan_54[["LBD"]] = c('N/A')
Gan_54[["Braak"]] = c('5')
Gan_54[["Thal"]] = c('N/A')
Gan_54[["TDP.43"]] = c('N/A')
Gan_54[["ClinicalDx"]] = c('Cerebrovascular disease')
Gan_54[["APOE"]] = c('E3E3')
Gan_54[["Age"]] = c('82')
Gan_54[["Age_Onset"]] = c('70')
Gan_54[["Sex"]] = c('M')
Gan_54[["PMI"]] = c('9')
Gan_54[["INDDID"]] = c('104332')
rm(Gan_54.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_54[["percent.mt"]] <- PercentageFeatureSet(object = Gan_54, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_55.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_55/outs/filtered_feature_bc_matrix")
Gan_55 <- CreateSeuratObject(counts = Gan_55.counts, project = "AD_R47H_E3E4_M_1", min.cells = 3, min.features = 200)
Gan_55[["Condition"]] = c('AD_R47H_E3E4_M')
Gan_55[["Condition_1"]] = c('R47H_E3E4')
Gan_55[["Condition_2"]] = c('R47H_M')
Gan_55[["Condition_3"]] = c('E3E4_M')
Gan_55[["TREM2"]] = c('R47H')
Gan_55[["Dx"]] = c('AD')
Gan_55[["LBD"]] = c('N/A')
Gan_55[["Braak"]] = c('6')
Gan_55[["Thal"]] = c('N/A')
Gan_55[["TDP.43"]] = c('N/A')
Gan_55[["ClinicalDx"]] = c('Probable Alzheimers Disease')
Gan_55[["APOE"]] = c('E3E4')
Gan_55[["Age"]] = c('77')
Gan_55[["Age_Onset"]] = c('68')
Gan_55[["Sex"]] = c('M')
Gan_55[["PMI"]] = c('14')
Gan_55[["INDDID"]] = c('113755')
rm(Gan_55.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_55[["percent.mt"]] <- PercentageFeatureSet(object = Gan_55, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_56.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_56/outs/filtered_feature_bc_matrix")
Gan_56 <- CreateSeuratObject(counts = Gan_56.counts, project = "AD_R47H_E3E4_M_2", min.cells = 3, min.features = 200)
Gan_56[["Condition"]] = c('AD_R47H_E3E4_M')
Gan_56[["Condition_1"]] = c('R47H_E3E4')
Gan_56[["Condition_2"]] = c('R47H_M')
Gan_56[["Condition_3"]] = c('E3E4_M')
Gan_56[["TREM2"]] = c('R47H')
Gan_56[["Dx"]] = c('AD')
Gan_56[["LBD"]] = c('N/A')
Gan_56[["Braak"]] = c('6')
Gan_56[["Thal"]] = c('N/A')
Gan_56[["TDP.43"]] = c('N/A')
Gan_56[["ClinicalDx"]] = c('PPA_Logopenic')
Gan_56[["APOE"]] = c('E3E4')
Gan_56[["Age"]] = c('60')
Gan_56[["Age_Onset"]] = c('56')
Gan_56[["Sex"]] = c('M')
Gan_56[["PMI"]] = c('12')
Gan_56[["INDDID"]] = c('100957')
rm(Gan_56.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_56[["percent.mt"]] <- PercentageFeatureSet(object = Gan_56, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_59.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_59/outs/filtered_feature_bc_matrix")
Gan_59 <- CreateSeuratObject(counts = Gan_59.counts, project = "Non_WT_E2E3_F_1", min.cells = 3, min.features = 200)
Gan_59[["Condition"]] = c('Non_WT_E2E3_F')
Gan_59[["Condition_1"]] = c('WT_E2E3')
Gan_59[["Condition_2"]] = c('WT_F')
Gan_59[["Condition_3"]] = c('E2E3_F')
Gan_59[["TREM2"]] = c('WT')
Gan_59[["Dx"]] = c('Normal')
Gan_59[["LBD"]] = c('N/A')
Gan_59[["Braak"]] = c('1')
Gan_59[["Thal"]] = c('N/A')
Gan_59[["TDP.43"]] = c('N/A')
Gan_59[["ClinicalDx"]] = c('N/A')
Gan_59[["APOE"]] = c('E2E3')
Gan_59[["Age"]] = c('72')
Gan_59[["Age_Onset"]] = c('N/A')
Gan_59[["Sex"]] = c('F')
Gan_59[["PMI"]] = c('6')
Gan_59[["INDDID"]] = c('119534')
rm(Gan_59.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_59[["percent.mt"]] <- PercentageFeatureSet(object = Gan_59, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_60.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_60/outs/filtered_feature_bc_matrix")
Gan_60 <- CreateSeuratObject(counts = Gan_60.counts, project = "Non_WT_E2E3_F_2", min.cells = 3, min.features = 200)
Gan_60[["Condition"]] = c('Non_WT_E2E3_F')
Gan_60[["Condition_1"]] = c('WT_E2E3')
Gan_60[["Condition_2"]] = c('WT_F')
Gan_60[["Condition_3"]] = c('E2E3_F')
Gan_60[["TREM2"]] = c('WT')
Gan_60[["Dx"]] = c('Normal')
Gan_60[["LBD"]] = c('N/A')
Gan_60[["Braak"]] = c('1')
Gan_60[["Thal"]] = c('N/A')
Gan_60[["TDP.43"]] = c('N/A')
Gan_60[["ClinicalDx"]] = c('N/A')
Gan_60[["APOE"]] = c('E2E3')
Gan_60[["Age"]] = c('83')
Gan_60[["Age_Onset"]] = c('N/A')
Gan_60[["Sex"]] = c('F')
Gan_60[["PMI"]] = c('3')
Gan_60[["INDDID"]] = c('112090')
rm(Gan_60.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_60[["percent.mt"]] <- PercentageFeatureSet(object = Gan_60, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_61.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_61/outs/filtered_feature_bc_matrix")
Gan_61 <- CreateSeuratObject(counts = Gan_61.counts, project = "Non_WT_NA_F", min.cells = 3, min.features = 200)
Gan_61[["Condition"]] = c('Non_WT_NA_F')
Gan_61[["Condition_1"]] = c('WT_NA')
Gan_61[["Condition_2"]] = c('WT_F')
Gan_61[["Condition_3"]] = c('NA_F')
Gan_61[["TREM2"]] = c('WT')
Gan_61[["Dx"]] = c('Normal')
Gan_61[["LBD"]] = c('N/A')
Gan_61[["Braak"]] = c('2')
Gan_61[["Thal"]] = c('N/A')
Gan_61[["TDP.43"]] = c('N/A')
Gan_61[["ClinicalDx"]] = c('N/A')
Gan_61[["APOE"]] = c('NA')
Gan_61[["Age"]] = c('75')
Gan_61[["Age_Onset"]] = c('N/A')
Gan_61[["Sex"]] = c('F')
Gan_61[["PMI"]] = c('12')
Gan_61[["INDDID"]] = c('125061')
rm(Gan_61.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_61[["percent.mt"]] <- PercentageFeatureSet(object = Gan_61, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_62.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_62/outs/filtered_feature_bc_matrix")
Gan_62 <- CreateSeuratObject(counts = Gan_62.counts, project = "Non_WT_E2E3_M_1", min.cells = 3, min.features = 200)
Gan_62[["Condition"]] = c('Non_WT_E2E3_M')
Gan_62[["Condition_1"]] = c('WT_E2E3')
Gan_62[["Condition_2"]] = c('WT_M')
Gan_62[["Condition_3"]] = c('E2E3_M')
Gan_62[["TREM2"]] = c('WT')
Gan_62[["Dx"]] = c('Normal')
Gan_62[["LBD"]] = c('N/A')
Gan_62[["Braak"]] = c('0')
Gan_62[["Thal"]] = c('N/A')
Gan_62[["TDP.43"]] = c('N/A')
Gan_62[["ClinicalDx"]] = c('N/A')
Gan_62[["APOE"]] = c('E2E3')
Gan_62[["Age"]] = c('68')
Gan_62[["Age_Onset"]] = c('N/A')
Gan_62[["Sex"]] = c('M')
Gan_62[["PMI"]] = c('14')
Gan_62[["INDDID"]] = c('118709')
rm(Gan_62.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_62[["percent.mt"]] <- PercentageFeatureSet(object = Gan_62, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_63.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_63/outs/filtered_feature_bc_matrix")
Gan_63 <- CreateSeuratObject(counts = Gan_63.counts, project = "Non_WT_E2E3_M_2", min.cells = 3, min.features = 200)
Gan_63[["Condition"]] = c('Non_WT_E2E3_M')
Gan_63[["Condition_1"]] = c('WT_E2E3')
Gan_63[["Condition_2"]] = c('WT_M')
Gan_63[["Condition_3"]] = c('E2E3_M')
Gan_63[["TREM2"]] = c('WT')
Gan_63[["Dx"]] = c('Normal')
Gan_63[["LBD"]] = c('N/A')
Gan_63[["Braak"]] = c('1')
Gan_63[["Thal"]] = c('N/A')
Gan_63[["TDP.43"]] = c('N/A')
Gan_63[["ClinicalDx"]] = c('N/A')
Gan_63[["APOE"]] = c('E2E3')
Gan_63[["Age"]] = c('72')
Gan_63[["Age_Onset"]] = c('N/A')
Gan_63[["Sex"]] = c('M')
Gan_63[["PMI"]] = c('17')
Gan_63[["INDDID"]] = c('120927')
rm(Gan_63.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_63[["percent.mt"]] <- PercentageFeatureSet(object = Gan_63, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_64.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_64/outs/filtered_feature_bc_matrix")
Gan_64 <- CreateSeuratObject(counts = Gan_64.counts, project = "Non_WT_E2E4_M", min.cells = 3, min.features = 200)
Gan_64[["Condition"]] = c('Non_WT_E2E4_M')
Gan_64[["Condition_1"]] = c('WT_E2E4')
Gan_64[["Condition_2"]] = c('WT_M')
Gan_64[["Condition_3"]] = c('E2E4_M')
Gan_64[["TREM2"]] = c('WT')
Gan_64[["Dx"]] = c('Normal')
Gan_64[["LBD"]] = c('N/A')
Gan_64[["Braak"]] = c('0')
Gan_64[["Thal"]] = c('N/A')
Gan_64[["TDP.43"]] = c('N/A')
Gan_64[["ClinicalDx"]] = c('N/A')
Gan_64[["APOE"]] = c('E2E4')
Gan_64[["Age"]] = c('61')
Gan_64[["Age_Onset"]] = c('N/A')
Gan_64[["Sex"]] = c('M')
Gan_64[["PMI"]] = c('6')
Gan_64[["INDDID"]] = c('100786')
rm(Gan_64.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_64[["percent.mt"]] <- PercentageFeatureSet(object = Gan_64, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_65.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_65/outs/filtered_feature_bc_matrix")
Gan_65 <- CreateSeuratObject(counts = Gan_65.counts, project = "Non_WT_E3E3_M_1", min.cells = 3, min.features = 200)
Gan_65[["Condition"]] = c('Non_WT_E3E3_M')
Gan_65[["Condition_1"]] = c('WT_E3E3')
Gan_65[["Condition_2"]] = c('WT_M')
Gan_65[["Condition_3"]] = c('E3E3_M')
Gan_65[["TREM2"]] = c('WT')
Gan_65[["Dx"]] = c('Normal')
Gan_65[["LBD"]] = c('N/A')
Gan_65[["Braak"]] = c('1')
Gan_65[["Thal"]] = c('N/A')
Gan_65[["TDP.43"]] = c('N/A')
Gan_65[["ClinicalDx"]] = c('N/A')
Gan_65[["APOE"]] = c('E3E3')
Gan_65[["Age"]] = c('75')
Gan_65[["Age_Onset"]] = c('N/A')
Gan_65[["Sex"]] = c('M')
Gan_65[["PMI"]] = c('17')
Gan_65[["INDDID"]] = c('113464')
rm(Gan_65.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_65[["percent.mt"]] <- PercentageFeatureSet(object = Gan_65, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
Gan_66.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_UPenn/cellranger/Gan_66/outs/filtered_feature_bc_matrix")
Gan_66 <- CreateSeuratObject(counts = Gan_66.counts, project = "Non_WT_E3E3_M_2", min.cells = 3, min.features = 200)
Gan_66[["Condition"]] = c('Non_WT_E3E3_M')
Gan_66[["Condition_1"]] = c('WT_E3E3')
Gan_66[["Condition_2"]] = c('WT_M')
Gan_66[["Condition_3"]] = c('E3E3_M')
Gan_66[["TREM2"]] = c('WT')
Gan_66[["Dx"]] = c('Normal')
Gan_66[["LBD"]] = c('N/A')
Gan_66[["Braak"]] = c('2')
Gan_66[["Thal"]] = c('N/A')
Gan_66[["TDP.43"]] = c('N/A')
Gan_66[["ClinicalDx"]] = c('N/A')
Gan_66[["APOE"]] = c('E3E3')
Gan_66[["Age"]] = c('83')
Gan_66[["Age_Onset"]] = c('N/A')
Gan_66[["Sex"]] = c('M')
Gan_66[["PMI"]] = c('6')
Gan_66[["INDDID"]] = c('119767')
rm(Gan_66.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
Gan_66[["percent.mt"]] <- PercentageFeatureSet(object = Gan_66, pattern = "^MT-") #recognize mitochondrial transcripts

setwd("/athena/ganlab/scratch/lif4001/Human_UPenn/data_analysis/DF_2ndRound")
###############################################################################################
all <- Gan_43
pdf("Gan_43_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_43_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_43_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Gan_43_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

#saveRDS(all,"Gan_43_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
#all<-readRDS("Gan_43_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_43_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6582) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.24, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.24_303", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_43_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_43_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.24_303" #visualizing the singlet vs doublet cells
pdf("Gan_43_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_43_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_43_singlets.rds")
singlets<-readRDS("Gan_43_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_43_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_43_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_43_singlets_PCA.rds")

###############################################################################################
all <- Gan_44
pdf("Gan_44_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_44_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_44_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Gan_44_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

#saveRDS(all,"Gan_44_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
#all<-readRDS("Gan_44_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_44_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6516) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_300", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_44_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_44_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_300" #visualizing the singlet vs doublet cells
pdf("Gan_44_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_44_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_44_singlets.rds")
singlets<-readRDS("Gan_44_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_44_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_44_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_44_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_45
pdf("Gan_45_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_45_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)

#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_45_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("Gan_45_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

#saveRDS(all,"Gan_45_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
#all<-readRDS("Gan_45_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_45_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6648) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.27, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.27, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.27_306", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_45_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_45_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.27_306" #visualizing the singlet vs doublet cells
pdf("Gan_45_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_45_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_45_singlets.rds")
singlets<-readRDS("Gan_45_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_45_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_45_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_45_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_46
pdf("Gan_46_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_46_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_46_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_46_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_46_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5649) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_220", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_46_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_46_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_220" #visualizing the singlet vs doublet cells
pdf("Gan_46_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_46_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_46_singlets.rds")
singlets<-readRDS("Gan_46_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_46_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_46_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_46_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_47
pdf("Gan_47_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_47_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_47_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_47_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_47_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8564) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.03, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.03_522", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_47_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_47_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.03_522" #visualizing the singlet vs doublet cells
pdf("Gan_47_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_47_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_47_singlets.rds")
singlets<-readRDS("Gan_47_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_47_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_47_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_47_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_48
pdf("Gan_48_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_48_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_48_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_48_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_48_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7751) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.16, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.16, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.16_419", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_48_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_48_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.16_419" #visualizing the singlet vs doublet cells
pdf("Gan_48_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_48_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_48_singlets.rds")
singlets<-readRDS("Gan_48_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_48_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_48_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_48_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_49
pdf("Gan_49_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_49_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_49_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_49_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_49_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5548) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.25, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.25, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.25_216", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_49_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_49_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.25_216" #visualizing the singlet vs doublet cells
pdf("Gan_49_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_49_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_49_singlets.rds")
singlets<-readRDS("Gan_49_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_49_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_49_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_49_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_50
pdf("Gan_50_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_50_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_50_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_50_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_50_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*5488) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_214", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_50_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_50_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_214" #visualizing the singlet vs doublet cells
pdf("Gan_50_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_50_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_50_singlets.rds")
singlets<-readRDS("Gan_50_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_50_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_50_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_50_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_51
pdf("Gan_51_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_51_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_51_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_51_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_51_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8028) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.11, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.11, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.11_490", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_51_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_51_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.11_490" #visualizing the singlet vs doublet cells
pdf("Gan_51_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_51_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_51_singlets.rds")
singlets<-readRDS("Gan_51_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_51_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_51_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_51_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_52
pdf("Gan_52_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_52_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_52_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_52_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_52_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7528) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_407", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_52_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_52_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_407" #visualizing the singlet vs doublet cells
pdf("Gan_52_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_52_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_52_singlets.rds")
singlets<-readRDS("Gan_52_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_52_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_52_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_52_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_54
pdf("Gan_54_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_54_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_54_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_54_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_54_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7219) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.22, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.22, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.22_390", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_54_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_54_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.22_390" #visualizing the singlet vs doublet cells
pdf("Gan_54_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_54_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_54_singlets.rds")
singlets<-readRDS("Gan_54_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_54_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_54_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_54_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_55
pdf("Gan_55_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_55_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_55_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_55_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_55_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7132) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_385", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_55_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_55_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_385" #visualizing the singlet vs doublet cells
pdf("Gan_55_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_55_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_55_singlets.rds")
singlets<-readRDS("Gan_55_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_55_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_55_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_55_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_56
pdf("Gan_56_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_56_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_56_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_56_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_56_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.023*3681) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.3, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_85", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_56_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_56_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.3_85" #visualizing the singlet vs doublet cells
pdf("Gan_56_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_56_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_56_singlets.rds")
singlets<-readRDS("Gan_56_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_56_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_56_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_56_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_59
pdf("Gan_59_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_59_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_59_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_59_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_59_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7225) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.24, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.24_390", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_59_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_59_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.24_390" #visualizing the singlet vs doublet cells
pdf("Gan_59_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_59_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_59_singlets.rds")
singlets<-readRDS("Gan_59_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_59_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_59_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_59_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_60
pdf("Gan_60_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_60_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_60_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_60_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_60_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8735) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_533", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_60_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_60_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_533" #visualizing the singlet vs doublet cells
pdf("Gan_60_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_60_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_60_singlets.rds")
singlets<-readRDS("Gan_60_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_60_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_60_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_60_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_61
pdf("Gan_61_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_61_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_61_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_61_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_61_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.023*3119) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_72", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_61_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_61_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_72" #visualizing the singlet vs doublet cells
pdf("Gan_61_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_61_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_61_singlets.rds")
singlets<-readRDS("Gan_61_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_61_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_61_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_61_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_62
pdf("Gan_62_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_62_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_62_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_62_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_62_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8825) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_538", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_62_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_62_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_538" #visualizing the singlet vs doublet cells
pdf("Gan_62_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_62_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_62_singlets.rds")
singlets<-readRDS("Gan_62_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_62_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_62_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_62_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_63
pdf("Gan_63_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_63_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_63_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_63_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_63_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6150) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_283", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_63_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_63_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_283" #visualizing the singlet vs doublet cells
pdf("Gan_63_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_63_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_63_singlets.rds")
singlets<-readRDS("Gan_63_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_63_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_63_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_63_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_64
pdf("Gan_64_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_64_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_64_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_64_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_64_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8047) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_491", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_64_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_64_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_491" #visualizing the singlet vs doublet cells
pdf("Gan_64_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_64_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_64_singlets.rds")
singlets<-readRDS("Gan_64_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_64_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_64_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_64_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_65
pdf("Gan_65_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_65_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_65_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_65_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_65_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7107) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_384", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_65_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_65_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_384" #visualizing the singlet vs doublet cells
pdf("Gan_65_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_65_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_65_singlets.rds")
singlets<-readRDS("Gan_65_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_65_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_65_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_65_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- Gan_66
pdf("Gan_66_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Gan_66_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Gan_66_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_66_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Gan_66_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8121) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.22, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.22, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.22_495", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Gan_66_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("Gan_66_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.22_495" #visualizing the singlet vs doublet cells
pdf("Gan_66_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"Gan_66_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"Gan_66_singlets.rds")
singlets<-readRDS("Gan_66_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("Gan_66_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("Gan_66_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"Gan_66_singlets_PCA.rds")

###############################################################################################
###############################################################################################
