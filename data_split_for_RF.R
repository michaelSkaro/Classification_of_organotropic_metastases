#!/usr/bin/env Rscript

# The objective of this script will be to put out data into our classifier function. This is 
# a feeder Rscript that will conduct all of the data pre-processing prior to 
# training. This will includer reading in raw data, annotating the data
# concatenating annotations as the labels, and finally prepairing these data
# for learning by writing them to the appropriate file directory. 


# invoke libraries. 

library("randomForest")
library("stringr")
library("tidyverse")
library("dplyr")
library("caret")
library("doMC")

# Projects i Need to add
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

# add the clinical annotations of progression

# read in the files in one by one and select the correct columns. 
# concatenate the files and end with one uniform progression file. 

tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv") %>% 
  tibble::column_to_rownames(var ="V1")


proj <- projects[1]

for(proj in projects){
setwd("~/storage/PanCancerAnalysis/ML_2019")


# read in data   
dat <- data.table::fread(str_glue("{proj}_learning.csv", header = TRUE))
colnames(dat)[1] <- "barcode"

# remove the normal samples, in this section we are only classifying the metastatic cases from 
# here on out

dat <- dat[dat$sample_type == "Primary_Tumor",]


# add column annotation for the metastatic staus of the sample's metastatic status

dat<- left_join(dat, tumor.samples, by="barcode")

# get rid of some of the unnecessary columns that were added, we want the majority
# of the data to be numeric with only a few feactorial variables.

dat <- dat %>%
  dplyr::select(-barcode,-sample_type.x,-sample_type.y, -project, -tumor_stage, -caseID, -fileID, -filename, 
                -barcode_short, -bcr_patient_barcode,-Metastatic_status, -LymphNodeStatus)


# change name of Met_loc to just labels


dat <- dat %>%
  dplyr::rename(lables = met_loc)


write.csv(dat, file = str_glue("~/storage/PanCancerAnalysis/ML_2019/RF_input/{proj}_met_loc.csv"))


}

