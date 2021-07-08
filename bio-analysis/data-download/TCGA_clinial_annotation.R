# Towards annotation of all TCGA metastatic loci
library(TCGAbiolinks)
library(xml2)
library(tidyverse)

projects <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC",
              "TCGA-CHOL","TCGA-COAD","TCGA-DLBC","TCGA-ESCA",
              "TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC",
              "TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC",
              "TCGA-LUAD","TCGA-LUSC","TCGA-MESO","TCGA-OV",
              "TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ",
              "TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT",
              "TCGA-THCA","TCGA-THYM","TCGA-UCEC","TCGA-UCS","TCGA-UVM")

setwd("/mnt/storage/mskaro1/Clinical_annotation")

ACC <- function(proj){
    dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), 
                             header = TRUE,na.strings=c("","NA"))
    dat1 <- dat %>%
      dplyr::select(fileID,ct_scan_findings,distant_metastasis_anatomic_site,
                    metastatic_neoplasm_initial_diagnosis_anatomic_site,
                  `metastatic_neoplasm_initial_diagnosis_anatomic_site[1]`,
                  `metastatic_neoplasm_initial_diagnosis_anatomic_site[2]`,
                  `metastatic_neoplasm_initial_diagnosis_anatomic_site[3]`,
                  new_neoplasm_event_occurrence_anatomic_site,
                  new_neoplasm_occurrence_anatomic_site_text,
                  number_of_lymphnodes_positive_by_he,other_malignancy_anatomic_site)
    write.csv(dat1,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
    met_anno <- data.table::fread("met_anno/Temp_TCGA-ACC_met_anno.txt")
    dat <- left_join(dat, met_anno, by="fileID")
    write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
} # Complete
BLCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE,na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(fileID,metastatic_site,
                  `metastatic_site[1]`,`metastatic_site[2]`,`metastatic_site[3]`,
                  new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text, 
                  other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,other_metastatic_site) 
  write.csv(dat1,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/Temp_{proj}_met_anno.txt")) %>%
    dplyr::select(c(fileID,Metastatic_site))
    
  
  dat <- left_join(dat,met_anno, by ="fileID")
  write.csv(df, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
  # when you read it into python to do the analysis you can use this line.
  #import pandas as pd
  #new_df = pd.concat([df.drop('Issues', 1), df['Issues'].str.get_dummies(sep=",")], 1)
  
  
} # Complete
BRCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"),header = TRUE,na.strings=c("","NA"))
  dat1 <- dat %>%
    dplyr::select(
      fileID,bcr_patient_barcode,metastatic_site_at_diagnosis,metastatic_site_at_diagnosis_other,
      `metastatic_site_at_diagnosis[1]`,
      `metastatic_site_at_diagnosis[2]`,
      `metastatic_site_at_diagnosis[3]`,
      `metastatic_site_at_diagnosis[4]`,
      new_neoplasm_event_occurrence_anatomic_site,
      new_tumor_event_after_initial_treatment,
      other_malignancy_anatomic_site
    ) %>%
    unite("Loci", metastatic_site_at_diagnosis:other_malignancy_anatomic_site, sep= ",", 
          remove = FALSE)
  dat1 <- dat1 %>% dplyr::select(fileID,Loci)
  write.csv(met_anno,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
  met_anno <- data.table::fread("met_anno/Temp_TCGA-BRCA_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  
  df <- left_join(dat,met_anno, by ="fileID")
  
  
  # When we get things more clean here
  #met_anno <- cbind(met_anno, qdapTools::mtabulate(strsplit(met_anno$Loci, ",")))
  write.csv(df,str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
  
  
} # Complete
CESC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE,na.strings=c("","NA"))
  dat1 <-  dat %>% dplyr::select(
    fileID,
    bcr_patient_barcode,
    diagnostic_ct_result_outcome,
    `diagnostic_ct_result_outcome[1]`,
    `diagnostic_ct_result_outcome[2]`,
    `diagnostic_ct_result_outcome[3]`,
    `diagnostic_ct_result_outcome[4]`,
    `diagnostic_ct_result_outcome[5]`,
    `diagnostic_ct_result_outcome[6]`,
    diagnostic_mri_result_outcome,
    `diagnostic_mri_result_outcome[1]`,
    `diagnostic_mri_result_outcome[2]`,
    `diagnostic_mri_result_outcome[3]`,
    `diagnostic_mri_result_outcome[4]`,
    `diagnostic_mri_result_outcome[5]`,
    `diagnostic_mri_result_outcome[6]`,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    new_neoplasm_event_type,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  write.csv(dat1,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  colnames(met_anno) <- c("fileID","Metastatic_site")
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  

} # Complete
CHOL <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), header = TRUE, na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  write.csv(dat1,str_glue("met_anno/{proj}_met_anno.txt"))
  df <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  dat <- left_join(dat,df, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
} # Too small of a data set
COAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    new_neoplasm_event_type,
    non_nodal_tumor_deposits,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  )
  write.csv(dat1,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
  df <- data.table::fread(str_glue("met_anno/Temp_{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- left_join(dat,df, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
} # Complete
DLBC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(
    fileID,
    bone_marrow_involvement,
    lymph_node_involvement_site,
    `lymph_node_involvement_site[1]`,
    `lymph_node_involvement_site[2]`,
    `lymph_node_involvement_site[3]`,
    `lymph_node_involvement_site[4]`,
    `lymph_node_involvement_site[5]`,
    `lymph_node_involvement_site[6]`,
    `lymph_node_involvement_site[7]`,
    `lymph_node_involvement_site[8]`,
    `lymph_node_involvement_site[9]`,
    `lymph_node_involvement_site[10]`,
    tumor_tissue_site,
    `tumor_tissue_site[1]`,
    `tumor_tissue_site[2]`,
    `tumor_tissue_site[3]`,
    `tumor_tissue_site[5]`,
    `tumor_tissue_site[6]`,
    `tumor_tissue_site[6]`
  )
  write.csv(dat1,str_glue("met_anno/Temp_{proj}_met_anno.txt"))
  met_anno <- data.table::fread(str_glue("met_anno/{proj}_met_anno.txt", na.strings=c("","NA"))) %>%
    dplyr::select(fileID, Metastasis)
  dat <- left_join(dat,met_anno, by ="fileID")
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
} # Complete
ESCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    number_of_lymphnodes_positive_by_he,
    number_of_lymphnodes_positive_by_ihc
  ))
  met_anno <- data.table::fread("met_anno/Temp_TCGA-ESCA_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by ="fileID")
  
  write.csv(dat, str_glue("/mnt/storage/mskaro1/Clinical_annotation/met_anno/Complete/Clinical_annotation_metastatic_locations_{proj}.csv"))
  
} # Not enough documented locations
GBM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>%
    dplyr::select(c(fileID, new_neoplasm_event_type,other_malignancy_anatomic_site,new_neoplasm_event_type,new_tumor_event_after_initial_treatment,tumor_tissue_site))
  # stop
} # Too small
HNSC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  dat1 <- dat %>%
    dplyr::select(c(fileID, new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,number_of_lymphnodes_positive_by_he,
                    number_of_lymphnodes_positive_by_ihc,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text))
  #write.csv(dat1, "met_anno/Temp_TCGA-HNSC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-HNSC_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")  
  write.csv(dat, "met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-HNSC.csv")
  
} # Complete
KICH <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  # stop 
  
} # Too small
KIRC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    #new_neoplasm_event_occurrence_anatomic_site,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    number_of_lymphnodes_positive
    #number_of_lymphnodes_positive_by_ihc
  ))
  #write.csv(dat1, "met_anno/Temp_TCGA-KIRC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-KIRC_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat, "met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-KIRC.csv")
  
} # Complete
KIRP <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    #new_neoplasm_event_occurrence_anatomic_site,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    number_of_lymphnodes_positive
    #number_of_lymphnodes_positive_by_ihc
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-KIRP_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-KIRP_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  write.csv(dat, "met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-KIRP.csv")
  
} # Complete
LAML <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  # Stop
  
} # Not annotation
LGG <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    #new_neoplasm_event_occurrence_anatomic_site,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    #number_of_lymphnodes_positive_by_ihc
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-LGG_met_anno.txt")
  # too few annotated metastases
  
  
} # Too few annotated Mets
LIHC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-LIHC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-LIHC_met_anno.txt") %>%
    dplyr::select(c(fileID, Metastatic_site))
  dat <- dplyr::left_join(dat, met_anno, by = "fileID")
  write.csv(dat1, "met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-LIHC.csv")
} # Complete
LUAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-LUAD_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-LUAD_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-LUAD.csv")
  
} # Complete
LUSC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-LUSC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-LUSC_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-LUSC.csv")
  
} # Complete
MESO <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-MESO_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-MESO_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-MESO.csv")
  
  
} # Complete
OV <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    other_malignancy_anatomic_site,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-OV_met_anno.txt")
  
  #Stop
  
  
  
} # Too small
PAAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-PAAD_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-PAAD_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-PAAD.csv")
  
  
} # Complete
PCPG <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-PCPG_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-PCPG_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-PCPG.csv")
  
  
} # Something Wrong
PRAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-PRAD_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-PRAD_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-PRAD.csv")
  
  
} # Something Wrong
READ <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    lymphatic_invasion,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-READ_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-READ_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-READ.csv")
  
  
} # Complete
SARC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
  fileID,
  new_neoplasm_event_occurrence_anatomic_site,
  new_neoplasm_occurrence_anatomic_site_text,
  other_malignancy_anatomic_site,
  other_malignancy_anatomic_site,
  other_malignancy_anatomic_site_text,
  `tumor_tissue_site[1]`,
  `tumor_tissue_site[2]`,
  `tumor_tissue_site[3]`
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-SARC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-SARC_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-SARC.csv")
  
} # Complete
SKCM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_tumor_metastasis_anatomic_site,
    new_tumor_metastasis_anatomic_site_other_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-SKCM_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-SKCM_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-SKCM.csv")
  
} # Complete
STAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    number_of_lymphnodes_positive_by_he,
    other_malignancy_anatomic_site,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-STAD_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-STAD_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-STAD.csv")
  
  
  
} # Complete
TGCT <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-TGCT_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-TGCT_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-TGCT.csv")
  
} # Complete
THCA <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    metastatic_neoplasm_confirmed_diagnosis_method_name,
    `metastatic_neoplasm_confirmed_diagnosis_method_name[1]`,
    `metastatic_neoplasm_confirmed_diagnosis_method_name[2]`,
    metastatic_site,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site_text,
    number_of_lymphnodes_positive_by_he,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-THCA_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-THCA_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-THCA.csv")
  
} # Complete
THYM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-THYM_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-THYM_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-THYM.csv")
  
} # Complete
UCEC <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-UCEC_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-UCEC_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-UCEC.csv")
} # Complete
UCS <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-UCS_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-UCS_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-UCS.csv")
  
} # Complete
UVM <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    pathologic_N
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-UVM_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-UVM_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-UVM.csv")
} # Complete
