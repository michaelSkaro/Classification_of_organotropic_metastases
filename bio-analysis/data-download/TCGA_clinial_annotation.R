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
  
} # Too small 
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
  
} # Complete
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
  
} # No annotation
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
  
  
} # Complete
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
  
  
  
} # Complete
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
  
  
} # Complete
PRAD <- function(proj){
  dat <- data.table::fread(str_glue("Clinical_annotation_{proj}.csv"), na.strings=c("","NA"))
  dat1 <- dat %>% dplyr::select(c(
    fileID,
    number_of_lymphnodes_positive_by_he,
    pathologic_N,
    new_neoplasm_event_occurrence_anatomic_site,
    new_neoplasm_occurrence_anatomic_site_text,
    other_malignancy_anatomic_site,
    other_malignancy_anatomic_site_text
  ))
  write.csv(dat1, "met_anno/Temp_TCGA-PRAD_met_anno.txt")
  met_anno <- data.table::fread("met_anno/Temp_TCGA-PRAD_met_anno.txt") %>%
    dplyr::select(c(fileID,Metastatic_site))
  dat <- dplyr::left_join(dat,met_anno, by = "fileID")
  write.csv(dat,"met_anno/Complete/Clinical_annotation_metastatic_locations_TCGA-PRAD.csv")
  
  
} # Complete
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


met_anno_all <- data.table::fread("Annotated_met_loc_all.txt")
df <- cbind(met_anno_all, mtabulate(strsplit(met_anno_all$Metastastic_site, ","))) %>%
  select(-V1)


met_anno_all_long <- met_anno_all %>%
  tidyr::separate_rows(Metastastic_site, sep = ",", convert = TRUE)
colnames(met_anno_all_long) <- c("fileID","Metastatic_site","CT")

counts <- ddply(met_anno_all_long, .(Metastatic_site, CT), nrow)
names(counts) <- c("Metastatic_site", "CT", "Freq")
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}

## transform all columns
counts <- counts %>% mutate_each(funs(empty_as_na)) 
counts$log.abundance <- log(counts$Freq) 
counts <- counts %>%
  dplyr::filter(CT != "OV") %>%
  dplyr::filter(CT !="LGG") %>%
  dplyr::filter(CT !="TGCT") %>%
  dplyr::filter(CT !="THYM") %>%
  dplyr::mutate(CT = replace(CT,CT == "READ", "COADREAD")) %>%
  dplyr::mutate(CT = replace(CT,CT == "COAD", "COADREAD")) %>%
  #dplyr::filter(Freq >8) %>%
  dplyr::arrange(CT) %>%
  dplyr::select(-log.abundance) %>%
  dplyr::group_by(Metastatic_site, CT)
counts <- counts[c(2,1,3)]
sites <- as.data.frame(unique(counts$Metastatic_site))
colnames(sites) <- "sites"
write.csv(sites, "sites.csv")

dat <- data.table::fread("Metastasis_by_loc_by_cancer.csv", sep = ",", header = TRUE) %>%
  dplyr::select(-V1) %>%
  pivot_longer(!Metastatic_site, names_to = "CT", values_to = "count") %>%
  dplyr::filter(CT != "OV") %>%
  dplyr::filter(CT !="LGG") %>%
  dplyr::filter(CT !="TGCT") %>%
  dplyr::filter(CT !="THYM") %>%
  dplyr::mutate(CT = replace(CT,CT == "READ", "COADREAD")) %>%
  dplyr::mutate(CT = replace(CT,CT == "COAD", "COADREAD"))
dat$log.abundance <- log(dat$count)
dat$log.abundance[dat$log.abundance<0] <- 0


p1 <- ggplot(data = dat, mapping = aes(x = CT,y = Metastatic_site, fill =log.abundance)) +
  geom_raster() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(angle = 30, size = 5, face = "bold")) +
  ylab(label = "Metatatic Sites") + 
  xlab(label  ="Cancer Types") +
  ggtitle("Site specific progression")+
  scale_fill_gradient(name = "Log Frequency",
                      low = "Black",
                      high = "Red")


dat <- data.table::fread("Metastasis_by_loc_by_cancer.csv", sep = ",", header = TRUE) %>%
  dplyr::select(-V1) %>%
  pivot_longer(!Metastatic_site, names_to = "CT", values_to = "count") %>%
  dplyr::filter(CT != "OV") %>%
  dplyr::filter(CT !="LGG") %>%
  dplyr::filter(CT !="TGCT") %>%
  dplyr::filter(CT !="THYM") %>%
  dplyr::mutate(CT = replace(CT,CT == "READ", "COADREAD")) %>%
  dplyr::mutate(CT = replace(CT,CT == "COAD", "COADREAD")) %>%
  #dplyr::filter(count>=8) %>%
  dplyr::filter(Metastatic_site %in% selected)
dat$log.abundance <- log(dat$count)
dat$log.abundance[dat$log.abundance<0] <- 0

ggplot(data = dat, mapping = aes(x = CT,y = Metastatic_site, fill =log.abundance)) +
  geom_raster() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.y = element_text(angle = 30, size = 5, face = "bold")) +
  ylab(label = "Metatatic Sites") + 
  xlab(label  ="Cancer Types") +
  ggtitle("Selected Locations")+
  scale_fill_gradient(name = "Log Frequency",
                      low = "Black",
                      high = "Red")





# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods  
# [10] base     
# 
# other attached packages:
#   [1] cowplot_1.1.1            DOSE_3.16.0              clusterProfiler_3.18.1  
# [4] msigdbr_7.2.1            org.Hs.eg.db_3.12.0      AnnotationDbi_1.52.0    
# [7] IRanges_2.24.1           S4Vectors_0.28.1         circlize_0.4.12         
# [10] ComplexHeatmap_2.6.2     xml2_1.3.2               TCGAbiolinks_2.18.0     
# [13] forcats_0.5.1            stringr_1.4.0            dplyr_1.0.5             
# [16] purrr_0.3.4              readr_1.4.0              tidyr_1.1.3             
# [19] tibble_3.1.1             tidyverse_1.3.1          flowAI_1.20.1           
# [22] MIMOSA_1.28.1            Biobase_2.50.0           reshape_0.8.8           
# [25] plyr_1.8.6               MASS_7.3-53              flowMerge_2.38.0        
# [28] snow_0.4-3               foreach_1.5.1            Rgraphviz_2.34.0        
# [31] feature_1.2.15           graph_1.68.0             BiocGenerics_0.36.1     
# [34] flowTrans_1.42.0         flowClust_3.28.0         COMPASS_1.28.0          
# [37] flowViz_1.54.0           lattice_0.20-41          cytolib_2.2.1           
# [40] openCyto_2.2.0           CytoML_2.2.2             flowStats_4.2.0         
# [43] ggcyto_1.18.0            ncdfFlow_2.36.0          BH_1.75.0-0             
# [46] RcppArmadillo_0.10.4.0.0 ggplot2_3.3.3            flowWorkspace_4.2.0     
# [49] flowCore_2.2.0          
# 
# loaded via a namespace (and not attached):
#   [1] rappdirs_0.3.3              SparseM_1.81                R.methodsS3_1.8.1          
# [4] coda_0.19-4                 bit64_4.0.5                 knitr_1.33                 
# [7] DelayedArray_0.16.3         R.utils_2.10.1              data.table_1.14.0          
# [10] rpart_4.1-15                RCurl_1.98-1.3              generics_0.1.0             
# [13] timeSeries_3062.100         RSQLite_2.2.7               shadowtext_0.0.8           
# [16] enrichplot_1.10.2           bit_4.0.4                   lubridate_1.7.10           
# [19] SummarizedExperiment_1.20.0 assertthat_0.2.1            viridis_0.6.0              
# [22] xfun_0.22                   fBasics_3042.89.1           hms_1.0.0                  
# [25] evaluate_0.14               DEoptimR_1.0-8              fansi_0.4.2                
# [28] progress_1.2.2              dbplyr_2.1.1                readxl_1.3.1               
# [31] geneplotter_1.68.0          igraph_1.2.6                DBI_1.1.1                  
# [34] tmvnsim_1.0-2               mcmc_0.9-7                  ellipsis_0.3.1             
# [37] ks_1.12.0                   backports_1.2.1             annotate_1.68.0            
# [40] MCMCpack_1.5-0              RcppParallel_5.1.2          biomaRt_2.46.3             
# [43] MatrixGenerics_1.2.1        vctrs_0.3.7                 quantreg_5.85              
# [46] Cairo_1.5-12.2              abind_1.4-5                 cachem_1.0.4               
# [49] withr_2.4.2                 ggforce_0.3.3               aws.signature_0.6.0        
# [52] robustbase_0.93-7           prettyunits_1.1.1           mclust_5.4.7               
# [55] mnormt_2.0.2                cluster_2.1.0               crayon_1.4.1               
# [58] genefilter_1.72.1           ellipse_0.4.2               pkgconfig_2.0.3            
# [61] tweenr_1.0.2                GenomeInfoDb_1.26.7         changepoint_2.2.2          
# [64] rlang_0.4.10                spatial_7.3-13              lifecycle_1.0.0            
# [67] MatrixModels_0.5-0          downloader_0.4              BiocFileCache_1.14.0       
# [70] modelr_0.1.8                polyclip_1.10-0             cellranger_1.1.0           
# [73] tcltk_4.0.3                 matrixStats_0.58.0          Matrix_1.3-2               
# [76] zoo_1.8-9                   reprex_2.0.0                base64enc_0.1-3            
# [79] GlobalOptions_0.1.2         viridisLite_0.4.0           png_0.1-7                  
# [82] rjson_0.2.20                stabledist_0.7-1            bitops_1.0-7               
# [85] R.oo_1.24.0                 KernSmooth_2.23-18          blob_1.2.1                 
# [88] shape_1.4.5                 qvalue_2.22.0               jpeg_0.1-8.1               
# [91] aws.s3_0.3.21               scales_1.1.1                memoise_2.0.0              
# [94] magrittr_2.0.1              hexbin_1.28.2               zlibbioc_1.36.0            
# [97] scatterpie_0.1.6            compiler_4.0.3              hdrcde_3.4                 
# [100] RColorBrewer_1.1-2          clue_0.3-59                 DESeq2_1.30.1              
# [103] rrcov_1.5-5                 cli_2.4.0                   XVector_0.30.0             
# [106] Formula_1.2-4               tidyselect_1.1.0            stringi_1.5.3              
# [109] TCGAbiolinksGUI.data_1.10.0 RProtoBufLib_2.2.0          yaml_2.2.1                 
# [112] GOSemSim_2.16.1             locfit_1.5-9.4              askpass_1.1                
# [115] ggrepel_0.9.1               latticeExtra_0.6-29         fastmatch_1.1-0            
# [118] tools_4.0.3                 rstudioapi_0.13             statip_0.2.3               
# [121] gridExtra_2.3               farver_2.1.0                ggraph_2.0.5               
# [124] stable_1.1.4                BiocManager_1.30.12         rvcheck_0.1.8              
# [127] digest_0.6.27               pracma_2.3.3                Rcpp_1.0.6                 
# [130] GenomicRanges_1.42.0        broom_0.7.6                 fda_5.1.9                  
# [133] httr_1.4.2                  IDPmisc_1.1.20              colorspace_2.0-0           
# [136] rvest_1.0.0                 XML_3.99-0.6                fs_1.5.0                   
# [139] pdist_1.2                   rainbow_3.6                 modeest_2.4.0              
# [142] splines_4.0.3               rmutil_1.1.5                RBGL_1.66.0                
# [145] conquer_1.0.2               graphlayouts_0.7.1          xtable_1.8-4               
# [148] jsonlite_1.7.2              fds_1.8                     tidygraph_1.2.0            
# [151] corpcor_1.6.9               timeDate_3043.102           UpSetR_1.4.0               
# [154] testthat_3.0.2              R6_2.5.0                    pillar_1.6.0               
# [157] htmltools_0.5.1.1           glue_1.4.2                  fastmap_1.1.0              
# [160] BiocParallel_1.24.1         codetools_0.2-18            fgsea_1.16.0               
# [163] pcaPP_1.9-74                mvtnorm_1.1-1               utf8_1.2.1                 
# [166] curl_4.3                    gtools_3.8.2                GO.db_3.12.1               
# [169] openssl_1.4.3               survival_3.2-7              rmarkdown_2.7              
# [172] munsell_0.5.0               DO.db_2.9                   GetoptLong_1.0.5           
# [175] GenomeInfoDbData_1.2.4      iterators_1.0.13            haven_2.4.1                
# [178] reshape2_1.4.4              gtable_0.3.0        
# 
# 
