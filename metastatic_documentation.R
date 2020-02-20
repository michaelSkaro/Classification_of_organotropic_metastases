# Clear your workspace
rm(list=ls(all=TRUE))
# Load some things
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(gsubfn)
library(dplyr)

# make a working directory to read all the files from
setwd("~/CSBL_shared/clinical/TCGA_xml") # most comprehensive


# lets decide which data sets we are going to work on
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

#projects <- c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")
# annotate weird gene symbols and get some clinical data
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")
refDat <- data.table::fread("~/storage/Metastatic_Organo_Tropism/Metastatic_database_project_information.csv")
refDat <- refDat[order(refDat$Sample_id),]
refDat<- refDat[match(unique(refDat$Sample_id), refDat$Sample_id),]
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
clinical$Sample_id <- substr(clinical$barcode, 0,16)





## define a helper function
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}



i <- projects[1]

# do the things

for(i in projects){
  
  if(i== "TCGA-BLCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    BLCA_met <- as.data.frame(dat %>% dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site, metastatic_site,`metastatic_site[1]`,
                                                    `metastatic_site[2]`,`metastatic_site[3]`,new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,
                                                    number_of_lymphnodes_positive_by_he)) %>%
      mutate(BLCA_met, LymphNodeStatus = ifelse(is.na(number_of_lymphnodes_positive_by_he), 0,
                                                ifelse(number_of_lymphnodes_positive_by_he == 0, 0, 1))) %>%
      tidyr::unite(met_loc, other_malignancy_anatomic_site:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE, sep = ",") %>%
      mutate_each(funs(empty_as_na))
      
      
      index <- BLCA_met$met_loc == ",,,,,," 
      BLCA_met$met_loc[index] <- NA
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,,,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      
      index <- BLCA_met$met_loc == ",None," 
      BLCA_met$met_loc[index] <- NA
      BLCA_met$met_loc <- stringr::str_replace(BLCA_met$met_loc, ",,", "")
      BLCA_met$met_loc <- stringr::str_replace(BLCA_met$met_loc, ",", "")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Other specify", "")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Other, specify", "")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "[|]", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "None,", "")
      
      index <- BLCA_met$malignancy_type == "Prior Malignancy|Synchronous Malignancy" 
      BLCA_met$Metastatic_status[index] <- 1
    
      
      mutate(BLCA_met, Metastatic_status = case_when(
        is.na(new_neoplasm_event_occurrence_anatomic_site) & is.na(new_neoplasm_event_type) & 
          is.na(new_neoplasm_occurrence_anatomic_site_text) | malignancy_type == " Prior Malignancy" ~ 0,
        TRUE ~1)) %>%
      
      
    BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "[||]", ",")
    
    index <- BLCA_met$met_loc == ",No New Tumor Event"
    BLCA_met$Metastatic_status[index] <- 0
    
    
    
    index <- BLCA_met$malignancy_type == "Synchronous Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$malignancy_type == "Synchronous Malignancy|Prior Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
                                                         
    index <- BLCA_met$malignancy_type == "Prior Malignancy|Synchronous Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
    
    
    write.csv(BLCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_status.csv"))
    
  
    print("BLCA Done")
    rm(BLCA_met)
    rm(dat)
  }
  
    if(i = "TCGA-BRCA"){
      
      dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv")))
      
      ## define a helper function
      empty_as_na <- function(x){
        if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
        ifelse(as.character(x)!="", x, NA)
      }
      
      ## transform all columns
      
      
      
      BRCA_met <- as.data.frame(dat %>%
                                  dplyr::select(bcr_patient_barcode, metastatic_site_at_diagnosis, number_of_lymphnodes_positive_by_he,
                                               `metastatic_site_at_diagnosis[1]`,`metastatic_site_at_diagnosis[2]`,
                                               `metastatic_site_at_diagnosis[3]`,`metastatic_site_at_diagnosis[4]`,
                                                new_neoplasm_event_occurrence_anatomic_site, 
                                                new_neoplasm_occurrence_anatomic_site_text))%>% 
        mutate_each(funs(empty_as_na)) %>%
        mutate(BRCA_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
        tidyr::unite(met_site_merge, `metastatic_site_at_diagnosis[1]`:`metastatic_site_at_diagnosis[4]`, na.rm =TRUE, sep = ",") %>%
        tidyr::unite(neoplasm_and_distant_met, met_site_merge:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
        dplyr::select(bcr_patient_barcode,metastatic_site_at_diagnosis,neoplasm_and_distant_met,number_of_lymphnodes_positive_by_he,LymphNodeStatus) %>%
        tidyr::unite(met_loc, metastatic_site_at_diagnosis:neoplasm_and_distant_met,na.rm =TRUE,sep = ",") %>%
        mutate_each(funs(empty_as_na))%>%
        mutate(BRCA_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0, 0, 1))
        
      index <- is.na(BRCA_met$LymphNodeStatus)
      BRCA_met$LymphNodeStatus[index] <- 0
      
      #lymph node status
      
      index <- is.na(BRCA_met$number_of_lymphnodes_positive_by_he)
      BRCA_met$number_of_lymphnodes_positive_by_he[index] <- 0
      
      
      
      # consolidate columns
      
      
      write.csv(BRCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
      
      print("BRCA Done")
      rm(BRCA_met)
    }
  
  if(i== "TCGA-COAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    COAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,new_neoplasm_event_type, other_malignancy_anatomic_site, 
                                              site_of_additional_surgery_new_tumor_event_mets,number_of_lymphnodes_positive_by_he)) %>%
                                tidyr::unite(met_loc,new_neoplasm_event_type:site_of_additional_surgery_new_tumor_event_mets, na.rm =TRUE,sep = ",") %>%
                                mutate_each(funs(empty_as_na)) %>%
                                mutate(COAD_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                mutate(COAD_met, Metastatic_status = ifelse(met_loc ==",," & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 
                                
      
                            
                              
    index <- is.na(COAD_met$LymphNodeStatus)
    COAD_met$LymphNodeStatus[index] <- 0
    
    
    index <- COAD_met$malignancy_type == "Prior Malignancy"
    COAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(COAD_met$Metastatic_status)
    COAD_met$Metastatic_status[index] <- 0
    
    
    index <- COAD_met$met_loc == ",,"
    COAD_met$met_loc[index] <- NA
    
    index <- is.na(COAD_met$number_of_lymphnodes_positive_by_he)
    COAD_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    write.csv(COAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("COAD Done")
    rm(COAD_met)
  }
  
  if(i== "TCGA-ESCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    ESCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive_by_he,new_neoplasm_event_type,
                                              new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,other_malignancy_anatomic_site)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                tidyr::unite(met_loc,new_neoplasm_event_type:other_malignancy_anatomic_site, na.rm =TRUE,sep = ",") %>%
                                mutate(ESCA_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(ESCA_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 
      
    index <- is.na(ESCA_met$LymphNodeStatus)
    ESCA_met$LymphNodeStatus[index] <- 0
    
    
    index <- ESCA_met$malignancy_type == "Prior Malignancy"
    ESCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(ESCA_met$Metastatic_status)
    ESCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(ESCA_met$number_of_lymphnodes_positive_by_he)
    ESCA_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
                                            
    
    write.csv(ESCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("ESCA Done")
    rm(ESCA_met)
  } 
  
  if(i== "TCGA-HNSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    HNSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive_by_he,new_neoplasm_event_type, other_malignancy_anatomic_site_text,
                                      new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,other_malignancy_anatomic_site)) %>%
                                      mutate_each(funs(empty_as_na))%>%
                                      tidyr::unite(met_loc,new_neoplasm_event_type:other_malignancy_anatomic_site, na.rm =TRUE,sep = ",") %>%
                                      mutate(HNSC_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                      mutate_each(funs(empty_as_na))%>%
                                      mutate(HNSC_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 

    index <- is.na(HNSC_met$LymphNodeStatus)
    HNSC_met$LymphNodeStatus[index] <- 0
    
    
    index <- HNSC_met$malignancy_type == "Prior Malignancy"
    HNSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(HNSC_met$Metastatic_status)
    HNSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(HNSC_met$number_of_lymphnodes_positive_by_he)
    HNSC_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    
    write.csv(HNSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("HNSC Done")
    rm(HNSC_met)
  } 
  
  if(i== "TCGA-KICH"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    KICH_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive, 
                                              other_malignancy_anatomic_site)) %>%
                                mutate(KICH_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(KICH_met, met_loc = ifelse(malignancy_type == "Synchronous Malignancy", malignancy_type, NA)) %>%
                                mutate(KICH_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0, 0, 1)) 
    
    write.csv(KICH_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KICH is bulshit but it's Done")
    rm(KICH_met)
  } 
  
  if(i== "TCGA-KIRC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive,malignancy_type,
                                              other_malignancy_anatomic_site, other_malignancy_anatomic_site_text)) %>% 
                              mutate_each(funs(empty_as_na))%>%
                              tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
                              mutate(KIRC_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
                              mutate(KIRC_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus >0 | is.na(LymphNodeStatus) | malignancy_type =="Prior Malignnacy" , 0, 1))
    
    index <- is.na(KIRC_met$LymphNodeStatus)
    KIRC_met$LymphNodeStatus[index] <- 0
    
    index <- KIRC_met$LymphNodeStatus >0
    KIRC_met$Metastatic_status[index] <- 1
    
    index <- KIRC_met$malignancy_type == "Prior Malignancy"
    KIRC_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRC_met$Metastatic_status)
    KIRC_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRC_met$number_of_lymphnodes)
    KIRC_met$number_of_lymphnodes_positive[index] <- 0
    
    write.csv(KIRC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRC Done")
    rm(KIRC_met)
  } 
  
  if(i== "TCGA-KIRP"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRP_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive,malignancy_type,
                                              other_malignancy_anatomic_site, other_malignancy_anatomic_site_text)) %>% 
      mutate_each(funs(empty_as_na))%>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate(KIRP_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
      mutate(KIRP_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus >0 | is.na(LymphNodeStatus) | malignancy_type =="Prior Malignnacy" , 0, 1))
    
    index <- is.na(KIRP_met$LymphNodeStatus)
    KIRP_met$LymphNodeStatus[index] <- 0
    
    index <- KIRP_met$LymphNodeStatus >0
    KIRP_met$Metastatic_status[index] <- 1
    
    index <- KIRP_met$malignancy_type == "Prior Malignancy"
    KIRP_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRP_met$Metastatic_status)
    KIRP_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRP_met$number_of_lymphnodes)
    KIRP_met$number_of_lymphnodes_positive[index] <- 0
    
    write.csv(KIRP_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRP Done")
    rm(KIRP_met)
  
  
  
  } 
  
  if(i== "TCGA-LIHC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LIHC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, malignancy_type, other_malignancy_anatomic_site, other_malignancy_anatomic_site_text,new_neoplasm_event_occurrence_anatomic_site)) %>% 
      tidyr::unite(met_loc,other_malignancy_anatomic_site:new_neoplasm_event_occurrence_anatomic_site, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(LIHC_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" , 1, 0))
    
    LIHC_met[LIHC_met == ",,"] <- NA
    
    
    index <- !is.na(LIHC_met$met_loc) & is.na(LIHC_met$Metastatic_status)
    LIHC_met$Metastatic_status[index] <- 1
    
    index <- LIHC_met$met_loc
    LIHC_met$Metastatic_status[index] <- 1
    
    index <- LIHC_met$malignancy_type == "Prior Malignancy"
    LIHC_met$Metastatic_status[index] <- 0
    
    index <- is.na(LIHC_met$Metastatic_status)
    LIHC_met$Metastatic_status[index] <- 0
    
    LIHC_met$met_loc <- str_replace_all(LIHC_met$met_loc, ",,","")
    
    
    write.csv(LIHC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LIHC Done")
    rm(LIHC_met)
    
    
    
  } 
  
  
  if(i== "TCGA-LUAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 

    LUAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,
                                              other_malignancy_anatomic_site, pathologic_N,pathologic_M))%>% 
      tidyr::unite(met_loc,other_malignancy_anatomic_site, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(LUAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(LUAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    LUAD_met[LUAD_met == ",,"] <- NA
    
    LUAD_met[LUAD_met == "MX"] <- "M0"
    
    index <- is.na(LUAD_met$pathologic_M)
    LUAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(LUAD_met$pathologic_N)
    LUAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(LUAD_met$met_loc) & is.na(LUAD_met$Metastatic_status)
    LUAD_met$Metastatic_status[index] <- 0
    
    index <- LUAD_met$malignancy_type == "Prior Malignancy"
    LUAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(LUAD_met$Metastatic_status)
    LUAD_met$Metastatic_status[index] <- 0
    
    LUAD_met$met_loc <- str_replace_all(LUAD_met$met_loc, ",,","")
    
    
  
    write.csv(LUAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUAD Done")
    rm(LUAD_met)
    
    
    
  } 
  
  if(i== "TCGA-LUSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LUSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
                                tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(LUSC_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                                            | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
                                mutate(LUSC_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    LUSC_met[LUSC_met == ","] <- NA
    
    LUSC_met[LUSC_met == "MX"] <- "M0"
    
    index <- is.na(LUSC_met$pathologic_M)
    LUSC_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(LUSC_met$pathologic_N)
    LUSC_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(LUSC_met$met_loc) & is.na(LUSC_met$Metastatic_status)
    LUSC_met$Metastatic_status[index] <- 0
    
    index <- LUSC_met$malignancy_type == "Prior Malignancy"
    LUSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(LUSC_met$Metastatic_status)
    LUSC_met$Metastatic_status[index] <- 0
    
    LUSC_met$met_loc <- str_replace_all(LUSC_met$met_loc, ",,","")
    
    write.csv(LUSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUSC Done")
    rm(LUSC_met)
    
    
    
  } 
  
  if(i== "TCGA-PRAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    PRAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(PRAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(PRAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    PRAD_met[PRAD_met == ","] <- NA
    
    PRAD_met[PRAD_met == "MX"] <- "M0"
    
    index <- is.na(PRAD_met$pathologic_M)
    PRAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(PRAD_met$pathologic_N)
    PRAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(PRAD_met$met_loc) & is.na(PRAD_met$Metastatic_status)
    PRAD_met$Metastatic_status[index] <- 0
    
    index <- PRAD_met$malignancy_type == "Prior Malignancy"
    PRAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(PRAD_met$Metastatic_status)
    PRAD_met$Metastatic_status[index] <- 0
    
    PRAD_met$met_loc <- str_replace_all(PRAD_met$met_loc, ",,","")
    
    
    index <- PRAD_met$LymphNodeStatus > 0
    PRAD_met$Metastatic_status[index] <- 1
    
    
    
    write.csv(PRAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("PRAD Done")
    rm(PRAD_met)
    
  } 
  
  if(i== "TCGA-STAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    STAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(STAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(STAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    STAD_met[STAD_met == ","] <- NA
    
    STAD_met[STAD_met == "MX"] <- "M0"
    
    index <- is.na(STAD_met$pathologic_M)
    STAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(STAD_met$pathologic_N)
    STAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(STAD_met$met_loc) & is.na(STAD_met$Metastatic_status)
    STAD_met$Metastatic_status[index] <- 0
    
    index <- STAD_met$malignancy_type == "Prior Malignancy"
    STAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(STAD_met$Metastatic_status)
    STAD_met$Metastatic_status[index] <- 0
    
    STAD_met$met_loc <- str_replace_all(STAD_met$met_loc, ",,","")
    
    
    index <- STAD_met$LymphNodeStatus > 0
    STAD_met$Metastatic_status[index] <- 1
    
    
 
    
    write.csv(STAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("STAD Done")
    rm(STAD_met)
    
    
    
  } 
  
  if(i== "TCGA-THCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    THCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(THCA_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(THCA_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    THCA_met[THCA_met == ","] <- NA
    
    THCA_met[THCA_met == "MX"] <- "M0"
    
    index <- is.na(THCA_met$pathologic_M)
    THCA_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(THCA_met$pathologic_N)
    THCA_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(THCA_met$met_loc) & is.na(THCA_met$Metastatic_status)
    THCA_met$Metastatic_status[index] <- 0
    
    index <- THCA_met$malignancy_type == "Prior Malignancy"
    THCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(THCA_met$Metastatic_status)
    THCA_met$Metastatic_status[index] <- 0
    
    THCA_met$met_loc <- str_replace_all(THCA_met$met_loc, ",,","")
    
    
    index <- THCA_met$LymphNodeStatus > 0
    THCA_met$Metastatic_status[index] <- 1
    
    write.csv(THCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("THCA Done")
    rm(THCA_met)
    
    
    
  } 
  
}

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pbapply_1.4-2                            EBImage_4.24.0                           remotes_2.1.0                            gsubfn_0.7                              
# [5] proto_1.0.0                              forcats_0.4.0                            purrr_0.3.3                              readr_1.3.1                             
# [9] tidyr_1.0.2                              tibble_2.1.3                             tidyverse_1.3.0                          targetscan.Mm.eg.db_0.6.1               
# [13] data.table_1.12.8                        RColorBrewer_1.1-2                       gplots_3.0.1.2                           Glimma_1.10.1                           
# [17] edgeR_3.24.3                             limma_3.38.3                             vsn_3.50.0                               pheatmap_1.0.12                         
# [21] apeglm_1.4.2                             EnsDb.Mmusculus.v79_2.99.0               ensembldb_2.6.8                          AnnotationFilter_1.6.0                  
# [25] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4 GenomicFeatures_1.34.8                   tximport_1.10.1                          org.Mm.eg.db_3.7.0                      
# [29] rpart_4.1-13                             randomForest_4.6-14                      keras_2.2.5.0                            kernlab_0.9-29                          
# [33] caret_6.0-85                             lattice_0.20-38                          plyr_1.8.5                               stringr_1.4.0                           
# [37] dplyr_0.8.4                              WGCNA_1.68                               fastcluster_1.1.25                       dynamicTreeCut_1.63-1                   
# [41] org.Hs.eg.db_3.7.0                       AnnotationDbi_1.44.0                     clusterProfiler_3.10.1                   EnhancedVolcano_1.0.1                   
# [45] ggrepel_0.8.1                            ggplot2_3.2.1                            DESeq2_1.22.2                            SummarizedExperiment_1.12.0             
# [49] DelayedArray_0.8.0                       BiocParallel_1.16.6                      matrixStats_0.55.0                       Biobase_2.42.0                          
# [53] GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                      IRanges_2.16.0                           S4Vectors_0.20.1                        
# [57] BiocGenerics_0.28.0                     
# 
# loaded via a namespace (and not attached):
# [1] rtracklayer_1.42.2       ModelMetrics_1.2.2.1     coda_0.19-3              acepack_1.4.1            bit64_0.9-7              knitr_1.27               RCurl_1.98-1.1           doParallel_1.0.15       
# [9] generics_0.0.2           preprocessCore_1.44.0    cowplot_1.0.0            RSQLite_2.2.0            europepmc_0.3            tensorflow_2.0.0         bit_1.1-15.1             enrichplot_1.2.0        
# [17] xml2_1.2.2               lubridate_1.7.4          assertthat_0.2.1         viridis_0.5.1            gower_0.2.1              xfun_0.12                hms_0.5.3                DEoptimR_1.0-8          
# [25] fansi_0.4.1              progress_1.2.2           readxl_1.3.1             dbplyr_1.4.2             caTools_1.17.1.3         igraph_1.2.4.2           DBI_1.1.0                geneplotter_1.60.0      
# [33] htmlwidgets_1.5.1        ellipsis_0.3.0           backports_1.1.5          annotate_1.60.1          biomaRt_2.38.0           vctrs_0.2.2              abind_1.4-5              withr_2.1.2             
# [41] ggforce_0.3.1            triebeard_0.3.0          robustbase_0.93-5        bdsmatrix_1.3-4          checkmate_1.9.4          GenomicAlignments_1.18.1 prettyunits_1.1.1        cluster_2.0.7-1         
# [49] DOSE_3.8.2               lazyeval_0.2.2           crayon_1.3.4             genefilter_1.64.0        recipes_0.1.9            pkgconfig_2.0.3          tweenr_1.0.1             nlme_3.1-137            
# [57] ProtGenerics_1.14.0      nnet_7.3-12              rlang_0.4.4              lifecycle_0.1.0          affyio_1.52.0            modelr_0.1.5             cellranger_1.1.0         polyclip_1.10-0         
# [65] tiff_0.1-5               Matrix_1.2-15            urltools_1.7.3           reprex_0.3.0             base64enc_0.1-3          whisker_0.4              ggridges_0.5.2           png_0.1-7               
# [73] viridisLite_0.3.0        bitops_1.0-6             KernSmooth_2.23-15       pROC_1.16.1              Biostrings_2.50.2        blob_1.2.1               qvalue_2.14.1            robust_0.4-18.2         
# [81] jpeg_0.1-8.1             gridGraphics_0.4-1       scales_1.1.0             memoise_1.1.0            magrittr_1.5             gdata_2.18.0             zlibbioc_1.28.0          compiler_3.5.1          
# [89] bbmle_1.0.23.1           rrcov_1.5-2              Rsamtools_1.34.1         cli_2.0.1                affy_1.60.0              XVector_0.22.0           htmlTable_1.13.3         Formula_1.2-3           
# [97] MASS_7.3-51.1            tidyselect_1.0.0         stringi_1.4.5            yaml_2.2.1               emdbook_1.3.11           GOSemSim_2.8.0           locfit_1.5-9.1           latticeExtra_0.6-28     
# [105] grid_3.5.1               fastmatch_1.1-0          tools_3.5.1              rstudioapi_0.10          foreach_1.4.7            foreign_0.8-71           gridExtra_2.3            prodlim_2019.11.13      
# [113] farver_2.0.3             ggraph_2.0.0             digest_0.6.23            rvcheck_0.1.7            BiocManager_1.30.10      lava_1.6.6               Rcpp_1.0.3               broom_0.5.4             
# [121] httr_1.4.1               colorspace_1.4-1         rvest_0.3.5              fs_1.3.1                 XML_3.99-0.3             reticulate_1.14          splines_3.5.1            graphlayouts_0.5.0      
# [129] ggplotify_0.0.4          fit.models_0.5-14        xtable_1.8-4             jsonlite_1.6.1           tidygraph_1.1.2          timeDate_3043.102        UpSetR_1.4.0             zeallot_0.1.0           
# [137] ipred_0.9-9              R6_2.4.1                 Hmisc_4.3-0              pillar_1.4.3             htmltools_0.4.0          glue_1.3.1               fftwtools_0.9-8          class_7.3-14            
# [145] codetools_0.2-15         fgsea_1.8.0              pcaPP_1.9-73             mvtnorm_1.0-12           numDeriv_2016.8-1.1      curl_4.3                 tfruns_1.4               gtools_3.8.1            
# [153] GO.db_3.7.0              survival_3.1-8           munsell_0.5.0            DO.db_2.9                GenomeInfoDbData_1.2.0   iterators_1.0.12         impute_1.56.0            haven_2.2.0             
# [161] reshape2_1.4.3           gtable_0.3.
# 
