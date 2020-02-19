# look at some clinical data
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(gsubfn)

setwd("~/CSBL_shared/clinical/TCGA_xml") # most comprehensive



projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

#projects <- c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")

annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")
refDat <- data.table::fread("~/storage/Metastatic_Organo_Tropism/Metastatic_database_project_information.csv")
refDat <- refDat[order(refDat$Sample_id),]
refDat<- refDat[match(unique(refDat$Sample_id), refDat$Sample_id),]
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
clinical$Sample_id <- substr(clinical$barcode, 0,16)



i <- projects[8]
for(i in projects){
  
  if(i== "TCGA-BLCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    BLCA_met <- as.data.frame(dat %>% dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive_by_he,
                                                    other_malignancy_anatomic_site,metastatic_site,new_tumor_event_after_initial_treatment,
                                                    new_neoplasm_event_type,new_neoplasm_event_occurrence_anatomic_site,
                                                    new_neoplasm_occurrence_anatomic_site_text,new_tumor_event_additional_surgery_procedure,
                                                    metastatic_site,`metastatic_site[1]`,`metastatic_site[2]`,`metastatic_site[3]`))
    #lymph node status
    #View(BLCA_met)
    # consolidate columns
    library(dplyr)
    
    ## define a helper function
    empty_as_na <- function(x){
      if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
      ifelse(as.character(x)!="", x, NA)
    }
    
    ## transform all columns
    
    BLCA_met<- BLCA_met %>% 
      mutate_each(funs(empty_as_na)) %>%
      mutate(BLCA_met, LymphNodeStatus = ifelse(is.na(number_of_lymphnodes_positive_by_he), 0,
                                                ifelse(number_of_lymphnodes_positive_by_he == 0, 0, 1))) %>%
      mutate(BLCA_met, Metastatic_status = case_when(
        is.na(new_neoplasm_event_occurrence_anatomic_site) & is.na(new_neoplasm_event_type) & 
          is.na(new_neoplasm_occurrence_anatomic_site_text) | malignancy_type == " Prior Malignancy" ~ 0,
        TRUE ~1)) %>%
      tidyr::unite(metastatic_site, `metastatic_site[1]`:`metastatic_site[3]`, na.rm =TRUE, sep = "|") %>%
      tidyr::unite(met_site_merge, new_neoplasm_event_type:new_neoplasm_occurrence_anatomic_site_text , na.rm =TRUE, sep = "|") %>%
      dplyr::select(bcr_patient_barcode, malignancy_type, metastatic_site, met_site_merge, other_malignancy_anatomic_site, 
                    number_of_lymphnodes_positive_by_he, LymphNodeStatus, Metastatic_status) %>%
      tidyr::unite(met_loc, metastatic_site:other_malignancy_anatomic_site , na.rm =TRUE, sep = "|") 
      
      
    
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
    LUAD_met$pathologic_M[index] <- "N0"
    
    
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
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment,
                                              new_neoplasm_event_type,
                                              other_malignancy_anatomic_site,
                                              `new_neoplasm_event_type[1]`,`new_neoplasm_event_type[2]`))
    
    write.csv(LUSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUSC Done")
    rm(LUSC_met)
    
    
    
  } 
  
  if(i== "TCGA-PRAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    PRAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment,new_tumor_event_after_initial_treatment,
                                              bone_scan_results,new_neoplasm_event_type,new_neoplasm_event_occurrence_anatomic_site,malignancy_type,
                                              other_malignancy_anatomic_site,histological_type_other,tumor_progression_post_ht,other_malignancy_anatomic_site_text,
                                              other_malignancy_histological_type, other_malignancy_histological_type_text,new_neoplasm_occurrence_anatomic_site_text
                                              ))
    
    write.csv(PRAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("PRAD Done")
    rm(PRAD_met)
    
    
    
  } 
  
  if(i== "TCGA-STAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    STAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,number_of_lymphnodes_positive_by_he,
                                              new_tumor_event_after_initial_treatment,new_neoplasm_event_type,
                                              new_neoplasm_event_occurrence_anatomic_site,additional_surgery_metastatic_procedure,
                                              malignancy_type,other_malignancy_anatomic_site,other_malignancy_histological_type,
                                              other_malignancy_histological_type_text,additional_surgery_locoregional_procedure,
                                              new_neoplasm_occurrence_anatomic_site_text,days_to_additional_surgery_locoregional_procedure,
                                              residual_disease_post_new_tumor_event_margin_status,`new_neoplasm_event_type[1]`,`new_neoplasm_event_type[2]`,
                                              other_malignancy_anatomic_site_text))
    
    write.csv(PRAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("STAD Done")
    rm(STAD_met)
    
    
    
  } 
  
  if(i== "TCGA-THCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    THCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,number_of_lymphnodes_positive_by_he,lymph_node_preoperative_assessment_diagnostic_imaging_type,
                                              new_tumor_event_after_initial_treatment,new_neoplasm_event_type,
                                              new_neoplasm_event_occurrence_anatomic_site,metastatic_site, other_metastatic_site,
                                              malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_histological_type,new_neoplasm_event_type,new_tumor_event_additional_surgery_procedure,
                                              new_neoplasm_event_occurrence_anatomic_site,other_malignancy_histological_type_text,metastatic_neoplasm_confirmed_diagnosis_method_name,
                                              `metastatic_neoplasm_confirmed_diagnosis_method_name[1]`,`metastatic_neoplasm_confirmed_diagnosis_method_name[2]`))
    
    write.csv(PRAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("THCA Done")
    rm(THCA_met)
    
    
    
  } 
  
  
  
  
  
  
  
}

