# look at some clinical data

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



i <- projects[1]
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
    
    
    index <- BLCA_met$malignancy_type == "Synchronous Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$malignancy_type == "Synchronous Malignancy|Prior Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
                                                         
    index <- BLCA_met$malignancy_type == "Prior Malignancy|Synchronous Malignancy" 
    BLCA_met$Metastatic_status[index] <- 1
    
    
    write.csv(BLCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_status_.csv"))
    
  
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
        tidyr::unite(neoplasm_and_distant_met, met_site_merge:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE,sep = ",")
        
        
        
        
        
        # mutate(BRCA_met, LymphNodeStatus = ifelse(is.na(number_of_lymphnodes_positive_by_he & number_of_lymphnodes_positive_by_ihc), 0,
        #                                           ifelse(number_of_lymphnodes_positive_by_he == 0 & number_of_lymphnodes_positive_by_ihc == 0, 0,
        #                                                  ifelse(number_of_lymphnodes_positive_by_he == 0 & is.na(number_of_lymphnodes_positive_by_ihc), 0,
        #                                                         ifelse(is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_ihc ==0, 0, 1)))))
        # 
      index <- is.na(BRCA_met$LymphNodeStatus)
      BRCA_met$LymphNodeStatus[index] <- 0
      
      #lymph node status
      
      
      
      
      
      
      # consolidate columns
      
      
      write.csv(BRCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
      
      print("BRCA Done")
      rm(BRCA_clinical)
      rm(BRCA_met)
    }
  
  if(i== "TCGA-COAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    COAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_neoplasm_event_type ,new_tumor_event_after_initial_treatment,other_malignancy_anatomic_site, other_malignancy_type,
                                              site_of_additional_surgery_new_tumor_event_mets,other_malignancy_laterality))
    
    write.csv(COAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("COAD Done")
    rm(COAD_met)
  }
  
  if(i== "TCGA-ESCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    ESCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,number_of_lymphnodes_positive_by_ihc,number_of_lymphnodes_positive_by_he,
                                              new_tumor_event_after_initial_treatment,new_neoplasm_event_type, new_tumor_event_additional_surgery_procedure,
                                              new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_event_occurrence_anatomic_site_text,other_malignancy_type,other_malignancy_anatomic_site))
    
    write.csv(ESCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("ESCA Done")
    rm(ESCA_met)
  } 
  
  if(i== "TCGA-HNSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    HNSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,number_of_lymphnodes_positive_by_ihc,number_of_lymphnodes_positive_by_he,
                                              new_tumor_event_after_initial_treatment, additional_surgery_metastatic_procedure,
                                              new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,
                                              new_neoplasm_event_type,new_tumor_event_additional_surgery_procedure,
                                              malignancy_type,other_malignancy_anatomic_site))
    
    
    write.csv(EHNSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("HNSC Done")
    rm(HNSC_met)
  } 
  
  if(i== "TCGA-KICH"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    KICH_met <- as.data.frame(dat %>% dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment, number_of_lymphnodes_positive,additional_radiation_therapy,
                                                       additional_pharmaceutical_therapy,
                                                       additional_surgery_metastatic_procedure,malignancy_type,
                                                       other_malignancy_anatomic_site,other_malignancy_histological_type,
                                                       other_malignancy_histological_type_text,other_malignancy_laterality))
    
    write.csv(KICH_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KICH is bulshit but it's Done")
    rm(KICH_met)
  } 
  
  if(i== "TCGA-KIRC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment, number_of_lymphnodes_positive,additional_radiation_therapy,
                                              additional_pharmaceutical_therapy,
                                              additional_surgery_metastatic_procedure,malignancy_type,
                                              other_malignancy_anatomic_site,other_malignancy_histological_type,
                                              other_malignancy_histological_type_text,other_malignancy_laterality))
    
    
    write.csv(KIRC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRC Done")
    rm(KIRC_met)
  } 
  
  if(i== "TCGA-KIRP"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRP_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment, number_of_lymphnodes_positive,additional_radiation_therapy,
                                              additional_pharmaceutical_therapy,
                                              additional_surgery_metastatic_procedure,malignancy_type,
                                              other_malignancy_anatomic_site,other_malignancy_histological_type,
                                              other_malignancy_histological_type_text,other_malignancy_laterality,other_malignancy_anatomic_site_text))
    
    
    write.csv(KIRP_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRP Done")
    rm(KIRP_met)
  
  
  
  } 
  
  if(i== "TCGA-LIHC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LIHC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment,new_neoplasm_event_type,
                                              new_neoplasm_event_occurrence_anatomic_site,new_tumor_event_additional_surgery_procedure,
                                              new_neoplasm_occurrence_anatomic_site_text))
    
    
    write.csv(LIHC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LIHC Done")
    rm(LIHC_met)
    
    
    
  } 
  
  
  if(i== "TCGA-LUAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LUAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment,
                                              new_neoplasm_event_type,
                                              additional_surgery_metastatic_procedure,
                                              additional_surgery_locoregional_procedure,
                                              other_malignancy_anatomic_site,location_in_lung_parenchyma,
                                              `new_neoplasm_event_type[1]`,`new_neoplasm_event_type[2]`,
                                              malignancy_type,other_malignancy_histological_type,other_malignancy_histological_type_text))

    write.csv(LUAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUAD Done")
    rm(LUAD_met)
    
    
    
  } 
  
  if(i== "TCGA-LUSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LUSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,new_tumor_event_after_initial_treatment,
                                              new_neoplasm_event_type,
                                              additional_surgery_metastatic_procedure,
                                              additional_surgery_locoregional_procedure,
                                              other_malignancy_anatomic_site,location_in_lung_parenchyma,
                                              `new_neoplasm_event_type[1]`,`new_neoplasm_event_type[2]`,
                                              malignancy_type,other_malignancy_histological_type,other_malignancy_histological_type_text))
    
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




