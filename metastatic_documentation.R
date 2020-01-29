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



i <- projects[5]
for(i in projects){
  dat <- as.data.frame(datatable::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv")))
  
  if(i== "TCGA-BLCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    BRCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive_by_he, number_of_lymphnodes_positive_by_ihc,
                                              metastatic_site_at_diagnosis, new_neoplasm_event_occurrence_anatomic_site,
                                              metastatic_site_at_diagnosis_other,other_malignancy_laterality,
                                              `metastatic_site_at_diagnosis[1]`,`metastatic_site_at_diagnosis[2]`,`metastatic_site_at_diagnosis[3]`,
                                              `metastatic_site_at_diagnosis[4]`,other_malignancy_anatomic_site_text))
    #lymph node status
    
    # consolidate columns
    
    
    
    
    write.csv(dat, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    dat$`metastatic_site_at_diagnosis[1]`
    
    
    print("BLCA Done")
    rm(BLCA_met)
    rm(dat)
  }
  
    if(i = "TCGA-BRCA"){
      
      dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
      
      BRCA_met <- as.data.frame(dat %>%
                                  dplyr::select(bcr_patient_barcode, metastatic_site, number_of_lymphnodes_positive_by_he, number_of_lymphnodes_positive_by_ihc,
                                                `metastatic_site[1]`,`metastatic_site[2]`,`metastatic_site[3]`,new_neoplasm_event_occurrence_anatomic_site))
      
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
    
    write.csv(ESCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("HNSC Done")
    rm(HNSC_met)
  } 
  
  
  
  
  
}




