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



i <- projects[8]
for(i in projects){
  dat <- as.data.frame(datatable::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv")))
  
  if(i== "TCGA-BLCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    BLCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive_by_he, number_of_lymphnodes_positive_by_ihc,
                                              metastatic_site_at_diagnosis, new_neoplasm_event_occurrence_anatomic_site,
                                              metastatic_site_at_diagnosis_other,other_malignancy_laterality,
                                              `metastatic_site_at_diagnosis[1]`,`metastatic_site_at_diagnosis[2]`,`metastatic_site_at_diagnosis[3]`,
                                              `metastatic_site_at_diagnosis[4]`,other_malignancy_anatomic_site_text))
    #lymph node status
    
    # consolidate columns
    
    
    
    
    write.csv(BLCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
  
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










