# working with curated TCGA data
setwd("~/storage/MAE_analysis")

library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)


projects <- c("BLCA","BRCA","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")

proj <- projects[1]
for(proj in projects){

    suppressMessages({
      dat <- curatedTCGAData(proj,
                             assays = c("miRNASeqGene", "Mutation", "RNASeq2GeneNorm", "Methylation"),
                             dry.run = FALSE)
    })
    
    # describe the data we have downloaded
    
    colData(dat)[1:4, 1:10]
    
    # return a list of the data types being investigated
    
    experiments(dat)
    
    # make sample map, we neep to map all of the objects and their values back to patient, their stage and progression info
    
    samples_map <- sampleMap(dat)
    write.csv(samples_map, file = str_glue("~/storage/MAE_analysis/Sample_Maps/{proj}_sample_map.csv"))
    
    
    # subsetting will be important for experiments where we need to merge the patients in each one of the classes 
    
    # multiassayexperiment[i = rownames, j = primary or colnames, k = assay]
      # i = Subsetting operations always return another MultiAssayExperiment. 
      # For example, the following will return any rows named “MAPK14” or “IGFBP2”, and remove any assays where no rows match:
      # dat[c("MAPK14", "IGFBP2"), , ]
    
      # j= The following will keep only patients of pathological stage iv, and all their associated assays:
      # dat[, miniACC$pathologic_stage == "stage iv", ]
    
      # k = And the following will keep only the RNA-seq dataset, and only patients for which this assay is available:
      # dat[, , "RNASeq2GeneNorm"]
    
    # export the data as matrix csv files
    
    clin <- colData(dat)   
    save(clin, file = str_glue("~/storage/MAE_analysis/clinical/{proj}_clin.RData"))
    
    # methylation
    
    dat.exp <- assays(dat)[[1]]
    
    write.csv(dat.exp, file = str_glue("~/storage/MAE_analysis/methylation/{proj}_methly.csv"))
    
    
    # miRNa
    dat.exp <- assays(dat)[[2]]
    
    write.csv(dat.exp, file = str_glue("~/storage/MAE_analysis/miRNA/{proj}_miRNA.csv"))
    
    
    # mutation
    
    dat.exp <- assays(dat)[[3]]
    
    write.csv(dat.exp, file = str_glue("~/storage/MAE_analysis/mutation/{proj}_mut.csv"))
    
    # expression
    
    dat.exp <- assays(dat)[[4]]
    
    write.csv(dat.exp, file = str_glue("~/storage/MAE_analysis/expression/{proj}_expr.csv"))
    
    # save raw objects if we need to return them
    
    save(dat, file = str_glue("~/storage/MAE_analysis/MAE_objects/{proj}_MAE.Rdata"))
    
  
  
}






projects <- c("BLCA","BRCA","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")

proj <- projects[1]
for(proj in projects){
  # prepare the data for the analysis of each data type:
  
  load(str_glue("~/storage/MAE_analysis/clinical/{proj}_clin.RData"))
  
  #View(as.data.frame(clin@listData))
  
  # expand to long format with the Primary tumors or the Helathy tissue Normals
  
  annot <- as.data.frame(clin@listData)
  
  # attach the labels for each of the following:
  
  # methylation
    # Cancer/Normal
    # Cancer type
    # pathologic TNM stage
    # site of progression if there is one
    # site of progression tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
    # map the samples to the UUID and map the progression to the sample
  
  df.exp <- data.table::fread(file = str_glue("/home/mskaro1/storage/MAE_analysis/methylation/{proj}_methly.csv"), header = TRUE) %>%
    as_tibble() %>% 
    column_to_rownames("V1")
  
  df.exp <- rownames_to_column(df.exp, var = "probe_ids")
  
  # convert colnames to a vector
  
  cols <- colnames(df.exp[,2:length(colnames(df.exp))])
  
  cols.translate <- barcodeToUUID(cols)
  
  # translate them to the caseIDs
  
  colnames(cols.translate) <- c("submitter_id", "caseID")
  
  df.grid<- left_join(cols.translate, tumor.samples, by=  "caseID")
  
  # the attach the case IDs and the clincal matches
  
  
  
    
}


#complete data annotation

setwd("~/storage/MAE_analysis/Sample_Maps")

projects <- c("BLCA","BRCA","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")


BLCA <- as_tibble(data.table::fread("BLCA_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
BRCA <- as_tibble(data.table::fread("BRCA_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
COAD <- as_tibble(data.table::fread("COAD_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
ESCA <- as_tibble(data.table::fread("ESCA_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
HNSC <- as_tibble(data.table::fread("HNSC_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
KIRC <- as_tibble(data.table::fread("KIRC_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
KIRP <- as_tibble(data.table::fread("KIRP_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
LIHC <- as_tibble(data.table::fread("LIHC_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
LUAD <- as_tibble(data.table::fread("LUAD_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
LUSC <- as_tibble(data.table::fread("LUSC_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
PRAD <- as_tibble(data.table::fread("PRAD_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
STAD <- as_tibble(data.table::fread("STAD_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
THCA <- as_tibble(data.table::fread("THCA_sample_map.csv", header = TRUE)) %>%
  column_to_rownames("V1")
dat <- as.data.frame(rbind(BLCA,BRCA,COAD,ESCA,HNSC,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,STAD,THCA))

write.csv(dat, "complete_sample_map.csv")
rm(BLCA)
rm(BRCA)
rm(COAD)
rm(ESCA)
rm(HNSC)
rm(KIRC)
rm(KIRP)
rm(LIHC)
rm(LUAD)
rm(LUSC)
rm(PRAD)
rm(STAD)
rm(THCA)

map_ids <- barcodeToUUID(dat$primary)
colnames(dat)[2] <- "submitter_id"

full_map <- left_join(dat, map_ids, by = "submitter_id")



projects <- c("BLCA","BRCA","COAD","ESCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")

proj <- projects[11]

for(proj in projects){
  
  # attach the labels for each of the following:
  
  # Cancer/Normal
  # Cancer type
  # pathologic TNM stage
  # site of progression if there is one
  # site of progression tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
  # map the samples to the UUID and map the progression to the sample
  
  
  # expression
  
  # Cancer/Normal
  # Cancer type
  # pathologic TNM stage
  # site of progression if there is one
  # site of progression tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
  # map the samples to the UUID and map the progression to the sample
  
  # miRNA
  # Cancer/Normal
  # Cancer type
  # pathologic TNM stage
  # site of progression if there is one
  # site of progression tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
  # map the samples to the UUID and map the progression to the sample
  
  # mutation
  # Cancer/Normal
  # Cancer type
  # pathologic TNM stage
  # site of progression if there is one
  # site of progression tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
  # map the samples to the UUID and map the progression to the sample
  
  load(str_glue("~/storage/MAE_analysis/clinical/{proj}_clin.RData"))
  
  clin <- as.data.frame(clin)
  
  patient <- clin %>%
    dplyr::select(c("patientID","pathology_T_stage","pathology_N_stage","pathology_M_stage","admin.file_uuid","patient.bcr_patient_uuid"))
  
  colnames(patient)[1] <- "submitter_id"
  
  current <- left_join(full_map, patient, by ="submitter_id")
  
  current <- current[complete.cases(current),]
  
  current$project <- stringr::str_extract(current$assay, str_glue("{proj}"))
  
  
  # for PRAD ONLY
  #patient <- clin %>%
  #  dplyr::select(c("patientID","pathology_T_stage","pathology_N_stage","admin.file_uuid","patient.bcr_patient_uuid"))
  #current$pathology_M_stage <- NA
  #current <- current[c(1,2,3,4,5,6,10,7,8,9)]
  
  write.csv(current, file =str_glue("~/storage/MAE_analysis/clinical/full_map/{proj}_MAE_clinical_data"))
  
  rm(clin)
  rm(current)
  rm(patient)
  
}

# merge it all together now
setwd("~/storage/MAE_analysis/clinical/full_map")

BLCA <- as_tibble(data.table::fread("BLCA_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
BRCA <- as_tibble(data.table::fread("BRCA_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
COAD <- as_tibble(data.table::fread("COAD_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
ESCA <- as_tibble(data.table::fread("ESCA_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
HNSC <- as_tibble(data.table::fread("HNSC_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
KIRC <- as_tibble(data.table::fread("KIRC_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
KIRP <- as_tibble(data.table::fread("KIRP_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
LIHC <- as_tibble(data.table::fread("LIHC_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
LUAD <- as_tibble(data.table::fread("LUAD_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
LUSC <- as_tibble(data.table::fread("LUSC_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
PRAD <- as_tibble(data.table::fread("PRAD_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
STAD <- as_tibble(data.table::fread("STAD_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
THCA <- as_tibble(data.table::fread("THCA_MAE_clinical_data", header = TRUE)) %>%
  column_to_rownames("V1")
dat_clin <- as.data.frame(rbind(BLCA,BRCA,COAD,ESCA,HNSC,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,STAD,THCA))

rm(BLCA)
rm(BRCA)
rm(COAD)
rm(ESCA)
rm(HNSC)
rm(KIRC)
rm(KIRP)
rm(LIHC)
rm(LUAD)
rm(LUSC)
rm(PRAD)
rm(STAD)
rm(THCA)

write.csv(dat_clin, "Complete_sample_map.csv")
