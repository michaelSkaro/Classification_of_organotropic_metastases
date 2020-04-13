# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)


# Get all of the TNM staging with the barcodes for each of the cancer types:
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-HNSC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")
projects <- c("BLCA","BRCA","COAD","HNSC","LUAD","LUSC","PRAD","STAD","THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

proj <- projects[1]
setwd("~/CSBL_shared/clinical/TCGAbiolinks")
for(proj in projects){
  
  # load the R data file
  
  load(str_glue("~/CSBL_shared/clinical/TCGAbiolinks/{proj}_clin.Rdata"))
  
  if(exists("BLCA_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    BLCA_clinical <- BLCA_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    BLCA_clinical$N_stage_MS <- str_extract(BLCA_clinical$stage_event_tnm_categories, "N\\w")
    BLCA_clinical$M_stage_MS <- str_extract(BLCA_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    BLCA_clinical <-  BLCA_clinical %>%
      mutate(BLCA_clinical, TNM_stage_MS = ifelse(BLCA_clinical$M_stage_MS == "M1", 
             "M", BLCA_clinical$N_stage_MS))
      
    BLCA_clinical <- BLCA_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
  # write out clincal file as a csv so we can change the project files  
   
  write.csv(BLCA_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
  
    
  #remove clincal file
  
  rm(BLCA_clinical)
  
  }
  
  if(exists("BRCA_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    BRCA_clinical <- BRCA_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    BRCA_clinical$N_stage_MS <- str_extract(BRCA_clinical$stage_event_tnm_categories, "N\\w")
    BRCA_clinical$M_stage_MS <- str_extract(BRCA_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    BRCA_clinical <-  BRCA_clinical %>%
      mutate(BRCA_clinical, TNM_stage_MS = ifelse(BRCA_clinical$M_stage_MS == "M1", 
                                                  "M", BRCA_clinical$N_stage_MS))
    
    BRCA_clinical <- BRCA_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(BRCA_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(BRCA_clinical)
    
  }
  
  if(exists("COAD_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    COAD_clinical <- COAD_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    COAD_clinical$N_stage_MS <- str_extract(COAD_clinical$stage_event_tnm_categories, "N\\w")
    COAD_clinical$M_stage_MS <- str_extract(COAD_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    COAD_clinical <-  COAD_clinical %>%
      mutate(COAD_clinical, TNM_stage_MS = ifelse(COAD_clinical$M_stage_MS == "M1", 
                                                  "M", COAD_clinical$N_stage_MS))
    
    COAD_clinical <- COAD_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files 
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(COAD_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(COAD_clinical)
    
    
  }
  
  if(exists("HNSC_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    HNSC_clinical <- HNSC_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    HNSC_clinical$N_stage_MS <- str_extract(HNSC_clinical$stage_event_tnm_categories, "N\\w")
    HNSC_clinical$M_stage_MS <- str_extract(HNSC_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    HNSC_clinical <-  HNSC_clinical %>%
      mutate(HNSC_clinical, TNM_stage_MS = ifelse(HNSC_clinical$M_stage_MS == "M1", 
                                                  "M", HNSC_clinical$N_stage_MS))
    
    HNSC_clinical <- HNSC_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(HNSC_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(HNSC_clinical)
    
    
  }
  
  if(exists("LUAD_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    LUAD_clinical <- LUAD_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    LUAD_clinical$N_stage_MS <- str_extract(LUAD_clinical$stage_event_tnm_categories, "N\\w")
    LUAD_clinical$M_stage_MS <- str_extract(LUAD_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    LUAD_clinical <-  LUAD_clinical %>%
      mutate(LUAD_clinical, TNM_stage_MS = ifelse(LUAD_clinical$M_stage_MS == "M1", 
                                                  "M", LUAD_clinical$N_stage_MS))
    
    LUAD_clinical <- LUAD_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(LUAD_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(LUAD_clinical)
  }
  
  if(exists("LUSC_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    LUSC_clinical <- LUSC_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    LUSC_clinical$N_stage_MS <- str_extract(LUSC_clinical$stage_event_tnm_categories, "N\\w")
    LUSC_clinical$M_stage_MS <- str_extract(LUSC_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    LUSC_clinical <-  LUSC_clinical %>%
      mutate(LUSC_clinical, TNM_stage_MS = ifelse(LUSC_clinical$M_stage_MS == "M1", 
                                                  "M", LUSC_clinical$N_stage_MS))
    
    LUSC_clinical <- LUSC_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(LUSC_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(LUSC_clinical)  
  }
  
  if(exists("PRAD_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    PRAD_clinical <- PRAD_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    PRAD_clinical$N_stage_MS <- str_extract(PRAD_clinical$stage_event_tnm_categories, "N\\w")
    PRAD_clinical$M_stage_MS <- str_extract(PRAD_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    PRAD_clinical <-  PRAD_clinical %>%
      mutate(PRAD_clinical, TNM_stage_MS = ifelse(PRAD_clinical$M_stage_MS == "M1", 
                                                  "M", PRAD_clinical$N_stage_MS))
    
    PRAD_clinical <- PRAD_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(PRAD_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(PRAD_clinical) 
  }
  
  if(exists("STAD_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    STAD_clinical <- STAD_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    STAD_clinical$N_stage_MS <- str_extract(STAD_clinical$stage_event_tnm_categories, "N\\w")
    STAD_clinical$M_stage_MS <- str_extract(STAD_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    STAD_clinical <-  STAD_clinical %>%
      mutate(STAD_clinical, TNM_stage_MS = ifelse(STAD_clinical$M_stage_MS == "M1", 
                                                  "M", STAD_clinical$N_stage_MS))
    
    STAD_clinical <- STAD_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(STAD_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(STAD_clinical) 
  }
  
  if(exists("THCA_clinical")){
    
    # extract barcode + event_stage_tnm_catagories using select command
    
    THCA_clinical <- THCA_clinical %>%
      dplyr::select(bcr_patient_barcode,stage_event_tnm_categories)
    
    
    # segregate the TNM staging into N0,N1,N2, N3 and M using mutate
    
    
    THCA_clinical$N_stage_MS <- str_extract(THCA_clinical$stage_event_tnm_categories, "N\\w")
    THCA_clinical$M_stage_MS <- str_extract(THCA_clinical$stage_event_tnm_categories, "M\\w")
    
    
    # add proper column TNM staging for each patient in N0, N1,N2,N3,M
    #NX means the nearby lymph nodes can't be assessed, for example, if they were previously removed.
    #N0 means nearby lymph nodes do not contain cancer.
    #N1, N2, N3: These numbers are based on the number of 
    #lymph nodes involved and how much cancer is found in them. 
    #The higher the N number, the greater the extent of the lymph node involvement.
    
    THCA_clinical <-  THCA_clinical %>%
      mutate(THCA_clinical, TNM_stage_MS = ifelse(THCA_clinical$M_stage_MS == "M1", 
                                                  "M", THCA_clinical$N_stage_MS))
    
    THCA_clinical <- THCA_clinical %>%
      dplyr::select(bcr_patient_barcode,TNM_stage_MS)
    
    
    # write out clincal file as a csv so we can change the project files  
    
    write.csv(THCA_clinical, file = str_glue("~/storage/Machine_Learning/TNM_stage/TCGA-{proj}_TNM.csv"))
    
    
    #remove clincal file
    
    rm(THCA_clinical) 
  }
  
}    
    

# Get all of the TNM staging with the barcodes for each of the cancer types:
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-HNSC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")


# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$sample_type == "Primary Tumor",]
tumor.samples <- tumor.samples[tumor.samples$project %in% projects,]

tumor.samples$bcr_patient_barcode <- substr(tumor.samples$barcode, 0, 12)
normal.samples$bcr_patient_barcode <- substr(normal.samples$barcode, 0, 12)
clinical$barcode_short <- substr(clinical$barcode, 0,16)




setwd("~/CSBL_shared/RNASeq/TCGA/counts")
#proj <- projects[4]
for(proj in projects){
  
    df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
      as_tibble() %>%
      tibble::column_to_rownames(var = "Ensembl")

    n<-dim(df.exp)[1]
    df.exp<-df.exp[1:(n-5),]

    df.exp.cancertype <- df.exp
    
    coldata.t<- tumor.samples[tumor.samples$project == proj,]
    coldata.n <- normal.samples[normal.samples$project == proj,]
  
    TNM <- data.table::fread(str_glue("~/storage/Machine_Learning/TNM_stage/{proj}_TNM.csv"), header = TRUE) %>%
      column_to_rownames("V1")
    
    coldata.t <- coldata.t %>%
      left_join(TNM, by= "bcr_patient_barcode")
    
    coldata.n$TNM_stage_MS <- "Normal"
  # read in TNM staging information and attach the TNM staging as a column
  
  coldata <- rbind(coldata.n, coldata.t)  
    
  coldata <- coldata[!is.na(coldata$TNM_stage_MS),]
  
  coldata.cancerType <- coldata
  
  All_stages <- unique(coldata$TNM_stage_MS)[2:7]
  
  
  
  for(TNM_stage in All_stages){
  
  coldata.t <- as.data.frame(coldata.cancerType[coldata.cancerType$TNM_stage_MS == TNM_stage,])
  
  # get rid of duplicates
  
  coldata.t <- coldata.t[!duplicated(coldata.t),]
  
  coldata <- as.data.frame(rbind(coldata.n, coldata.t))
  
  df.exp <- df.exp.cancertype
  
  df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
  
  
  rownames(coldata) <- coldata$barcode
  coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)
  
  rownames(coldata) <- sort(rownames(coldata))
  colnames(df.exp) <- sort(colnames(df.exp))
  
  dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ sample_type)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
  
  dds <- DESeq(dds)
  
  # conduct DE analysis for each one of the classes against Normal for each cancer type
  # save dds objects, use old code to make Pathway enrichment
  
  save(dds, file = str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/{proj}_{TNM_stage}_vs.Normal_DE_.RData"))
  
  res <- results(dds)
  res <- as.data.frame(res)
  resOrdered <- res[order(res$pvalue),]
  
  write.csv(resOrdered, file = str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/{proj}_{TNM_stage}_DE.csv"))
  
  rm(dds)
  rm(res)
  rm(resOrdered)
  
  }
  
}

# put together a list of the DGEs

# use PCA to see separation between clusters


# conduct DE analysis for each one of the classes against Normal for each cancer type
# save dds objects, use old code to make Pathway enrichment

