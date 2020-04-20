# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)
library(TCGAbiolinks)
library(clusterProfiler)

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
stages <- c("NX", "N1", "N2", "N3", "M")

stage <- stages[5]
proj <- projects[9]
for(proj in projects){
  for(stage in stages){
  
  if(file.exists(str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/res/{stage}/{proj}_{stage}_DE.csv"))){
    
    dat <- data.table::fread(str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/res/{stage}/{proj}_{stage}_DE.csv"))
    dat <- as.data.frame(dat)
    dat <- column_to_rownames(dat, "V1")
    dat$ENSEMBL <- substr(rownames(dat), 1, 15)
    
    #dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
    
    dat <- dat %>%
      dplyr::filter(padj < 0.05) 
  
    #library(TCGAbiolinks)
    # Enrichment Analysis EA
    # Gene Ontology (GO) and Pathway enrichment by DEGs list
    #bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
    Genelist <- clusterProfiler::bitr(dat$ENSEMBL, OrgDb = org.Hs.eg.db ,fromType = "ENSEMBL", toType = "SYMBOL", drop = TRUE)
    
    
    if(length(Genelist$ENSEMBL) >0){
    system.time(ansEA <- TCGAanalyze_EAcomplete(TFname=str_glue("DEA genes {proj} {stage} Vs Normal"),Genelist$SYMBOL))
    
    # Enrichment Analysis EA (TCGAVisualize)
    # Gene Ontology (GO) and Pathway enrichment barPlot
    savR<- str_glue("{proj}_{stage}_PE_")
    TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                            GOBPTab = ansEA$ResBP,
                            GOCCTab = ansEA$ResCC,
                            GOMFTab = ansEA$ResMF,
                            PathTab = ansEA$ResPat,
                            nRGTab = Genelist, 
                            nBar = 10,
                            filename=paste0("~/storage/Machine_Learning/DE_TNM_by_Class/res/Viz/",
                                            savR,"TCGA_viz",".pdf"))
    }
    else
      next
    
    }
    else
      next
  
  }
 
}

# filtered DGE

stage <- stages[1]
proj <- projects[1]

for(proj in projects){
  foo <- dat[FALSE,]
  print(proj)
  for(stage in stages){
    print(stage)
    if(proj =="TCGA-LUSC" || stage == "N3"){
      next
    } else
    if(file.exists(str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/res/{stage}/{proj}_{stage}_DE.csv"))){
      
      dat <- data.table::fread(str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/res/{stage}/{proj}_{stage}_DE.csv"))
      dat <- as.data.frame(dat)
      dat <- column_to_rownames(dat, "V1")
      dat$ENSEMBL <- substr(rownames(dat), 1, 15)
      
      #dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
      
      dat <- dat %>%
        dplyr::filter(padj < 0.05)

      dat$stage <- stage
      
      foo <- rbind(dat,foo)
      dat <- foo
      
    }
    else 
      next
    
  }
  write.csv(dat, file = str_glue("~/storage/Machine_Learning/DE_TNM_by_Class/res/Filtered_DGE/Combined_Filtered_{proj}_DE.csv")) 
  
}

# use PCA to see separation between clusters


  

# conduct DE analysis for each one of the classes against Normal for each cancer type
# save dds objects, use old code to make Pathway enrichment

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] candisc_0.8-0               heplots_1.3-5               car_3.0-2                   carData_3.0-2               corrplot_0.84              
# [6] wesanderson_0.3.6           limma_3.38.3                psych_1.8.12                scales_1.0.0                GEOquery_2.50.5            
# [11] hexbin_1.27.2               vsn_3.50.0                  pheatmap_1.0.12             EnhancedVolcano_1.0.1       ggrepel_0.8.0              
# [16] org.Hs.eg.db_3.7.0          AnnotationDbi_1.44.0        DESeq2_1.22.2               SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
# [21] BiocParallel_1.16.6         matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
# [26] IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0         clusterProfiler_3.10.1      dendsort_0.3.3             
# [31] ggpubr_0.2                  magrittr_1.5                WGCNA_1.66                  fastcluster_1.1.25          dynamicTreeCut_1.63-1      
# [36] forcats_0.4.0               stringr_1.4.0               dplyr_0.8.0.1               purrr_0.3.2                 readr_1.3.1                
# [41] tidyr_0.8.3                 tibble_2.1.1                ggplot2_3.1.0               tidyverse_1.2.1            
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_0.2.5       robust_0.4-18          RSQLite_2.1.1          htmlwidgets_1.3        grid_3.5.1             munsell_0.5.0         
# [7] codetools_0.2-15       preprocessCore_1.44.0  withr_2.1.2            colorspace_1.4-1       GOSemSim_2.8.0         knitr_1.22            
# [13] rstudioapi_0.10        robustbase_0.93-4      DOSE_3.8.2             urltools_1.7.2         GenomeInfoDbData_1.2.0 mnormt_1.5-5          
# [19] polyclip_1.10-0        bit64_0.9-7            farver_1.1.0           generics_0.0.2         xfun_0.5               R6_2.4.0              
# [25] doParallel_1.0.14      locfit_1.5-9.1         bitops_1.0-6           fgsea_1.8.0            gridGraphics_0.3-0     assertthat_0.2.1      
# [31] ggraph_1.0.2           nnet_7.3-12            enrichplot_1.2.0       gtable_0.3.0           affy_1.60.0            rlang_0.3.3           
# [37] genefilter_1.64.0      splines_3.5.1          lazyeval_0.2.2         acepack_1.4.1          impute_1.56.0          broom_0.5.1           
# [43] europepmc_0.3          checkmate_1.9.1        yaml_2.2.0             BiocManager_1.30.4     reshape2_1.4.3         abind_1.4-5           
# [49] modelr_0.1.4           backports_1.1.3        qvalue_2.14.1          Hmisc_4.2-0            tools_3.5.1            ggplotify_0.0.3       
# [55] affyio_1.52.0          RColorBrewer_1.1-2     ggridges_0.5.1         Rcpp_1.0.1             plyr_1.8.4             base64enc_0.1-3       
# [61] progress_1.2.0         zlibbioc_1.28.0        RCurl_1.95-4.12        prettyunits_1.0.2      rpart_4.1-13           viridis_0.5.1         
# [67] cowplot_0.9.4          haven_2.1.0            cluster_2.0.7-1        data.table_1.12.0      DO.db_2.9              openxlsx_4.1.0        
# [73] triebeard_0.3.0        mvtnorm_1.0-10         hms_0.4.2              xtable_1.8-3           XML_3.98-1.19          rio_0.5.16            
# [79] readxl_1.3.1           gridExtra_2.3          compiler_3.5.1         crayon_1.3.4           htmltools_0.3.6        pcaPP_1.9-73          
# [85] Formula_1.2-3          geneplotter_1.60.0     rrcov_1.4-7            lubridate_1.7.4        DBI_1.0.0              tweenr_1.0.1          
# [91] MASS_7.3-51.1          Matrix_1.2-15          cli_1.1.0              igraph_1.2.4           pkgconfig_2.0.2        fit.models_0.5-14     
# [97] rvcheck_0.1.3          foreign_0.8-71         xml2_1.2.0             foreach_1.4.4          annotate_1.60.1        XVector_0.22.0        
# [103] rvest_0.3.2            digest_0.6.18          cellranger_1.1.0       fastmatch_1.1-0        htmlTable_1.13.1       curl_3.3              
# [109] nlme_3.1-137           jsonlite_1.6           viridisLite_0.3.0      pillar_1.3.1           lattice_0.20-38        httr_1.4.0            
# [115] DEoptimR_1.0-8         survival_2.43-3        GO.db_3.7.0            glue_1.3.1             zip_2.0.1              UpSetR_1.3.3          
# [121] iterators_1.0.10       bit_1.1-14             ggforce_0.2.1          stringi_1.4.3          blob_1.1.1             latticeExtra_0.6-28   
# [127] memoise_1.1.0
