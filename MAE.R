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
  
  
  
  
  
    
}


# methylation


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








