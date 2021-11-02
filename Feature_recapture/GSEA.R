###############################################
###############################################
###############################################
library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(cowplot)
# read in the res values 
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")


# gene set enrichment analysis for the DE genes
# hallmark human cancer gene set
all_gene_sets = msigdbr(species = "Homo sapiens")
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

for(proj in projects){
  print(str_glue("Processing: {proj}"))
  for(org in organs){
    print(str_glue("Analyzing: {org}"))
    if(file.exists(str_glue("{proj}_DE_met_{org}.csv"))){
      dat <- data.table::fread(str_glue("{proj}_DE_met_{org}.csv", header = TRUE)) %>%
        tibble::column_to_rownames("V1") %>%
        # dplyr::filter(padj <0.05) %>% This was done previous steps
        # dplyr::filter(abs(log2FoldChange) > 0.5) %>%
        tibble::rownames_to_column("ENSEMBL")
      # remove NAs if they exist
      dat <- dat[complete.cases(dat),]
      dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
      
      geneList = dat[,1]
      
      ## feature 2: named vector
      names(geneList) = dat[,3]
      
      ## feature 3: decreasing orde
      geneList = sort(geneList, decreasing = TRUE)
      
      
      ego <-clusterProfiler::enrichGO(gene  = geneList,
                                      OrgDb         = "org.Hs.eg.db",
                                      keyType       = 'ENSEMBL',
                                      ont           = "MF",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
      
      #Get only  the significant results from enrichment
      res <- ego@result
      res <- res[res$pvalue<0.05,]
      write.csv(res, file = str_glue("MF_ego_{proj}_{org}.csv"))
      
      
    }
  }
}



# Modified pathways and the overlaps
library(tibble)
library(dplyr)
library(plyr)
compare_list <- as.data.frame(data.table::fread("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/compare_list.csv", header = TRUE))

i <- 1
out <- as.data.frame(t(c("CancerType1"= NA,
                         "CancerType2"= NA,
                         "Seeding_Location" = NA,
                         "odds.ratio" = NA,
                         "Enriched Processes CT1" = NA,
                         "Enriched Processes CT2" = NA,
                         "intersection"= NA,
                         "p.value" = NA)))

for(i in 1:length(compare_list$CancerType1)){
  print(i)
  print(paste0("BP_ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  print(paste0("BP_ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  c1 <- data.table::fread(paste0("BP_ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  c2 <- data.table::fread(paste0("BP_ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  
  # add all of the metric to a dataframe. 
  
  
  go.obj <- newGeneOverlap(listA = c1$ID, listB = c2$ID, genome.size = 23393)
  go.obj <- testGeneOverlap(go.obj)
  df <- as.data.frame(t(c("CancerType1"= compare_list[i,1],
                          "CancerType2"= compare_list[i,2],
                          "Seeding_Location" = compare_list[i,3],
                          "Enriched Processes CT1" = length(c1$ID),
                          "Enriched Processes CT2" = length(c2$ID),
                          "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                          "intersection"= length(go.obj@intersection),
                          "p.value" = go.obj@pval)))
  out <- rbind(out,df)
  
  
}
out <- out[2:length(out$CancerType1),]
out <- out[order(out$p.value),]  
write.csv(out, file = "/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/res/Fisher_exact_test_MF_semantic.csv")
