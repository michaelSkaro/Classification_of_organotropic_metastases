# Remember to use the draw function.
# The editor note:Complexheatmap::draw()

#If you put Heatmap() inside a function or a for/if/while chunk, 
#you won’t see the heatmap after executing Heatmap(). In this case, 
#you need to use draw() function explicitly as follows. 


library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
setwd("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets")
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
files_in_dir <- list.files()
org <- organs[1]
library(ComplexHeatmap)
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")

proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]
for(proj in projects){
  
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
      # complete data set of IG selected features. The train test splits will be merged 
      val <- str_sub(proj, 6,9)
      dat <- rbind(train, test)
      organ_labels <- dat %>%dplyr::select(org)
      
      dat <- dat %>%
        dplyr::select(-org)
      dat <- as.data.frame(t(log2(dat)))
      names(dat) <- organ_labels[,1]
      
      dat[mapply(is.infinite, dat)] <- 0
      basemean = rowMeans(dat)
      type = gsub("s\\d+_", "", colnames(dat))
      dat <- as.matrix(dat)
      ha = HeatmapAnnotation(type = type, annotation_name_side ="right")
      pdf(str_glue("heatmap_{proj}_{org}_log2_rasteration.pdf"), width = 5, height = 5)
      ht_list = Heatmap(dat, name = " Log2 counts",
                        border = TRUE,
                        cluster_columns =TRUE,
                        cluster_rows = TRUE,
                        row_dend_reorder = TRUE,
                        column_dend_reorder = TRUE,
                        show_heatmap_legend = TRUE,
                        show_column_names = FALSE, show_row_names = FALSE,
                        show_column_dend = TRUE, show_row_dend = FALSE, 
                        use_raster = TRUE, raster_resize = TRUE,
                        top_annotation = ha,
                        column_title = paste0(val," ", "Metastasize to ",str_glue("{org}")),
                        row_title = "Selected features")
      #top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), height = unit(2, "mm"))),
      #width = unit(15, "mm"))
      
      draw(ht_list)
      
      
      dev.off()
      
    }
  }
}

library(UpSetR)
library(stringr)
library(GeneOverlap)


for(proj in projects){
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
      # complete data set of IG selected features. The train test splits will be merged 
      val <- str_sub(proj, 6,9)
      dat <- rbind(train, test)
      organ_labels <- dat %>%dplyr::select(org)
      dat <- dat %>%
        dplyr::select(-org)
      
      # cancer type, metastatic location organ specific gene set enrichment 
      
      ego <-clusterProfiler::enrichGO(gene  = substring(colnames(dat),1,15),
                                      OrgDb         = "org.Hs.eg.db",
                                      keyType       = 'ENSEMBL',
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
      
      #Get only  the significant results from enrichment
      res <- ego@result
      #res <- res[res$p.adjust<0.05,]
      res <- res[res$pvalue<0.05,]
      write.csv2(res, file = str_glue("ego_{proj}_{org}.csv"))
      
    }
  }
}

# semantic analysis of pathways and the overlaps
library(tibble)
library(dplyr)
library(plyr)
compare_list <- as.data.frame(data.table::fread("compare_list.csv", header = TRUE))

i <- 1
out <- as.data.frame(t(c("CancerType1"= NA,
                         "CancerType2"= NA,
                         "Seeding_Location" = NA,
                         "odds.ratio" = NA,
                         "Enriched Processes CT1" = NA,
                         "Enriched Processes CT2" = NA,
                         "intersection"= NA,
                         "p.value" = NA)))
for(i in 1:length(CT1)){
  print(i)
  print(paste0("ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  print(paste0("ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  c1 <- data.table::fread(paste0("ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  c2 <- data.table::fread(paste0("ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  
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
out <- out[sort(out$p.value),]  
write.csv2(out, file = "/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/Fisher_exact_test_semantic.csv")
