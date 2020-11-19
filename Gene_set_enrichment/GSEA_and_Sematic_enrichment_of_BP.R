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
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph Node", "Pelvis", "Prostate")
files_in_dir <- list.files()
org <- organs[1]
# show expression Profiler differences
for(proj in projects){
  
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
    train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
    test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
    # complete data set of IG selected features. The train test splits will be merged 
    
    dat <- rbind(train, test)
    organ_labels <- dat %>%dplyr::select(org)
    dat <- dat %>%
      dplyr::select(-org)
    dat <- as.data.frame(t(log2(dat)))
    
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
                      column_title = str_glue("{proj} Metastasize to {org}"),
                      row_title = "Selected features")
    #top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), height = unit(2, "mm"))),
    #width = unit(15, "mm"))
    
    draw(ht_list)
    
    
    
    }
  }
}

# Gene set enrichment analysis

for(proj in projects){
  
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
    train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
    test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
    # complete data set of IG selected features. The train test splits will be merged 
    
    dat <- rbind(train, test)
    organ_labels <- dat %>%dplyr::select(org)
    dat <- dat %>%
      dplyr::select(-org)

    # cancer type, metastatic location organ specific gene set enrichment 
    
    ego <-clusterProfiler::enrichGO(gene  = substring(rownames(dat),1,15),
                                OrgDb         = "org.Hs.eg.db",
                                keyType       = 'ENSEMBL',
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
    
    #Get only  the significant results from enrichment
    res <- ego@result
    res <- res[res$p.adjust<0.05,]
    write.csv2(result, file = str_glue("ego_{proj}_{org}.csv")
    }
  }
}






