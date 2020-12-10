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
  
sessionInfo()

library(simplifyEnrichment)
library(stringr)
library(tibble)
dat <- data.table::fread("ego_TCGA-HNSC_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")

mat = GO_similarity(dat$ID, ont = "BP")
simplifyGO(mat, method = "binary_cut",column_title = "Biological Processes clustered by Sematic Similarity", plot = TRUE)

dat <- data.table::fread("ego_TCGA-BLCA_Lymph_Node.csv",header = TRUE) %>%
  column_to_rownames("V1")
mat = GO_similarity(dat$ID, ont = "BP")
simplifyGO(mat, method = "binary_cut",column_title = "Biological Processes clustered by Sematic Similarity", plot = TRUE)

# make upsetR plot to finish it off

library(UpSetR)
seed_loc <- c("Lymph_Node","Lung","Bone","Lung")

#Bone
BLCA_Bone <- data.table::fread("ego_TCGA-BLCA_Bone.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA_Bone <- data.table::fread("ego_TCGA-BRCA_Bone.csv", header = TRUE) %>%
  column_to_rownames("V1")
listInput <- list('BLCA Bone' = BLCA_Bone$ID ,'BRCA Bone' = BRCA_Bone$ID)




# Liver 
BLCA_Liver <- data.table::fread("ego_TCGA-BLCA_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA_Liver <- data.table::fread("ego_TCGA-BRCA_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
LIHC_Liver <- data.table::fread("ego_TCGA-LIHC_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
COAD_Liver <- data.table::fread("ego_TCGA-COAD_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")


library(UpSetR)
listInput <- list('BLCA Liver' = BLCA_Liver$ID ,'BRCA Liver' = BRCA_Liver$ID ,
                  'LIHC Liver' = LIHC_Liver$ID,'COAD Liver' = COAD_Liver$ID)

UpSetR::upset(fromList(listInput), order.by ="freq")


# Lung
BLCA<- data.table::fread("ego_TCGA-BLCA_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA<- data.table::fread("ego_TCGA-BRCA_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
LIHC<- data.table::fread("ego_TCGA-LIHC_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
HNSC<- data.table::fread("ego_TCGA-HNSC_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
LUAD<- data.table::fread("ego_TCGA-LUAD_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")

library(UpSetR)
listInput <- list('BLCA Lung' = BLCA$ID ,'BRCA Lung' = BRCA$ID , "HNSC Lung"=HNSC$ID,
                  'LIHC Lung' = LIHC$ID,'LUAD Liver' = LUAD$ID)


#Lymph_Node
BLCA_Lymph_Node <- data.table::fread("ego_TCGA-BLCA_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")
HNSC_Lymph_Node <- data.table::fread("ego_TCGA-HNSC_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")
listInput <- list('BLCA Lymph_Node' = BLCA_Lymph_Node$ID ,'HNSC Lymph_Node' = HNSC_Lymph_Node$ID)

UpSetR::upset(fromList(listInput), order.by ="freq")





# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] UpSetR_1.4.0              ComplexHeatmap_2.7.1      plyr_1.8.6                stringr_1.4.0            
# [5] GeneOverlap_1.24.0        org.Hs.eg.db_3.11.4       AnnotationDbi_1.50.3      IRanges_2.22.2           
# [9] S4Vectors_0.26.1          Biobase_2.48.0            simplifyEnrichment_0.99.5 BiocGenerics_0.34.0      
# [13] tibble_3.0.4              dplyr_1.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.14.0           colorspace_2.0-0       rjson_0.2.20           ellipsis_0.3.1         ggridges_0.5.2        
# [6] circlize_0.4.11        qvalue_2.20.0          GlobalOptions_0.1.2    clue_0.3-57            rstudioapi_0.13       
# [11] farver_2.0.3           urltools_1.7.3         graphlayouts_0.7.1     ggrepel_0.8.2          bit64_4.0.5           
# [16] fansi_0.4.1            scatterpie_0.1.5       xml2_1.3.2             splines_4.0.2          GOSemSim_2.14.2       
# [21] knitr_1.30             polyclip_1.10-0        jsonlite_1.7.1         Cairo_1.5-12.2         cluster_2.1.0         
# [26] GO.db_3.11.4           png_0.1-7              ggforce_0.3.2          BiocManager_1.30.10    compiler_4.0.2        
# [31] httr_1.4.2             rvcheck_0.1.8          assertthat_0.2.1       Matrix_1.2-18          cli_2.1.0             
# [36] tweenr_1.0.1           htmltools_0.5.0        prettyunits_1.1.1      tools_4.0.2            igraph_1.2.6          
# [41] NLP_0.2-1              gtable_0.3.0           glue_1.4.2             reshape2_1.4.4         DO.db_2.9             
# [46] tinytex_0.27           fastmatch_1.1-0        Rcpp_1.0.5             enrichplot_1.8.1       slam_0.1-47           
# [51] vctrs_0.3.4            ggraph_2.0.3           xfun_0.19              lifecycle_0.2.0        clusterProfiler_3.16.1
# [56] gtools_3.8.2           DOSE_3.14.0            europepmc_0.4          MASS_7.3-53            scales_1.1.1          
# [61] tidygraph_1.2.0        hms_0.5.3              RColorBrewer_1.1-2     yaml_2.2.1             memoise_1.1.0         
# [66] gridExtra_2.3          ggplot2_3.3.2          downloader_0.4         triebeard_0.3.0        stringi_1.5.3         
# [71] RSQLite_2.2.1          caTools_1.18.0         BiocParallel_1.22.0    shape_1.4.5            bitops_1.0-6          
# [76] rlang_0.4.8            pkgconfig_2.0.3        matrixStats_0.57.0     evaluate_0.14          lattice_0.20-41       
# [81] purrr_0.3.4            cowplot_1.1.0          bit_4.0.4              tidyselect_1.1.0       magrittr_1.5          
# [86] R6_2.5.0               gplots_3.1.0           generics_0.1.0         DBI_1.1.0              pillar_1.4.6          
# [91] proxyC_0.1.5           crayon_1.3.4           KernSmooth_2.23-17     rmarkdown_2.5          viridis_0.5.1         
# [96] GetoptLong_1.0.4       progress_1.2.2         data.table_1.13.2      blob_1.2.1             digest_0.6.27         
# [101] tm_0.7-7               tidyr_1.1.2            gridGraphics_0.5-0     RcppParallel_5.0.2     munsell_0.5.0         
# [106] viridisLite_0.3.0      ggplotify_0.0.5 
