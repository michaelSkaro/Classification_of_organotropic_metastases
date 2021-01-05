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

# Analyzing the selected features for the MOT project in context to
 # the selected features

library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tibble)
library(DESeq2)
library(fgsea)
setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels")
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[4]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph Node", "Pelvis", "Prostate")
files_in_dir <- list.files()
org <- organs[5]
library(ComplexHeatmap)


annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")


normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv", stringsAsFactors = TRUE) %>%
  tibble::column_to_rownames("V1")

# Oncogenic pathways for enrichment analysis
m_t2g <- msigdbr(species = "Homo sapiens", category = "C4") %>% 
  dplyr::select(gs_name, entrez_gene)



for(proj in projects){
    if(file.exists(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_metastatic_data_RNAseq.csv"))){
      
      met_annot <- data.table::fread(str_glue("{proj}_metastatic_data_RNAseq.csv"))
      
      coldata.m <- met_annot %>%
        dplyr::select(c(barcode,colnames(met_annot)[colnames(met_annot) %in% organs]))
      
      # read in counts data, subset counts data to the met tumors for each comparison
      df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
        as_tibble() %>%
        tibble::column_to_rownames(var = "Ensembl")
      
      n<-dim(df.exp)[1]
      df.exp<-df.exp[1:(n-5),]
      
      coldata.t<- tumor.samples[tumor.samples$project == proj,]
      coldata.n <- normal.samples[normal.samples$project == proj,]
      coldata <- rbind(coldata.n, coldata.t, fill = TRUE)
      coldata <- left_join(coldata, coldata.m, by = "barcode")
      coldata[is.na(coldata)] <- 0
      
      
      df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
      rownames(coldata) <- coldata$barcode
      coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)
      
      rownames(coldata) <- sort(rownames(coldata))
      colnames(df.exp) <- sort(colnames(df.exp))
      
      
      met_locs <- colnames(coldata.m)[-1]
      
      
      for(i in 1:length(met_locs)){
        
        # Customize df.exp and cioldata
        print(proj)
        print(met_locs[i])
        org_til<- met_locs[2]
        
        coldata.m.l <- coldata %>%
          filter(
            .data[[met_locs[[1]]]] == 1)
        coldata.n <- coldata %>%
         filter(sample_type == "Solid_Tissue_Normal")
        
        coldata.n.m. <- rbind(coldata.n, coldata.m.l)
        
        df.exp.t <- df.exp[ ,colnames(df.exp) %in% coldata.n.m.$barcode]
        rownames(coldata.n.m.) <- sort(rownames(coldata.n.m.))
        colnames(df.exp.t) <- sort(colnames(df.exp.t))
        
        dds <- DESeqDataSetFromMatrix(countData = df.exp.t, colData = coldata.n.m., design = ~ sample_type)
        
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]
        dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
        
        dds <- DESeq(dds)
        
        res <- results(dds)
        
        resOrdered <- res[order(res$pvalue),]
        resOrdered <- as.data.frame(resOrdered)
        
        
        
        res <- as.data.frame(res)
        
        write.csv(resOrdered, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_met_{org_til}_RNAseq.csv"))
        save(dds, file = str_glue("{proj}_{org_til}_DE_met.RData"))
        
      }
    }
}

# Build a superRes object that add the Project and Org columns as annotation. 
# This will allow us to iteratively complete all of the intersections, both by 
# project and by seeding locaiton 

setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels")
proj <- projects[1]
org <- organs[1]
organs[5] <- "Lymph_Node"
# make a dataframe with the same columns as the res object + 2 columns will have
# add an empty line that will be removed later on
#read in a res then just fall into the for loop
superRes <- res[1,]
superRes[1,] <- NA

for(proj in projects){
  for(org in organs){
  if(file.exists(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_met_{org}_RNAseq.csv"))){
    res <- data.table::fread(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_met_{org}_RNAseq.csv")) %>%
      tibble::column_to_rownames("V1")
    
    # filter for FC ≥ |1.0|
  
    res <- res %>%
      filter(abs(log2FoldChange) >.5) %>%
      filter(padj <= .05) %>%
      rownames_to_column("ENSEMBL")
    res$ENSEMBL <- substr(res$ENSEMBL,1,15)
    
    res$proj <- rep(x = proj, times = length(res$ENSEMBL))
    res$org <- rep(x = org, times = length(res$ENSEMBL))
    # add the columns to the res object as annotation
    
    # rbind the res object to the SuperRes.
    
    superRes <- as.data.frame(rbind(superRes,res))
    superRes <- superRes[complete.cases(superRes),]
    
    }
  }  
}




# complete the intersections.
# iterate over the projects
# make each of the of the comparisons
library(tibble)
library(dplyr)
library(plyr)
compare_list <- as.data.frame(data.table::fread("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/compare_list.csv", header = TRUE))

out <- as.data.frame(t(c("CancerType1"= NA,
                         "CancerType2"= NA,
                         "Seeding_Location" = NA,
                         "odds.ratio" = NA,
                         "DE genes CT1" = NA,
                         "DE genes CT2" = NA,
                         "intersection"= NA,
                         "p.value" = NA)))


library(UpSetR)
library(stringr)
library(GeneOverlap)
for(i in 1:length(compare_list$CancerType1)){
  
  comp <- superRes %>%
    filter(proj %in% compare_list[i,1:2]) %>%
    filter(org %in% compare_list[i,3])
  c1 <- comp[comp$proj==compare_list[i,1],]
  c2 <- comp[comp$proj==compare_list[i,2],]
  go.obj <- newGeneOverlap(listA = c1$ENSEMBL, 
                           listB = c2$ENSEMBL, 
                           genome.size = 60483)
  go.obj <- testGeneOverlap(go.obj)
  df <- as.data.frame(t(c("CancerType1"= compare_list[i,1],
                          "CancerType2"= compare_list[i,2],
                          "Seeding_Location" = compare_list[i,3],
                          "DE genes CT1" = length(c1$ENSEMBL),
                          "DE genes CT2" = length(c2$ENSEMBL),
                          "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                          "intersection"= length(go.obj@intersection),
                          "p.value" = go.obj@pval)))
  out <- rbind(out,df)
  out <- out[complete.cases(out),]
  
}

out$p.value <- as.numeric(out$p.value)
out <- out[order(out$p.value),] 
write.csv2(out, file = "/mnt/storage/mskaro1/Machine_Learning/Fisher_Exact_testDE_genes_tumors_met_same_loc.csv")


# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] circlize_0.4.11           ComplexHeatmap_2.7.1      forcats_0.5.0            
# [4] dplyr_1.0.2               purrr_0.3.4               readr_1.4.0              
# [7] tidyr_1.1.2               ggplot2_3.3.2             tidyverse_1.3.0          
# [10] UpSetR_1.4.0              tibble_3.0.4              stringr_1.4.0            
# [13] org.Hs.eg.db_3.11.4       AnnotationDbi_1.50.3      IRanges_2.22.2           
# [16] S4Vectors_0.26.1          Biobase_2.48.0            simplifyEnrichment_0.99.5
# [19] BiocGenerics_0.34.0      
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1           backports_1.2.0        fastmatch_1.1-0        plyr_1.8.6            
# [5] igraph_1.2.6           proxyC_0.1.5           splines_4.0.2          BiocParallel_1.22.0   
# [9] urltools_1.7.3         digest_0.6.27          htmltools_0.5.0        GOSemSim_2.14.2       
# [13] viridis_0.5.1          GO.db_3.11.4           fansi_0.4.1            magrittr_1.5          
# [17] memoise_1.1.0          tm_0.7-7               cluster_2.1.0          graphlayouts_0.7.1    
# [21] modelr_0.1.8           RcppParallel_5.0.2     matrixStats_0.57.0     enrichplot_1.8.1      
# [25] prettyunits_1.1.1      colorspace_2.0-0       blob_1.2.1             rvest_0.3.6           
# [29] ggrepel_0.8.2          haven_2.3.1            xfun_0.19              crayon_1.3.4          
# [33] jsonlite_1.7.1         scatterpie_0.1.5       glue_1.4.2             polyclip_1.10-0       
# [37] gtable_0.3.0           GetoptLong_1.0.4       shape_1.4.5            scales_1.1.1          
# [41] DOSE_3.14.0            DBI_1.1.0              Rcpp_1.0.5             viridisLite_0.3.0     
# [45] progress_1.2.2         clue_0.3-57            gridGraphics_0.5-0     bit_4.0.4             
# [49] europepmc_0.4          httr_1.4.2             fgsea_1.14.0           RColorBrewer_1.1-2    
# [53] ellipsis_0.3.1         pkgconfig_2.0.3        farver_2.0.3           dbplyr_2.0.0          
# [57] ggplotify_0.0.5        tidyselect_1.1.0       rlang_0.4.8            reshape2_1.4.4        
# [61] munsell_0.5.0          cellranger_1.1.0       tools_4.0.2            downloader_0.4        
# [65] cli_2.1.0              generics_0.1.0         RSQLite_2.2.1          broom_0.7.2           
# [69] ggridges_0.5.2         evaluate_0.14          yaml_2.2.1             knitr_1.30            
# [73] bit64_4.0.5            fs_1.5.0               tidygraph_1.2.0        ggraph_2.0.3          
# [77] slam_0.1-47            DO.db_2.9              xml2_1.3.2             compiler_4.0.2        
# [81] rstudioapi_0.13        png_0.1-7              reprex_0.3.0           tweenr_1.0.1          
# [85] stringi_1.5.3          lattice_0.20-41        Matrix_1.2-18          vctrs_0.3.4           
# [89] pillar_1.4.6           lifecycle_0.2.0        BiocManager_1.30.10    triebeard_0.3.0       
# [93] GlobalOptions_0.1.2    data.table_1.13.2      cowplot_1.1.0          qvalue_2.20.0         
# [97] R6_2.5.0               gridExtra_2.3          MASS_7.3-53            assertthat_0.2.1      
# [101] rjson_0.2.20           withr_2.3.0            hms_0.5.3              clusterProfiler_3.16.1
# [105] rmarkdown_2.5          rvcheck_0.1.8          Cairo_1.5-12.2         ggforce_0.3.2         
# [109] NLP_0.2-1   
  
