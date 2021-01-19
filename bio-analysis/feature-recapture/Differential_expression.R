# Analyzing the selected features for the MOT project in context to
 # the selected features


library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tibble)
library(DESeq2)
setwd("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets")
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
files_in_dir <- list.files()
org <- organs[1]
library(ComplexHeatmap)
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")

annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")


normal.samples <- clinical[sample_type == "Solid Tissue Normal"]

tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv", stringsAsFactors = TRUE) %>%
  tibble::column_to_rownames("V1")


proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]
for(proj in projects){
  
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      
      
      
      # just need the features
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      
      # read in the values
      dat <- data.table::fread((str_glue("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_metastatic_data_RNAseq.csv")))
      
     
      
      # select barcode and the organ columns to make the tumor design matrix
      
     
      
      # load the normal data for the project 
      
      normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
      normal.samples <- normal.samples[project == proj]
      # read in the counts data, filter for just normal data and the selected genes
      df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
        as_tibble() %>%
        tibble::column_to_rownames(var = "Ensembl")
      coldata.n <- normal.samples[normal.samples$project == proj,]
      df.exp <- df.exp[ ,colnames(df.exp) %in% coldata.n$barcode]
      df.exp <- df.exp[rownames(df.exp) %in% rownames(dat),]
      
      # make the design matrix
      coldata.n <- as.data.frame(t(coldata.n[,1:2]))
      header.true <- function(df) {
        names(df) <- as.character(unlist(df[2,]))
        df[-2,]
      }
      
      coldata.n <- header.true(coldata.n)
      metastatic_status <- colnames(coldata.n)
      coldata.n <- as.data.frame(sort(coldata.n))
      coldata.n <- as.data.frame(t(coldata.n))
      coldata.n$Metastatic_status <- metastatic_status
      coldata.n <- coldata.n %>%
        select(-project)
      design.mat <- rbind(design.mat,coldata.n)
      
      design.mat$disease_status <- ifelse(design.mat$Metastatic_status == "Solid Tissue Normal", "Normal", 
                                          "Tumor") 
      
      # make the df.exp r bind the normal samples to the tumors and make sure the columns and rownames match
      
      colnames(df.exp) <- rownames(coldata.n)
      
      # rownames to column called ENSEMBL to join on
      
      df.exp <- df.exp %>%
        tibble::rownames_to_column("ENSEMBL")
      dat <- dat %>%
        tibble::rownames_to_column("ENSEMBL")
      
      df.exp <- dplyr::left_join(df.exp, dat, by = "ENSEMBL") %>%
        column_to_rownames("ENSEMBL")
      
      # DEanalysis with DEseq2, use contrasts to find up and down regulated transcripts
      rownames(design.mat) <- str_replace_all(rownames(design.mat), " ", "")
      rownames(design.mat) <- str_replace_all(rownames(design.mat), "\\.", "")
      colnames(df.exp) <- str_replace_all(colnames(df.exp), " ", "")
      colnames(df.exp) <- str_replace_all(colnames(df.exp), "\\.", "")
      
      df.exp <- df.exp[sort(colnames(df.exp))]
      design.mat <- design.mat[sort(rownames(design.mat)),]
      
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = df.exp, colData = design.mat, design = ~ Metastatic_status + disease_status)

      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep,]
      
      dds$Metastatic_status <- relevel(dds$Metastatic_status, ref = "0")
      
      dds <- DESeq(dds)
      
      save(dds, file = str_glue("{proj}_DE_met.RData"))
      
      
      # write the Rdata file,the csv from results
      
      # look for bone, liver and lung and Lymph Node tropisms in each cancer type
      
  
      # GSEA with transcripts and ranked by fold change
      
      }
    }
  }
  
#dat <- data.table::fread("~/storage/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BRCA_metastatic_data_RNAseq.csv")


projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]
# DE  of selected features

# Up and down regulated pathways in positive tumors
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
      dat <- ceiling(dat)
      # Differential expression for just selected features
      
      # Read expression data and select in normal data
      
      # select row in the row names in the dat file
      
      # DEseq2
      
      # volcano, enrichment for the DE genes
      
      # save the DE genes in selected features
      
      # compare the overlapping up and down regulated pathways
      
      # remove NAs if they exist
      res <- res[complete.cases(res),]
      
      # add a column of NAs
      res$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "DOWN"
      
      # make the 
      res$delabel <- NA
      res$delabel[res$diffexpressed != "NO"] <- res$ENSEMBL[res$diffexpressed != "NO"]
      
      ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
        geom_point() + 
        theme_minimal() +
        geom_text()
      
      
      # Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
      # load library
      library(ggrepel)
      # plot adding up all layers we have seen so far
      p1 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "gray", "red")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      
      p2 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("purple", "black", "orange")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/Volcano/{proj}_{org}"),"_DE_vol_",".pdf"),
             plot = p1, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/Volcano/{proj}_{org}_2"),"_DE_vol_",".pdf"),
             plot = p2, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      
      
      
    }
  }
}





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
