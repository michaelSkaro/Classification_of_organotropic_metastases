# In this investigation we will coduct these analyses:


# Differenital expression between normal tissue and primary tumors 
# Differenital expression between primary tumors that have progressed and those that have not yet. *****DONE***** 
# Differenital expression between seeding loci from one primary location
# Differenital expression between primary tumors arising in different locations but seed in the same location

# Construct an adjacency matrix using WGCNA
# Extract the adjacency matricies. Select matricies with significant of coverage of samples and 
# use these as input matricies for learning

# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)



#projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

# add the clinical annotations of progression

# read in the files in one by one and select the correct columns. 
# concatenate the files and end with one uniform progression file. 

proj <- projects[1]

foo <- prog_file

for(proj in projects){
  
  prog_file <- data.table::fread(str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{proj}_metastatic_status.csv"), 
                                 stringsAsFactors = TRUE, header = TRUE) %>%
    dplyr::select(bcr_patient_barcode,met_loc,LymphNodeStatus,Metastatic_status)

  foo <- rbind(foo, prog_file)
  
}


foo <- foo[order(foo$bcr_patient_barcode),]

foo <- foo[!duplicated(foo$bcr_patient_barcode),]

prog_file <- foo




View(prog_file)

normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$sample_type == "Primary Tumor",]
tumor.samples <- tumor.samples[tumor.samples$project %in% projects,]

tumor.samples$bcr_patient_barcode <- substr(tumor.samples$barcode_short, 0, 12)
clinical$barcode_short <- substr(clinical$barcode, 0,16)




tumor.samples <- left_join(tumor.samples, prog_file, by = "bcr_patient_barcode")



write.csv(tumor.samples, file = str_glue("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv"))


proj <- projects[1]

for (proj in projects) {
  
  df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
    as_tibble() %>%
    tibble::column_to_rownames(var = "Ensembl")
  
  n<-dim(df.exp)[1]
  df.exp<-df.exp[1:(n-5),]
  
  coldata<- tumor.samples[tumor.samples$project == proj,]
  #coldata.n <- normal.samples[normal.samples$project == proj,]
  
  df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
  rownames(coldata) <- coldata$barcode
  coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)
  
  rownames(coldata) <- sort(rownames(coldata))
  colnames(df.exp) <- sort(colnames(df.exp))

  coldata$Metastatic_status <- as.factor(coldata$Metastatic_status)
  
  index <- is.na(coldata$Metastatic_status)
  coldata$Metastatic_status[index] <- 0
  
  
  dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ Metastatic_status)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds$Metastatic_status <- relevel(dds$Metastatic_status, ref = "0")

  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$pvalue),]
  resOrdered <- as.data.frame(resOrdered)
  res <- as.data.frame(res)
  
  write.csv(resOrdered, file = str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/{proj}_DE.csv"))
}

# Learning

# cran doesn't work for this file. You must load it into your working directory which is some what of a pain in the ass

library(RCurl)
## paste URL to make it easier to read code (cosmetic!)
dat_url <- paste0("https://github.com/cran/mlDNA/blob/master/data/mlDNA.rda")
f <- getBinaryURL()
L <- load(rawConnection(f))



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
# 
# 
