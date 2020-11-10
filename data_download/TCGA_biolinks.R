library(TCGAbiolinks)
library(stringr)
library(tidyverse)

#get the counts data
projects <- c(
  c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA"),
  c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")
)

  # get clinical 
for (proj in projects){
  
  query <- GDCquery(project = proj,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  
  GDCdownload(query)
  
  
  #save file for manipulation
  data <- GDCprepare(query, save = TRUE, 
                     save.filename = str_glue("{proj}.rda"),
                     remove.files.prepared = TRUE)
  }
  

#Clin data

library(TCGAbiolinks)
library(xml2)
library(tidyverse)

setwd("~/storage/data/TCGA/clinical/TCGAbiolinks")

# Download xml clinical data files -------------------------------------------
projects <- c(
  c("TCGA-BRCA", "TCGA-COAD", "TCGA-HNSC", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LUAD", "TCGA-STAD", "TCGA-THCA"),
  c("TCGA-BLCA", "TCGA-ESCA", "TCGA-KICH", "TCGA-LIHC", "TCGA-LUSC")
)
for (proj in projects) {
  query <- GDCquery(project = proj,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
}

# Parse xml files for each cancer type ---------------------------------------
projects <- list.files("./GDCdata/")

for (proj in projects.downloaded) {
  patients <- list.files(
    stringr::str_glue("./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/"),
    pattern = "xml$", recursive = T)
  filenames <- stringr::str_glue(
    "./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/{patients}"
  )
  patients <- map_chr(str_split(patients, "/"), function(x) x[[1]])
  
  res <- NULL
  for (i in seq_along(filenames)) {
    df <- read_xml(filenames[i])
    barcode <- df %>%
      xml_find_first("//shared:bcr_patient_barcode") %>%
      xml_text()
    
    # In project:patient, look for all fields with text
    df <- xml_find_all(df, "//*[local-name(.) = 'patient']//*[text()]")
    df <- tibble(
      fileID = patients[i],
      field = map_chr(df, xml_path),
      value = map_chr(df, xml_text)
    ) %>%
      mutate(field = str_extract(field, "[^:]+$"))
    res <- bind_rows(res, df)
  }
  fields <- c("fileID", unique(res$field))
  res <- res %>%
    distinct() %>%
    group_by(fileID, field) %>%
    summarise(value = paste(value, collapse = "|")) %>%
    ungroup() %>%
    spread(key = field, value = value) %>%
    select(one_of(fields))
  data.table::fwrite(res, str_glue("./{proj}.csv"))
}

seession.info()

                      # R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] plyr_1.8.6                  maftools_2.4.12             TCGAbiolinks_2.16.4         EnhancedVolcano_1.6.0       ggrepel_0.8.2              
# [6] org.Hs.eg.db_3.11.4         AnnotationDbi_1.50.3        clusterProfiler_3.16.1      DESeq2_1.28.1               SummarizedExperiment_1.18.2
# [11] DelayedArray_0.14.1         matrixStats_0.57.0          Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2        
# [16] IRanges_2.22.2              S4Vectors_0.26.1            BiocGenerics_0.34.0         forcats_0.5.0               dplyr_1.0.2                
# [21] purrr_0.3.4                 readr_1.4.0                 tidyr_1.1.2                 ggplot2_3.3.2               tidyverse_1.3.0            
# [26] tibble_3.0.4                stringr_1.4.0               RRHO_1.28.0                 UpSetR_1.4.0                GeneOverlap_1.24.0         
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1           backports_1.1.10       fastmatch_1.1-0        BiocFileCache_1.12.1   igraph_1.2.6           splines_4.0.2         
# [7] BiocParallel_1.22.0    urltools_1.7.3         digest_0.6.25          htmltools_0.5.0        GOSemSim_2.14.2        viridis_0.5.1         
# [13] GO.db_3.11.4           fansi_0.4.1            magrittr_1.5           memoise_1.1.0          annotate_1.66.0        graphlayouts_0.7.0    
# [19] modelr_0.1.8           R.utils_2.10.1         askpass_1.1            enrichplot_1.8.1       prettyunits_1.1.1      colorspace_1.4-1      
# [25] rappdirs_0.3.1         blob_1.2.1             rvest_0.3.6            xfun_0.18              haven_2.3.1            crayon_1.3.4          
# [31] RCurl_1.98-1.2         jsonlite_1.7.1         scatterpie_0.1.5       genefilter_1.70.0      survival_3.2-7         glue_1.4.2            
# [37] polyclip_1.10-0        gtable_0.3.0           zlibbioc_1.34.0        XVector_0.28.0         scales_1.1.1           DOSE_3.14.0           
# [43] futile.options_1.0.1   DBI_1.1.0              Rcpp_1.0.5             viridisLite_0.3.0      xtable_1.8-4           progress_1.2.2        
# [49] gridGraphics_0.5-0     bit_4.0.4              europepmc_0.4          httr_1.4.2             fgsea_1.14.0           gplots_3.1.0          
# [55] RColorBrewer_1.1-2     ellipsis_0.3.1         R.methodsS3_1.8.1      pkgconfig_2.0.3        XML_3.99-0.5           farver_2.0.3          
# [61] dbplyr_1.4.4           locfit_1.5-9.4         ggplotify_0.0.5        tidyselect_1.1.0       rlang_0.4.8            reshape2_1.4.4        
# [67] munsell_0.5.0          cellranger_1.1.0       tools_4.0.2            downloader_0.4         cli_2.1.0              generics_0.0.2        
# [73] RSQLite_2.2.1          broom_0.7.1            ggridges_0.5.2         evaluate_0.14          yaml_2.2.1             knitr_1.30            
# [79] bit64_4.0.5            fs_1.5.0               tidygraph_1.2.0        caTools_1.18.0         ggraph_2.0.3           formatR_1.7           
# [85] R.oo_1.24.0            DO.db_2.9              xml2_1.3.2             biomaRt_2.44.4         compiler_4.0.2         rstudioapi_0.11       
# [91] curl_4.3               reprex_0.3.0           tweenr_1.0.1           geneplotter_1.66.0     stringi_1.5.3          futile.logger_1.4.3   
# [97] lattice_0.20-41        Matrix_1.2-18          vctrs_0.3.4            pillar_1.4.6           lifecycle_0.2.0        BiocManager_1.30.10   
# [103] triebeard_0.3.0        data.table_1.13.0      cowplot_1.1.0          bitops_1.0-6           qvalue_2.20.0          R6_2.4.1              
# [109] KernSmooth_2.23-17     gridExtra_2.3          lambda.r_1.2.4         MASS_7.3-53            gtools_3.8.2           assertthat_0.2.1      
# [115] openssl_1.4.3          withr_2.3.0            GenomeInfoDbData_1.2.3 hms_0.5.3              VennDiagram_1.6.20     rmarkdown_2.4         
# [121] rvcheck_0.1.8          ggforce_0.3.2          lubridate_1.7.9        tinytex_0.26
# 

                      
                      
