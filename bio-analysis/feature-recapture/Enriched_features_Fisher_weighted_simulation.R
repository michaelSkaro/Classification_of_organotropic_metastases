# Enrichment of Ranked features for each seeding locations
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(dplyr)
library(tidyverse)

BiocManager::install("GeneOverlap")
library(GeneOverlap)

go_overlap <- function(x,y,background){
  
  go.obj <- newGeneOverlap(x,y,genome.size = background)
  go.obj <- testGeneOverlap(go.obj)
  data.frame( "CancerType1" = deparse(substitute(x)),
              "CancerType2" = deparse(substitute(y)),
              "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
             "intersection"= length(go.obj@intersection),
             "p.value" = go.obj@pval)
}


#BiocManager::install("GeneOverlap")


# make map fo weighted values to stratify GO BPs

# read in RNAseq data
df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
  as_tibble() %>%
  tibble::column_to_rownames(var = "Ensembl")
n<-dim(df.exp)[1]
df.exp<-df.exp[1:(n-5),]

# map the values of the IDs to GO IDs

goIDmap <- clusterProfiler::bitr(geneID = substring(rownames(df.exp),1,15), fromType = "ENSEMBL", toType = "GO", OrgDb = org.Hs.eg.db) %>%
  filter(ONTOLOGY == "BP")

# do not over penalize multi-mapping due to multiple values of evidence. Filter for first occurence of ENSEMBL:GOID matches
# cat the ENSEMBLID-GO:ID values, filter for first occurrence

goIDmap$cat <- paste0(goIDmap$ENSEMBL, "-", goIDmap$GO)

# filter filter for the first occuence in the dataframe

goIDmap <- goIDmap[!duplicated(goIDmap$cat),]
# add a count (by x weight) matrix value
weight <- as.data.frame(table(goIDmap$ENSEMBL))
colnames(weight) <- c("ENSEMBL","count")

goIDmap <- left_join(goIDmap, weight, by = "ENSEMBL")

# stratify the score by occurence rate distribution

plot(density(goIDmap$count))

# calculate the mean

mean_mapping_occurence <- mean(goIDmap$count)
# mean = ~24


# normalize the (by x weight) to the mean to make normalized 
# weight scoring mechanism as mean/count

goIDmap$weight <- mean_mapping_occurence/goIDmap$count
plot(density(goIDmap$weight))
mapvar <- var(goIDmap$weight)
mapsd <- sd(goIDmap$weight)

# simulate the feature selection

background_features <- 60483
background_BP <- 23393


# make databank for results

out_BP <- as.data.frame(t(c("SimulatedL1"= NA,
                            "SimulatedL2"= NA,
                            "odds.ratio" = NA,
                            "simulated.intersection"= NA,
                            "simulated_p.value" = NA))) 

for(i in 1:50000){
  
  
  # extract 1000 random features
  features1 <- sample(goIDmap$ENSEMBL, 1000)
  
  features2 <- sample(goIDmap$ENSEMBL, 1000)
  
  ego1 <- clusterProfiler::enrichGO(gene  = substring(features1,1,15),
                                    OrgDb         = "org.Hs.eg.db",
                                    keyType       = 'ENSEMBL',
                                    ont           = "BP",
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.01,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
  
 
  ego2 <- clusterProfiler::enrichGO(gene  = substring(features2,1,15),
                                     OrgDb         = "org.Hs.eg.db",
                                     keyType       = 'ENSEMBL',
                                     ont           = "BP",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01,
                                     qvalueCutoff  = 0.05,
                                     readable      = TRUE)
  # overlap the two sets, filter nonsignificant results
  
  res1 <- ego1@result
  res1 <- res1[res1$p.adjust<0.05,]
  res1 <- res1[res1$pvalue<0.05,]
  
  res2 <- ego2@result
  res2 <- res2[res2$p.adjust<0.05,]
  res2 <- res2[res2$pvalue<0.05,]
  
  # penalize high chance matches from the two sets, filter for matches > 1sd over mean
  
  res1 <- left_join(res1, goIDmap, by = "GO") %>%
    filter(weight > (mean_mapping_occurence + mapsd))
  
  res2 <- left_join(res1, goIDmap, by = "GO") %>%
    filter(weight > (mean_mapping_occurence + mapsd))
  
  go.obj <- newGeneOverlap(L1, L2, genome.size = background_BP)
  go.obj <- testGeneOverlap(go.obj)
  df <- as.data.frame(t(c("SimulatedL1"= length(L1),
                          "SimulatedL2"= length(L2),
                          "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                          "simulated.intersection"= length(go.obj@intersection),
                          "simulated_p.value" = go.obj@pval)))
  out_BP <- rbind(out_BP,df)
  
}

write.csv(out_BP,str_glue("Simulated_random_feature_selection_and_GO_BP.csv"))
rm(out_BP)



background_features <- 60483
background_BP <- 23393


# make databank for results
out_F <- as.data.frame(t(c("SimulatedL1"= NA,
                           "SimulatedL2"= NA,
                           "odds.ratio" = NA,
                           "simulated.intersection"= NA,
                           "simulated_p.value" = NA))) 

out_BP <- as.data.frame(t(c("SimulatedL1"= NA,
                            "SimulatedL2"= NA,
                            "odds.ratio" = NA,
                            "simulated.intersection"= NA,
                            "simulated_p.value" = NA))) 


# simulate the sampling
library(gtools)
comp <- combinations(n = 5, r = 2, v = c(100,200,300,400,500), repeats.allowed = TRUE)

for(i in c(1:15)){
  compL1 <- comp[i,1]
  compL2 <- comp[i,2]
  
  out_F <- as.data.frame(t(c("SimulatedL1"= NA,
                             "SimulatedL2"= NA,
                             "odds.ratio" = NA,
                             "simulated.intersection"= NA,
                             "simulated_p.value" = NA))) 
  
  out_BP <- as.data.frame(t(c("SimulatedL1"= NA,
                              "SimulatedL2"= NA,
                              "odds.ratio" = NA,
                              "simulated.intersection"= NA,
                              "simulated_p.value" = NA)))
  for(j in 1:50000){
    
    # simulate 100 length list
    L1 <- sample(x = c(1:background_features),size = comp[i,1],replace = FALSE)
    L2 <- sample(x = c(1:background_features),size = comp[i,2],replace = FALSE)
    
    go.obj <- newGeneOverlap(L1, L2, genome.size = background_features)
    go.obj <- testGeneOverlap(go.obj)
    df <- as.data.frame(t(c("SimulatedL1"= length(L1),
                            "SimulatedL2"= length(L2),
                            "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                            "simulated.intersection"= length(go.obj@intersection),
                            "simulated_p.value" = go.obj@pval)))
    out_F <- rbind(out_F,df)
    
    L1 <- sample(x = c(1:background_BP),size = comp[i,1],replace = FALSE)
    L2 <- sample(x = c(1:background_BP),size = comp[i,2],replace = FALSE)
    
    go.obj <- newGeneOverlap(L1, L2, genome.size = background_BP)
    go.obj <- testGeneOverlap(go.obj)
    df <- as.data.frame(t(c("SimulatedL1"= length(L1),
                            "SimulatedL2"= length(L2),
                            "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                            "simulated.intersection"= length(go.obj@intersection),
                            "simulated_p.value" = go.obj@pval)))
    out_BP <- rbind(out_BP,df)
  }
  
  write.csv(out_F,str_glue("Simulated_{compL1}_{compL2}_Features.csv"))
  write.csv(out_BP,str_glue("Simulated_{compL1}_{compL2}_BP.csv"))
  #rm(out_BP)
  #rm(out_F)
}

out_BP <- out_BP[2:750001,]
out_BP <- out_BP[order(out_BP$simulated.intersection),]
write.csv(out_BP, "simulated_results.csv")

# visualize the weighted simulation analysis
ggplot(data=out_BP, aes(x=1/(simulated.intersection), y=simulated_p.value, group=comparison, color=comparison)) +
  xlab("Simulated Intersection") +
  ylab("Simulated pvalue") +
  geom_line() +
  theme_classic()

















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
