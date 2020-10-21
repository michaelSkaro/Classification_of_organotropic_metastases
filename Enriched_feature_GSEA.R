# Enrichment of Ranked features for each seeding locations
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts")


projects <- c("BLCA", "BRCA", "COAD","HNSC","LIHC","LUAD")
proj <- projects[1]


for(proj in projects){
  dat <- data.table::fread(str_glue("{proj}.txt", header = TRUE))
  organs <- names(dat)
  for(org in organs){
    if(org %in% names(dat)){
      df.exp <- dat%>%
        dplyr::select(org)
      names(df.exp) <- "V1"
      df.exp$V1 <- substr(df.exp$V1, 1, 15)
    
    
    Genelist <- clusterProfiler::bitr(df.exp$V1, OrgDb = org.Hs.eg.db ,fromType = "ENSEMBL", toType = "SYMBOL", drop = TRUE)
    system.time(ansEA <- TCGAanalyze_EAcomplete(TFname=str_glue("Highly Enirched {proj} progression to {org}"),Genelist$SYMBOL))
    library(TCGAbiolinks)
    TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                            GOBPTab = ansEA$ResBP,
                            GOCCTab = ansEA$ResCC,
                            GOMFTab = ansEA$ResMF,
                            PathTab = ansEA$ResPat,
                            nRGTab = Genelist, 
                            nBar = 15,
                            filename=paste0("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts/",
                                            proj,"_",org,"_","TCGA_viz",".pdf"))
    
      }
    } 
  }
  
#Fisher's test within Cancer

BiocManager::install("GeneOverlap")
library(GeneOverlap)

BLCA <- data.table::fread("BLCA.txt", header = TRUE)%>%
  dplyr::select(-Rank)
BRCA <- data.table::fread("BRCA.txt", header = TRUE)%>%
  dplyr::select(-Rank)
COAD <- data.table::fread("COAD.txt", header = TRUE)%>%
  dplyr::select(-Rank)
LIHC <- data.table::fread("LIHC.txt", header = TRUE)%>%
  dplyr::select(-Rank)
HNSC <- data.table::fread("HNSC.txt", header = TRUE)%>%
  dplyr::select(-Rank)
LUAD <- data.table::fread("LUAD.txt", header = TRUE)%>%
  dplyr::select(-Rank)

BLCA <- as.data.frame(BLCA)
BRCA <- as.data.frame(BRCA)
COAD <- as.data.frame(COAD)
HNSC <- as.data.frame(HNSC)
LIHC <- as.data.frame(LIHC)
LUAD <- as.data.frame(LUAD)



projects <- c("BLCA", "BRCA", "COAD","HNSC","LIHC","LUAD")
proj <- projects[1]

# within cancer overlap: transcripts driving distant metastasis
for(proj in projects){
  dat <- as.data.frame(data.table::fread(str_glue("{proj}.txt", header = TRUE))) %>%
    dplyr::select(-Rank)
  combos <- combn(ncol(dat),2)
  out <- adply(combos, 2, function(x) {
    go.obj <- newGeneOverlap(dat[,x[1]],dat[,x[2]],genome.size = 60483)
    go.obj <- testGeneOverlap(go.obj)
    data.frame("CancerType"=proj,
                    "listA"=names(dat)[x[1]],
                    "listB"=names(dat)[x[2]],
                    "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                    "intersection"= length(go.obj@intersection),
                    "p.value" = go.obj@pval)
  })
  out <- out %>%
    dplyr::select(-"X1")
  write.csv(out, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/Fisher_exact_test_within_cancer_{proj}.csv"))
} 

# between cancer types comparisons

# within cancer overlap: transcripts driving distant metastasis

f <- function(x,y){
  
  go.obj <- newGeneOverlap(x,y,genome.size = 60483)
  go.obj <- testGeneOverlap(go.obj)
  data.frame( "CancerType1" = deparse(substitute(x)),
              "CancerType2" = deparse(substitute(y)),
              "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
             "intersection"= length(go.obj@intersection),
             "p.value" = go.obj@pval)
}
# BLCA vs BRCA
out<- t(sapply(intersect(colnames(BLCA),colnames(BRCA)), function(x) f(BLCA[,x], BRCA[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BLCA_BRCA_FE.csv")
# BLCA vs COAD
out<- t(sapply(intersect(colnames(BLCA),colnames(COAD)), function(x) f(BLCA[,x], COAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BLCA_COAD_FE.csv")
# BLCA vs HNSC
out<- t(sapply(intersect(colnames(BLCA),colnames(HNSC)), function(x) f(BLCA[,x], HNSC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BLCA_HNSC_FE.csv")
# BLCA vs LIHC
out<- t(sapply(intersect(colnames(BLCA),colnames(LIHC)), function(x) f(BLCA[,x], LIHC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BLCA_LIHC_FE.csv")
# BLCA VS LUAD
out<- t(sapply(intersect(colnames(BLCA),colnames(LUAD)), function(x) f(BLCA[,x], LUAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BLCA_LUAD_FE.csv")
# BRCA vs COAD
out<- t(sapply(intersect(colnames(BRCA),colnames(COAD)), function(x) f(BRCA[,x], COAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BRCA_COAD_FE.csv")
# BRCA vs HNSC
out<- t(sapply(intersect(colnames(BRCA),colnames(HNSC)), function(x) f(BRCA[,x], HNSC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BRCA_HNSC_FE.csv")
# BRCA vs LIHC
out<- t(sapply(intersect(colnames(BRCA),colnames(LIHC)), function(x) f(BRCA[,x], LIHC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BRCA_LIHC_FE.csv")
# BRCA vs LUAD
out<- t(sapply(intersect(colnames(BRCA),colnames(LUAD)), function(x) f(BRCA[,x], LUAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/BRCA_LUAD_FE.csv")
# COAD vs HNSC
out<- t(sapply(intersect(colnames(COAD),colnames(HNSC)), function(x) f(COAD[,x], HNSC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/COAD_HNSC_FE.csv")
# COAD vs LIHC
out<- t(sapply(intersect(colnames(COAD),colnames(LIHC)), function(x) f(COAD[,x], LIHC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/COAD_LIHC_FE.csv")
# COAD vs LUAD
out<- t(sapply(intersect(colnames(COAD),colnames(LUAD)), function(x) f(COAD[,x], LUAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/COAD_LUAD_FE.csv")
# HNSC vs LIHC
out<- t(sapply(intersect(colnames(HNSC),colnames(LIHC)), function(x) f(HNSC[,x], LIHC[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/HNSC_LIHC_FE.csv")
# HNSC vs LUAD 
out<- t(sapply(intersect(colnames(HNSC),colnames(LUAD)), function(x) f(HNSC[,x], LUAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/HNSC_LUAD_FE.csv")
# LIHC vs LUAD
out<- t(sapply(intersect(colnames(HNSC),colnames(LUAD)), function(x) f(HNSC[,x], LUAD[,x])))
write.csv(out, file = "/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts_10_20_20/LIHC_LUAD_FE.csv")






listInput <- list(BRCA_Bone = na.omit(BRCA$Bone),
                BLCA_Bone = na.omit(BLCA$Bone))
upset(fromList(listInput), order.by = "freq",nsets = 2)



# lets make some heatmaps

library(ComplexHeatmap)

# extract the gene sets that are important for each of the data sets
met_samples <- data.table::fread("/mnt/storage/mskaro1/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
met_samples <- met_samples[,2:13]

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[clinical$sample_type == "Solid Tissue Normal",]
tumor.samples <- clinical[clinical$sample_type != "Solid Tissue Normal",]


# objective: pull the samples from the files in the CSBL shared folder and
# combine them into a file for ML classification

projects <- unique(met_samples$project)
projects <- sort(projects)
proj <- projects[1]
organs <- tail(colnames(patients), n=7)
for(proj in projects){
  # read in the genes list table
  p <- substring(proj,6,9)
  genes <- data.table::fread(paste0("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts/",p,".txt"))
  
  # read in the patients
  patients <- data.table::fread(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_metastatic_data_RNAseq.csv"), header = TRUE) %>% 
    tibble::column_to_rownames("V1")
  organs <- tail(colnames(patients), n=7)
  labels <- patients %>%
    dplyr::select(barcode,organs)
  # sub set the rows based on the genes in enriched transcripts
  patients <- dplyr::select(-c(barcode,organs))
  
  # sub set the columns for patients seeding in each location
  
  # plot heatmap, with the rug for each location
  
}

# Lets do some fishers exact tests for the data sets






# ranked analysis then clustering

BiocManager::install("RRHO")

setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts/Features_sorted_for_each_class")

file_list <- list.files()
met_samples <- data.table::fread("/mnt/storage/mskaro1/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
met_samples <- met_samples[,2:13]
projects <- unique(met_samples$project)
projects <- sort(projects)
proj <- projects[1]
fil <- 
for(proj in projects){
  dat <- data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_important_features.csv"))
  View(dat)
  
  if(proj =="TCGA-BLCA"){
    Bone <- as.data.frame(dat) %>%
      dplyr::select("Bone","Bone Importance Score")
    # give a rank
    Bone$rank <- NA
    order.scores<-order(Bone$`Bone Importance Score`,Bone$Bone, decreasing = TRUE)
    Bone$rank[order.scores] <- 1:nrow(dat)
    
    # sort by score importance
    Bone <- Bone[order(Bone$`Bone Importance Score`, decreasing = TRUE),]
    
    
      
  }
  
  
}



# library(RRHO)
# BLCA_bone <- unique(as.data.frame(BLCA$Bone))
# colnames(BLCA_bone) <- "Bone"
# BRCA_bone <- unique(as.data.frame(BRCA$Bone))
# colnames(BRCA_bone) <- "Bone"
# 
# RRHO.Bone <- RRHO(BRCA_bone$Bone[1:250], BLCA_bone$Bone[1:250], BY =TRUE, alternative='enrichment')
# 
# 
# list.length <- 100
# list.names <- paste('Gene',1:list.length, sep='')
# gene.list1<- data.frame(list.names, sample(100))
# gene.list2<- data.frame(list.names, sample(100))
# #Compute overlap and significance
# RRHO.example <- RRHO(gene.list1, gene.list2,BY=TRUE, alternative='enrichment')



BiocManager::install("GeneOverlap")
library(GeneOverlap)



BLCA <- data.table::fread("BLCA.txt", header = TRUE)%>%
  dplyr::select(-Rank)
BRCA <- data.table::fread("BRCA.txt", header = TRUE)%>%
  dplyr::select(-Rank)
COAD <- data.table::fread("COAD.txt", header = TRUE)%>%
  dplyr::select(-Rank)
LIHC <- data.table::fread("LIHC.txt", header = TRUE)%>%
  dplyr::select(-Rank)
HNSC <- data.table::fread("HNSC.txt", header = TRUE)%>%
  dplyr::select(-Rank)

