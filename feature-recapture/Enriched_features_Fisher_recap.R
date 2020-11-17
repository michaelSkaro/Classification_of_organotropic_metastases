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


