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
  
# lets make some venn diagrams

# Load library
library(VennDiagram)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)

setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts")
projects <- c("BLCA", "BRCA", "COAD","HNSC","LIHC","LUAD")
proj <- "BLCA"


dat1 <- data.table::fread(str_glue("BLCA.txt", header = TRUE))
organs1 <- names(dat1)
proj2 <- "BRCA"
dat2 <- data.table::fread(str_glue("BRCA.txt", header = TRUE))

dat1 <- dat1 %>%
  dplyr::select(names(dat2))

BLCA_Bone <- as.list(na.omit(dat1$Bone))
BLCA_Liver <- as.list(na.omit(dat1$Liver))
BLCA_Lung <- as.list(na.omit(dat1$Lung))

BRCA_Bone <- as.list(na.omit(dat2$Bone))
BRCA_Liver <- as.list(na.omit(dat2$Liver))
BRCA_Lung <- as.list(na.omit(dat2$Lung))

# Load library
library(VennDiagram)

# Generate 3 sets of 200 words
BLCA_Bone
BLCA_Liver
BLCA_Lung
BRCA_Bone
BRCA_Liver
BRCA_Lung

#package VennDiagram
install.packages("VennDiagram") # only necessary 1st time
library(VennDiagram)
overlap <- intersect(BLCA_Bone, BRCA_Bone)  # easy/useful - 625 proteins
#draw the Venn diagram using the venn.plot() function
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(BLCA_Bone),#no MCP set
                                area2 = length(BRCA_Bone),#no JPR set
                                cross.area = length(overlap),
                                c("BLCA_Bone", "BRCA_Bone"), scaled = TRUE,
                                fill = c("green", "blue"),
                                cex = 1.5,
                                cat.cex = 1.5,
                                cat.pos = c(320, 25),
                                cat.dist = .05) 


overlap <- intersect(BLCA_Liver, BRCA_Liver)  # easy/useful - 625 proteins
#draw the Venn diagram using the venn.plot() function
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(BLCA_Liver),#no MCP set
                                area2 = length(BRCA_Liver),#no JPR set
                                cross.area = length(overlap),
                                c("BLCA_Liver", "BRCA_Liver"), scaled = TRUE,
                                fill = c("green", "blue"),
                                cex = 1.5,
                                cat.cex = 1.5,
                                cat.pos = c(320, 25),
                                cat.dist = .05) 

overlap <- intersect(BLCA_Lung, BRCA_Lung)  # easy/useful - 625 proteins
#draw the Venn diagram using the venn.plot() function
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = length(BLCA_Lung),#no MCP set
                                area2 = length(BRCA_Lung),#no JPR set
                                cross.area = length(overlap),
                                c("BLCA_Lung", "BRCA_Lung"), scaled = TRUE,
                                fill = c("green", "blue"),
                                cex = 1.5,
                                cat.cex = 1.5,
                                cat.pos = c(320, 25),
                                cat.dist = .05) 

#OR Upset may work best

listInput <- list(BRCA_Bone = na.omit(dat2$Bone),BRCA_Lung = na.omit(dat2$Lung),BRCA_Liver = na.omit(dat2$Liver),
                  BLCA_Bone = na.omit(dat1$Bone), BLCA_Lung = na.omit(dat1$Lung), 
                  BLCA_Liver = na.omit(dat1$Liver))
upset(fromList(listInput), order.by = "freq",nsets = 6)

listInput <- list(BRCA_Bone = na.omit(dat2$Bone),
                  BLCA_Bone = na.omit(dat1$Bone))
upset(fromList(listInput), order.by = "freq",nsets = 2)


setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/ranked_features/Ranked_transcripts")
projects <- c("BLCA", "BRCA", "COAD","HNSC","LIHC","LUAD")

dat1 <- data.table::fread("BLCA.txt", header = TRUE)
dat2 <- data.table::fread("COAD.txt", header = TRUE)

listInput <- list(BLCA_Prostate = na.omit(dat1$Prostate),COAD_Postate = na.omit(dat2$Prostate))
upset(fromList(listInput), order.by = "freq",nsets = 2)

dat1 <- data.table::fread("BLCA.txt", header = TRUE)
dat2 <- data.table::fread("BRCA.txt", header = TRUE)
dat3 <- data.table::fread("LIHC.txt", header = TRUE)

listInput <- list(BLCA_Liver = na.omit(dat1$Liver),BRCA_Liver = na.omit(dat2$Liver), 
                  LIHC_Liver = na.omit(dat3$Liver))
upset(fromList(listInput), order.by = "freq",nsets = 3)

dat1 <- data.table::fread("BLCA.txt", header = TRUE)
dat2 <- data.table::fread("BRCA.txt", header = TRUE)
dat3 <- data.table::fread("HNSC.txt", header = TRUE)
dat4 <- data.table::fread("LIHC.txt", header = TRUE)
dat5 <- data.table::fread("LUAD.txt", header = TRUE)

listInput <- list(BLCA_Lung = na.omit(dat1$Lung), BRCA_Lung = na.omit(dat2$Lung), 
HNSC_Lung = na.omit(dat3$Lung), LIHC_Lung = na.omit(dat4$Lung), 
LUAD_Lung = na.omit(dat5$Lung))
upset(fromList(listInput), order.by = "freq")


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

# pull the correct case ids for miRNAseq learning
# pull the correct case ids somatic variation learning







# plot the scores, plot the heatmaps for each of the genes

#drafts
set.seed(123)
mat1 = matrix(rnorm(80, 2), 8, 10)
mat1 = rbind(mat1, matrix(rnorm(40, -2), 4, 10))
rownames(mat1) = paste0("R", 1:12)
colnames(mat1) = paste0("C", 1:10)

mat2 = matrix(runif(60, max = 3, min = 1), 6, 10)
mat2 = rbind(mat2, matrix(runif(60, max = 2, min = 0), 6, 10))
rownames(mat2) = paste0("R", 1:12)
colnames(mat2) = paste0("C", 1:10)

le = sample(letters[1:3], 12, replace = TRUE)
names(le) = paste0("R", 1:12)

ind = sample(12, 12)
mat1 = mat1[ind, ]
mat2 = mat2[ind, ]
le = le[ind]

ht1 = Heatmap(mat1, name = "rnorm")
ht2 = Heatmap(mat2, name = "runif")
ht3 = Heatmap(le, name = "letters")

