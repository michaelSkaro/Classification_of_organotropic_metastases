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




