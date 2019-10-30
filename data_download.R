setwd("~/storage/Metastatic_Organo_Tropism")
#Metastatic Disease datasheet

refDat <- data.table::fread("Metastatic_database_project_information.csv")

# view the appropriate headers we need to iterate over

colnames(refDat)

# [1] "Experiment_id"     "Pubmed_id"         "Dataset_id"        "Platform_id"       "Standard"          "Class_id"          "Class_name"        "Sample_id"        
# [9] "Cancer_type"       "Cancer_subtype"    "Metastasis_status" "Primary_site"      "Metastasis_site"   "Sample_label" 

# use the experiment id in our i for loop, then use the datasetID and the platformID for j and k

exps <- unique(refDat$Experiment_id)
gseIDs <- unique(refDat$Dataset_id)
GSMs <- unique(refDat$Sample_id)
platforms <- unique(refDat$Platform_id)

plat <- platforms[1]
for(plat in platforms){
  coldata <- refDat %>%
    filter(Platform_id == plat)
  GSMs <- unique(coldata$Sample_id)
  dat <- NULL
  for(genesSets in GSMs){
    ##download gsms and concatenate
    genesSets <- GSMs[1] # for practice
    # download samples
    # gds <- getGEO("GSM11805")
    gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE, GSEMatrix = TRUE)
    # Metadata
    View(head(Meta(gsm)))
    # Values
    Table(gsm)[1:5,]
    # Look at Column descriptions:
    tibble(Columns(gsm))
    
    
    View(Table(gsm))
    View(Columns(gsm))
    
    eset <- GDS2eSet(gsm,do.log2=TRUE)
    var <- as.data.frame(eset@assayData$exprs) %>%
      tibble::rownames_to_column("Probes")
    dat <- var
    annotGPL <- pData(eset)
  }
  
}

# get each model GSE
##################################################################
#GSE10893
#GPL1390

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE10893", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("03000000000030000000300000011000000330000000004000",
               "10000000000000000000000000000302222222220000000000",
               "0004440000000000000000000000000000000000000000X100",
               "00000000000000000000000001110000011")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G4-G0, G1-G0, G2-G1, G3-G2, G4-G3, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE10893", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("03000000000030000000300000011000000330000000004000",
               "10000000000000000000000000000302222222220000000000",
               "0004440000000000000000000000000000000000000000X100",
               "00000000000000000000000001110000011")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
#set group names
sml <- paste("G", sml, sep="")  

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Breast+Cancer","Normal+Breast","Treated+Breast+Cancer","Breast+Cancer+Brain+Met","Breast+Cancer+Lymph+Met")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4","#dfeaf4","#f4dfdf","#f2cb98", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE10893", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")








# get the libraries
library(Biobase)
library(oligoClasses)
library(knitr)
library(BiocStyle)
library(oligo)
library(geneplotter)
library(ggplot2)
library(dplyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(ArrayExpress)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(topGO)
library(genefilter)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pheatmap)
library(mvtnorm)
library(DAAG)
library(multcomp)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(openxlsx)
library(devtools)
library(biomaRt)
library(ArrayExpress)


