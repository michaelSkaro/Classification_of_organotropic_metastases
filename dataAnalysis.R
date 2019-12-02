# data analysis of 61 GPLs on 124 experiments

# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Dec 2 16:35:54 EST 2019

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO GSE10893, GPL1390

gset <- getGEO("GSE10893", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("03000000000030000000300000022000000330000000003000",
               "20000000X0000000000000000000040XXXXXXXXX0000000000",
               "00053500000000000000000000000000000000000000006000",
               "00000000000000000000000002220022222")
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
cont.matrix <- makeContrasts(G6-G0, G2-G0, G3-G2, G4-G3, G5-G4, G6-G5, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = 22576)

setwd("~/storage/Metastatic_Organo_Tropism/results")
write.csv(tT,"GSE10893_GPL1390.csv")
write.csv(cont.matrix,"Design_Matrix_GSE10893_GPL1390.csv")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE10893", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("03000000000030000000300000022000000330000000003000",
               "20000000X0000000000000000000040XXXXXXXXX0000000000",
               "00053500000000000000000000000000000000000000006000",
               "00000000000000000000000002220022222")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  set group names

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Primary+Breast+Cancer","Brain+Metastasis+from+Breast+Tumor","Healthy+Control","Primary+Breast+tumor+with+Brain+Metastasis","Brain+Metastasis","Breast+Tumor+With+Lymph+Node+Metastasis","Skin+metastasis+from+Breast+Primary+tumor")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5","#dff4e4","#f4dff4","#c7ff9d", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE10893", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")





