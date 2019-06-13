# clear environment
rm(list=ls())

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
#library(EnrichmentBrowser)
set.seed(777)
raw_data_dir <- file.path(getwd(), "rawDataMAWorkflow")



# Make project directory

if(!dir.exists(raw_data_dir)){
  dir.create(raw_data_dir)
}

# get the data we need
anno_AE <- getAE("E-MTAB-2967", path=raw_data_dir, type="raw")

SDRF <- read.delim(
  url("http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2967/E-MTAB-2967.sdrf.txt"))

# annotate the data
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)


raw_data <- read.celfiles(as.character(SDRF$Array.Data.File), verbose = FALSE, phenoData = SDRF)

# Check I read in values correctly 
validObject(raw_data)

head(pData(raw_data))
head(exprs(raw_data))

# make the design matrix to split the groups. This step 
# can also be completed with the gather and spread commands
# in tidyr
pData(raw_data) <- pData(raw_data)[, c("Source.Name",
                                       "Characteristics.individual.",
                                       "Factor.Value.disease.",
                                       "Factor.Value.phenotype.")]


# QC

exp_raw <- log2(exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale = FALSE)

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = pData(raw_data)$Factor.Value.disease.,
                     Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                     Individual = pData(raw_data)$Characteristics.individual.)

# plot prinicipal component to show differences in data are biological and not from source of data
(qplot(PC1, PC2, data = dataGG, color = Disease,
       main = "PCA plot of the raw data (log-transformed)", size = I(2),
       asp = 1.0, geom = "text",
       label = Individual)
  + scale_colour_brewer(palette = "Set2"))


# shw the box plot of the data to control for outliers
boxplot(raw_data, target = "core",
        main = "Boxplots of log2-intensities for the raw data")

# this has non-zero, cannot seem to install correctly, skip this step
# arrayQualityMetrics(expressionset = raw_data,
#                     outdir = "Report_for_Palmieri_raw",
#                     force = TRUE, do.logtransform = TRUE,
#                     intgroup = c("Factor.Value.disease." , "Factor.Value.phenotype."))



head(ls("package:hugene10sttranscriptcluster.db"))


# Happy background normalization and RMA scaling :)

palmieri_eset <- oligo::rma(raw_data, target="core")


# now we can use the normalized data to mae a new PCA and then a happy heatmap

exp_palmieri <- exprs(palmieri_eset)
PCA <- prcomp(t(exp_palmieri), scale = FALSE)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease = pData(palmieri_eset)$Factor.Value.disease.,
                     Phenotype = pData(palmieri_eset)$Factor.Value.phenotype.)
(qplot(PC1, PC2, data = dataGG, color = Disease, shape = Phenotype,
       main = "PCA plot of the calibrated data", size = I(2), asp = 1.0)
  + scale_colour_brewer(palette = "Set2"))



dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))
colnames(dists) <- NULL
diag(dists) <- NA
rownames(dists) <- pData(palmieri_eset)$Factor.Value.phenotype.
hmcol <- colorRampPalette(rev(brewer.pal(9, "PuOr")))(255)
pheatmap(dists, col = rev(hmcol), clustering_distance_rows = "manhattan",
         clustering_distance_cols = "manhattan")




# filter genes

no_of_samples <- table(paste0(pData(palmieri_eset)$Factor.Value.disease., "_",
                              pData(palmieri_eset)$Factor.Value.phenotype.))

palmieri_medians <- rowMedians(exprs(palmieri_eset))
hist_res <- hist(palmieri_medians, 100, col="#e7efd8", freq = FALSE,
                 main = "Histogram of the median intensities",
                 xlab = "Median intensities")
emp_mu <- hist_res$breaks[which.max(hist_res$density)]
emp_sd <- mad(palmieri_medians)/2
prop_cental <- 0.50
lines(sort(palmieri_medians), prop_cental*dnorm(sort(palmieri_medians),
                                                mean = emp_mu , sd = emp_sd),col = "grey10", lwd = 4)


cut_val <- 0.05 / prop_cental
thresh_median <- qnorm(0.05 / prop_cental, emp_mu, emp_sd)

# cutt samples under threshold
samples_cutoff <- min(no_of_samples)

idx_thresh_median <- apply(exprs(palmieri_eset), 1, function(x){
  sum(x > thresh_median) >= samples_cutoff})

# subset based on true values

palmieri_filtered <- subset(palmieri_eset, idx_thresh_median)


# annotate your transcripts

anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys=(featureNames(palmieri_filtered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype="PROBEID")
probe_stats <- anno_palmieri %>%
  group_by(PROBEID) %>%
  summarize(no_of_matches = n_distinct(SYMBOL)) %>%
  filter(no_of_matches > 1)


# exclude some ids that are ambiguous

ids_to_exlude <- ((featureNames(palmieri_filtered) %in% probe_stats$PROBEID) |
                    featureNames(palmieri_filtered) %in% subset(anno_palmieri ,
                                                                is.na(SYMBOL))$PROBEID)
table(ids_to_exlude)
palmieri_final <- subset(palmieri_filtered, !ids_to_exlude)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))
fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)


rownames(fData(palmieri_final)) <-fData(palmieri_final)$PROBEID







# design matrix for DE 

individual <- as.character(pData(palmieri_final)$Characteristics.individual.)

# get ride of the spaces
tissue <- str_replace_all(pData(palmieri_final)$Factor.Value.phenotype., " ", "_")

# flip the tissues from non inflamed to NI if its not NI then its an I
tissue <- ifelse(tissue == "non-inflamed_colonic_mucosa", "nI", "I")


disease <- str_replace_all(pData(palmieri_final)$Factor.Value.disease., " ", "_")
disease <- ifelse(disease == "Crohn's_disease", "CD", "UC")


i <- individual[disease == "CD"]
design_palmieri_CD <- model.matrix(~ 0 + tissue[disease == "CD"] + i)
colnames(design_palmieri_CD)[1:2] <- c("I", "nI")

i <- individual[disease == "UC"]
design_palmieri_UC<- model.matrix(~ 0 + tissue[disease == "UC"] + i)
colnames(design_palmieri_UC)[1:2] <- c("I", "nI")




# Bayes model to fit the disease and tissue: Differential expression

contrast_matrix_CD <- makeContrasts(I-nI, levels = design_palmieri_CD)
palmieri_fit_CD <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "CD"],
                                              design = design_palmieri_CD),
                                        contrast_matrix_CD))
contrast_matrix_UC <- makeContrasts(I-nI, levels = design_palmieri_UC)
palmieri_fit_UC <- eBayes(contrasts.fit(lmFit(palmieri_final[,disease == "UC"],
                                              design = design_palmieri_UC),
                                        contrast_matrix_UC))
# chrohn disease
table_CD <- topTable(palmieri_fit_CD, number = Inf)
head(table_CD)
table(table_CD$adj.P.Val < 0.05)
table(table_CD$P.Value < 0.001)

hist(table_CD$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "inflammed vs non-imflamed - Crohn’s disease", xlab = "p-values")

# ulcerative collitus


table_UC <- topTable(palmieri_fit_UC, number = Inf)
head(table_UC)

table(table_UC$adj.P.Val < 0.05)
table(table_UC$P.Value < 0.001)

hist(table_UC$P.Value, col = brewer.pal(3, name = "Set2")[2],
     main = "inflammed vs non-imflamed - Ulcerative colitis", xlab = "p-values")



# compare the work to the already published results:

palmieri_DE_res <- sapply(1:4, function(i) read.xlsx(cols = 1,
                                                     "dat.xlsx",
                                                     sheet = i, startRow = 4))

names(palmieri_DE_res) <- c("CD_UP", "CD_DOWN", "UC_UP", "UC_DOWN")
palmieri_DE_res <- lapply(palmieri_DE_res, as.character)
paper_DE_genes_CD <- Reduce("c", palmieri_DE_res[1:2])
paper_DE_genes_UC <- Reduce("c", palmieri_DE_res[3:4])


overlap_CD <- length(intersect(subset(table_CD, P.Value < 0.001)$SYMBOL,
                               paper_DE_genes_CD)) / length(paper_DE_genes_CD)
overlap_UC <- length(intersect(subset(table_UC, P.Value < 0.001)$SYMBOL,
                               paper_DE_genes_UC)) / length(paper_DE_genes_UC)



DE_genes_CD <- subset(table_CD, adj.P.Val < 0.1)$PROBEID
back_genes_idx <- genefinder(palmieri_final, as.character(DE_genes_CD),
                             method="manhattan", scale="none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <-featureNames(palmieri_final)[back_genes_idx]
back_genes <- setdiff(back_genes, DE_genes_CD)
intersect(back_genes, DE_genes_CD)

multidensity(list(
  all= table_CD[,"AveExpr"] ,
  fore= table_CD[DE_genes_CD , "AveExpr"],
  back= table_CD[rownames(table_CD) %in% back_genes, "AveExpr"]),
  col = c("#e46981", "#ae7ee2", "#a7ad4a"),
  xlab="mean expression",
  main = "DE genes for CD - background - matching")




# now run the top go on DE genes

# in this section we are looking at enriched genes in the CD section

gene_IDs <- rownames(table_CD)
in_universe <- gene_IDs %in% c(DE_genes_CD , back_genes)
inSelection <- gene_IDs %in% DE_genes_CD
all_genes <- factor(as.integer(inSelection[in_universe]))
names(all_genes) <- gene_IDs[in_universe]


ont <- "BP"
top_GO_data <- new("topGOdata", ontology = ont, allGenes = all_genes,
                   nodeSize = 10, annot=annFUN.db, affyLib = "hugene10sttranscriptcluster.db")

result_top_GO_elim <- runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")


res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 100)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "hugene10sttranscriptcluster.db", geneCutOff = 1000)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), collapse = "")
})


session_info()



# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value                       
# version  R version 3.5.1 (2018-07-02)
# os       Debian GNU/Linux 9 (stretch)
# system   x86_64, linux-gnu           
# ui       RStudio                     
# language (EN)                        
# collate  en_US.UTF-8                 
# ctype    en_US.UTF-8                 
# tz       Etc/UTC                     
# date     2019-06-13                  
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                        * version   date       lib source        
# affxparser                       1.54.0    2018-10-30 [1] Bioconductor  
# affyio                           1.52.0    2018-10-30 [1] Bioconductor  
# annotate                       * 1.60.1    2019-03-07 [1] Bioconductor  
# AnnotationDbi                  * 1.44.0    2018-10-30 [1] Bioconductor  
# ArrayExpress                   * 1.42.0    2018-10-30 [1] Bioconductor  
# assertthat                       0.2.1     2019-03-21 [1] CRAN (R 3.5.1)
# backports                        1.1.4     2019-04-10 [1] CRAN (R 3.5.1)
# Biobase                        * 2.42.0    2018-10-30 [1] Bioconductor  
# BiocGenerics                   * 0.28.0    2018-10-30 [1] Bioconductor  
# BiocManager                      1.30.4    2018-11-13 [1] CRAN (R 3.5.1)
# BiocParallel                     1.16.6    2019-02-10 [1] Bioconductor  
# BiocStyle                      * 2.10.0    2018-10-30 [1] Bioconductor  
# biomaRt                        * 2.38.0    2018-10-30 [1] Bioconductor  
# Biostrings                     * 2.50.2    2019-01-03 [1] Bioconductor  
# bit                              1.1-14    2018-05-29 [1] CRAN (R 3.5.1)
# bit64                            0.9-7     2017-05-08 [1] CRAN (R 3.5.1)
# bitops                           1.0-6     2013-08-17 [1] CRAN (R 3.5.1)
# blob                             1.1.1     2018-03-25 [1] CRAN (R 3.5.1)
# callr                            3.2.0     2019-03-15 [1] CRAN (R 3.5.1)
# caTools                          1.17.1.2  2019-03-06 [1] CRAN (R 3.5.1)
# checkmate                        1.9.3     2019-05-03 [1] CRAN (R 3.5.1)
# cli                              1.1.0     2019-03-19 [1] CRAN (R 3.5.1)
# clusterProfiler                * 3.10.1    2018-12-20 [1] Bioconductor  
# codetools                        0.2-15    2016-10-05 [2] CRAN (R 3.5.1)
# colorspace                       1.4-1     2019-03-18 [1] CRAN (R 3.5.1)
# cowplot                          0.9.4     2019-01-08 [1] CRAN (R 3.5.1)
# crayon                           1.3.4     2017-09-16 [1] CRAN (R 3.5.1)
# DAAG                           * 1.22.1    2019-03-02 [1] CRAN (R 3.5.1)
# data.table                       1.12.2    2019-04-07 [1] CRAN (R 3.5.1)
# DBI                            * 1.0.0     2018-05-02 [1] CRAN (R 3.5.1)
# DelayedArray                     0.8.0     2018-10-30 [1] Bioconductor  
# desc                             1.2.0     2018-05-01 [1] CRAN (R 3.5.1)
# devtools                       * 2.0.2     2019-04-08 [1] CRAN (R 3.5.1)
# digest                           0.6.19    2019-05-20 [1] CRAN (R 3.5.1)
# DO.db                            2.9       2018-11-20 [1] Bioconductor  
# DOSE                             3.8.2     2019-01-14 [1] Bioconductor  
# dplyr                          * 0.8.1     2019-05-14 [1] CRAN (R 3.5.1)
# enrichplot                       1.2.0     2018-10-30 [1] Bioconductor  
# europepmc                        0.3       2018-04-20 [1] CRAN (R 3.5.1)
# evaluate                         0.14      2019-05-28 [1] CRAN (R 3.5.1)
# fansi                            0.4.0     2018-10-05 [1] CRAN (R 3.5.1)
# farver                           1.1.0     2018-11-20 [1] CRAN (R 3.5.1)
# fastmatch                        1.1-0     2017-01-28 [1] CRAN (R 3.5.1)
# ff                               2.2-14    2018-05-15 [1] CRAN (R 3.5.1)
# fgsea                            1.8.0     2018-10-30 [1] Bioconductor  
# foreach                          1.4.4     2017-12-12 [1] CRAN (R 3.5.1)
# fs                               1.3.1     2019-05-06 [1] CRAN (R 3.5.1)
# gdata                            2.18.0    2017-06-06 [1] CRAN (R 3.5.1)
# genefilter                     * 1.64.0    2018-10-30 [1] Bioconductor  
# geneplotter                    * 1.60.0    2018-10-30 [1] Bioconductor  
# GenomeInfoDb                     1.18.2    2019-02-12 [1] Bioconductor  
# GenomeInfoDbData                 1.2.0     2018-11-20 [1] Bioconductor  
# GenomicRanges                    1.34.0    2018-10-30 [1] Bioconductor  
# ggforce                          0.2.2     2019-04-23 [1] CRAN (R 3.5.1)
# ggplot2                        * 3.1.1     2019-04-07 [1] CRAN (R 3.5.1)
# ggplotify                        0.0.3     2018-08-03 [1] CRAN (R 3.5.1)
# ggraph                           1.0.2     2018-07-07 [1] CRAN (R 3.5.1)
# ggrepel                          0.8.1     2019-05-07 [1] CRAN (R 3.5.1)
# ggridges                         0.5.1     2018-09-27 [1] CRAN (R 3.5.1)
# glue                             1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
# GO.db                          * 3.7.0     2019-03-11 [1] Bioconductor  
# GOSemSim                         2.8.0     2018-10-30 [1] Bioconductor  
# gplots                         * 3.0.1.1   2019-01-27 [1] CRAN (R 3.5.1)
# graph                          * 1.60.0    2018-10-30 [1] Bioconductor  
# graphite                         1.28.2    2019-01-18 [1] Bioconductor  
# gridExtra                        2.3       2017-09-09 [1] CRAN (R 3.5.1)
# gridGraphics                     0.4-1     2019-05-20 [1] CRAN (R 3.5.1)
# gtable                           0.3.0     2019-03-25 [1] CRAN (R 3.5.1)
# gtools                           3.8.1     2018-06-26 [1] CRAN (R 3.5.1)
# highr                            0.8       2019-03-20 [1] CRAN (R 3.5.1)
# hms                              0.4.2     2018-03-10 [1] CRAN (R 3.5.1)
# htmltools                        0.3.6     2017-04-28 [1] CRAN (R 3.5.1)
# httr                             1.4.0     2018-12-11 [1] CRAN (R 3.5.1)
# hugene10sttranscriptcluster.db * 8.7.0     2019-06-13 [1] Bioconductor  
# igraph                           1.2.4.1   2019-04-22 [1] CRAN (R 3.5.1)
# IRanges                        * 2.16.0    2018-10-30 [1] Bioconductor  
# iterators                        1.0.10    2018-07-13 [1] CRAN (R 3.5.1)
# jsonlite                         1.6       2018-12-07 [1] CRAN (R 3.5.1)
# KernSmooth                       2.23-15   2015-06-29 [2] CRAN (R 3.5.1)
# knitr                          * 1.23      2019-05-18 [1] CRAN (R 3.5.1)
# labeling                         0.3       2014-08-23 [1] CRAN (R 3.5.1)
# lattice                        * 0.20-38   2018-11-04 [2] CRAN (R 3.5.1)
# latticeExtra                     0.6-28    2016-02-09 [1] CRAN (R 3.5.1)
# lazyeval                         0.2.2     2019-03-15 [1] CRAN (R 3.5.1)
# limma                          * 3.38.3    2018-12-02 [1] Bioconductor  
# LSD                            * 4.0-0     2018-01-26 [1] CRAN (R 3.5.1)
# magrittr                         1.5       2014-11-22 [1] CRAN (R 3.5.1)
# MASS                           * 7.3-51.1  2018-11-01 [2] CRAN (R 3.5.1)
# Matrix                           1.2-15    2018-11-01 [2] CRAN (R 3.5.1)
# matrixStats                    * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
# memoise                          1.1.0     2017-04-21 [1] CRAN (R 3.5.1)
# multcomp                       * 1.4-10    2019-03-05 [1] CRAN (R 3.5.1)
# munsell                          0.5.0     2018-06-12 [1] CRAN (R 3.5.1)
# mvtnorm                        * 1.0-10    2019-03-05 [1] CRAN (R 3.5.1)
# oligo                          * 1.46.0    2018-10-30 [1] Bioconductor  
# oligoClasses                   * 1.44.0    2018-10-30 [1] Bioconductor  
# openxlsx                       * 4.1.0.1   2019-05-28 [1] CRAN (R 3.5.1)
# org.Hs.eg.db                   * 3.7.0     2019-01-16 [1] Bioconductor  
# pd.hugene.1.0.st.v1            * 3.14.1    2019-06-13 [1] Bioconductor  
# pheatmap                       * 1.0.12    2019-01-04 [1] CRAN (R 3.5.1)
# pillar                           1.4.1     2019-05-28 [1] CRAN (R 3.5.1)
# pkgbuild                         1.0.3     2019-03-20 [1] CRAN (R 3.5.1)
# pkgconfig                        2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
# pkgload                          1.0.2     2018-10-29 [1] CRAN (R 3.5.1)
# plyr                             1.8.4     2016-06-08 [1] CRAN (R 3.5.1)
# polyclip                         1.10-0    2019-03-14 [1] CRAN (R 3.5.1)
# preprocessCore                   1.44.0    2018-10-30 [1] Bioconductor  
# prettyunits                      1.0.2     2015-07-13 [1] CRAN (R 3.5.1)
# processx                         3.3.1     2019-05-08 [1] CRAN (R 3.5.1)
# progress                         1.2.2     2019-05-16 [1] CRAN (R 3.5.1)
# ps                               1.3.0     2018-12-21 [1] CRAN (R 3.5.1)
# purrr                            0.3.2     2019-03-15 [1] CRAN (R 3.5.1)
# qvalue                           2.14.1    2019-01-10 [1] Bioconductor  
# R6                               2.4.0     2019-02-14 [1] CRAN (R 3.5.1)
# rappdirs                         0.3.1     2016-03-28 [1] CRAN (R 3.5.1)
# RColorBrewer                   * 1.1-2     2014-12-07 [1] CRAN (R 3.5.1)
# Rcpp                             1.0.1     2019-03-17 [1] CRAN (R 3.5.1)
# RCurl                            1.95-4.12 2019-03-04 [1] CRAN (R 3.5.1)
# reactome.db                      1.66.0    2019-01-23 [1] Bioconductor  
# ReactomePA                     * 1.26.0    2018-10-30 [1] Bioconductor  
# remotes                          2.0.4     2019-04-10 [1] CRAN (R 3.5.1)
# reshape2                         1.4.3     2017-12-11 [1] CRAN (R 3.5.1)
# rlang                            0.3.4     2019-04-07 [1] CRAN (R 3.5.1)
# rmarkdown                        1.13      2019-05-22 [1] CRAN (R 3.5.1)
# rprojroot                        1.3-2     2018-01-03 [1] CRAN (R 3.5.1)
# RSQLite                        * 2.1.1     2018-05-06 [1] CRAN (R 3.5.1)
# rstudioapi                       0.10      2019-03-19 [1] CRAN (R 3.5.1)
# rvcheck                          0.1.3     2018-12-06 [1] CRAN (R 3.5.1)
# S4Vectors                      * 0.20.1    2018-11-09 [1] Bioconductor  
# sandwich                         2.5-1     2019-04-06 [1] CRAN (R 3.5.1)
# scales                           1.0.0     2018-08-09 [1] CRAN (R 3.5.1)
# sessioninfo                      1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
# SparseM                        * 1.77      2017-04-23 [1] CRAN (R 3.5.1)
# stringi                          1.4.3     2019-03-12 [1] CRAN (R 3.5.1)
# stringr                        * 1.4.0     2019-02-10 [1] CRAN (R 3.5.1)
# SummarizedExperiment             1.12.0    2018-10-30 [1] Bioconductor  
# survival                       * 2.44-1.1  2019-04-01 [1] CRAN (R 3.5.1)
# TH.data                        * 1.0-10    2019-01-21 [1] CRAN (R 3.5.1)
# tibble                           2.1.3     2019-06-06 [1] CRAN (R 3.5.1)
# tidyr                            0.8.3     2019-03-01 [1] CRAN (R 3.5.1)
# tidyselect                       0.2.5     2018-10-11 [1] CRAN (R 3.5.1)
# topGO                          * 2.34.0    2018-10-30 [1] Bioconductor  
# triebeard                        0.3.0     2016-08-04 [1] CRAN (R 3.5.1)
# tweenr                           1.0.1     2018-12-14 [1] CRAN (R 3.5.1)
# UpSetR                           1.4.0     2019-05-22 [1] CRAN (R 3.5.1)
# urltools                         1.7.3     2019-04-14 [1] CRAN (R 3.5.1)
# usethis                        * 1.5.0     2019-04-07 [1] CRAN (R 3.5.1)
# utf8                             1.1.4     2018-05-24 [1] CRAN (R 3.5.1)
# vctrs                            0.1.0     2018-11-29 [1] CRAN (R 3.5.1)
# viridis                          0.5.1     2018-03-29 [1] CRAN (R 3.5.1)
# viridisLite                      0.3.0     2018-02-01 [1] CRAN (R 3.5.1)
# withr                            2.1.2     2018-03-15 [1] CRAN (R 3.5.1)
# xfun                             0.7       2019-05-14 [1] CRAN (R 3.5.1)
# XML                            * 3.98-1.20 2019-06-06 [1] CRAN (R 3.5.1)
# xml2                             1.2.0     2018-01-24 [1] CRAN (R 3.5.1)
# xtable                           1.8-4     2019-04-21 [1] CRAN (R 3.5.1)
# XVector                        * 0.22.0    2018-10-30 [1] Bioconductor  
# yaml                             2.2.0     2018-07-25 [1] CRAN (R 3.5.1)
# zeallot                          0.1.0     2018-01-28 [1] CRAN (R 3.5.1)
# zip                              2.0.2     2019-05-13 [1] CRAN (R 3.5.1)
# zlibbioc                         1.28.0    2018-10-30 [1] Bioconductor  
# zoo                              1.8-6     2019-05-28 [1] CRAN (R 3.5.1)











