# In this analysis we will dig into the mafft package in R.
# The mafft package is designed to analyze the TCGABiolinks Somatic 
# variant data. We will look to separate the data into Healthy vs Cancer,
# the second level of course will be cancer type specific variation
# followed by TNM staging correlation analysis and finally the 
# progression analysis for PanCancer progression. 
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)
library(TCGAbiolinks)
library(dplyr)
library(data.table)

query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")


maf <- GDCquery_Maf("BRCA", pipelines = "muse") %>% read.maf(clinicalData = clin)

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = maf, top = 20, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)


#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1), fontSize = .4)


drugInteractions(maf = maf, fontSize = .75)
OncogenicPathways(maf = maf)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BSgenome")
library(BSgenome)
brca.tnm = trinucleotideMatrix(maf = maf, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
BSgenome::available.genomes()

