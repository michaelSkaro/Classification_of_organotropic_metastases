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

maf <- GDCquery_Maf("BRCA", pipelines = "muse", ) %>% read.maf

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

oncoplot(maf = maf, top = 20, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)


#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1), fontSize = .4)


drugInteractions(maf = maf, fontSize = .8)
OncogenicPathways(maf = maf)

devtools::install_github(repo = "PoisonAlien/TCGAmutations")

# complete the analysis for all of the projects :)
#cite https://doi.org/10.1016/j.cels.2018.03.00

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

proj <- projects[1]

for(proj in projects){
  qmaf <- substr(proj, 6,10)
  
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  
  ggplot2::ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/somatic_variation_analysis/{proj}"),"_somatic_variation_summary",".pdf"),
         plot = plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE), device = "pdf", width = 11.5, height =8.5, units = "in")
  
}


for(proj in projects){
  qmaf <- substr(proj, 6,10)
  
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1), fontSize = .4)
  
  
}

for(proj in projects){
  qmaf <- substr(proj, 6,10)
  
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  drugInteractions(maf = maf, fontSize = .8)
  
  
}

for(proj in projects){
  qmaf <- substr(proj, 6,10)
  
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  ggplot2::ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/somatic_variation_analysis/{proj}"),"_OncoGenicPathways",".pdf"),
                  plot =  OncogenicPathways(maf = maf), device = "pdf", width = 11, height =8.5, units = "in")
  
}

projects
proj <- projects[5]
qmaf <- substr(proj, 6,10)
maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
oncoplot(maf = maf, top = 20, removeNonMutated = TRUE)  

WhatsLeftChap1{
  #MAFsplit(TNM-CancerType):{
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #1:summaryPlot()
  #2:OncoPlot()
  #3:OncogeneticPathways()
  #4:druggableCataogires()
  #
  #}
  #
  #
  #PanCancerMAF(TNM):{
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #1:summaryPlot()
  #2:OncoPlot()
  #3:OncogeneticPathways()
  #4:druggableCataogires()
  #
  #}
  #
  #Methylation(CancerType-DMR-Up){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #}
  #
  #Methylation(CancerType-DMR-Down){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #}
  #
  #Methylation(PanCancer-DMR-Up){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #}
  #
  #Methylation(PanCancer-DMR-Down){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #}
  #
  #eQTL(CancerType){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #
  #
  #}
  #
  #eQTL(PanCancer){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #
  #
  #}
  #
  #PathwayAnalysisParadigm(CancerType){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #
  #}
  #
  #PathwayAnalysisParadigm(PanCancer){
  #	T1
  #	T2
  #	T3
  #	T4
  #		N1
  #		N2
  #		N3
  #			M1
  #			M0
  #
  #}
  #}
