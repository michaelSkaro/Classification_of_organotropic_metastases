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

for(proj in projects){
  qmaf <- substr(proj, 6,10)
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
  ggplot2::ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/somatic_variation_analysis/{proj}"),"_Transition_vs_Transversion",".pdf"),
                  plot =  plotTiTv(res = titv), device = "pdf", width = 11, height =8.5, units = "in")

}

projects
proj <- projects[13]
qmaf <- substr(proj, 6,10)
maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)

maf@clinical.data$T_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "T\\w")
maf@clinical.data$N_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "N\\w")
maf@clinical.data$M_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "M\\w")

maf@clinical.data <-  maf@clinical.data %>%
  mutate(maf@clinical.data, M_stage = ifelse(maf@clinical.data$M_stage == "M1", 
                                         "M1", "M0"))

maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$T_stage != "TX")

maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$T_stage != "Ti")

maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$N_stage != "NX")

maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$M_stage != "MX")


#Color coding for patholic TNM stage classification
stage_cols = RColorBrewer::brewer.pal(n = 4,name = 'Spectral')
names(stage_cols) = c("T1","T2","T3","T4")
stage_cols = list(T_stage = stage_cols)



oncoplot(
  maf = maf,
  clinicalFeatures = c('T_stage','N_stage','M_stage'),
  sortByAnnotation = TRUE,
  annotationColor = stage_cols
)

proj

# remember to clear the plot!!!


# split the MAF objects by their clinical stages and conuct the same analyses.


projects
proj <- projects[1]
qmaf <- substr(proj, 6,10)

maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
maf_clin <- getClinicalData(maf)
maf_clin$T_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "T\\w")
T1.samples <- maf_clin[maf_clin$T_stage == 'T1',]
T1.samples <- T1.samples$T_stage

T1.maf <- maftools::subsetMaf(maf = maf, tsb = T1.samples, mafObj = TRUE)


#Color coding for patholic TNM stage classification
stage_cols = RColorBrewer::brewer.pal(n = 1,name = 'Spectral')
names(stage_cols) = c("T1")
stage_cols = list(T_stage = stage_cols)



oncoplot(
  maf = T1_maf,
  clinicalFeatures = c('T_stage'),
  sortByAnnotation = TRUE
)

proj

#WhatsLeftChap1{
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

# we will get a list of the samples in each of samples in each of the stages. Conduct the MAF splits. This will produce the associated outputs. We will be
# looking for different genes enriched in the changes in the transition to transversion rates, changes in the genes in the oncoplot, changes in the pathways
# in the OncogeneticPathways or changes in the druggable targets in the druggableCatagories.


projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

proj <- projects[1]

for(proj in projects){
  qmaf <- substr(proj, 6,10)
  
  maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
  
  # establish the correct splits
  
  TNM <- data.table::fread(str_glue("~/storage/Machine_Learning/TNM_stage/TNM_T/{proj}_TNM.csv"), header = TRUE) %>%
    dplyr::select(-"V1")
  
  
  
  
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAutils")
BiocManager::install("curatedTCGAData")

library(TCGAutils)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(RTCGAToolbox)
library(BiocFileCache)
library(rtracklayer)
library(R.utils)


# Cancer specific TNM staging
proj
projects
proj <- projects[12]
qmaf <- substr(proj, 6,10)
maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
maf@clinical.data$T_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "T\\w")
maf@clinical.data$N_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "N\\w")
maf@clinical.data$M_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "M\\w")
maf@clinical.data <-  maf@clinical.data %>%
  mutate(maf@clinical.data, M_stage = ifelse(maf@clinical.data$M_stage == "M1",
                                             "M1", "M0"))
maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$T_stage != "TX")
maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$T_stage != "Ti")
maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$N_stage != "NX")
maf@clinical.data <- maf@clinical.data %>%
  dplyr::filter(maf@clinical.data$M_stage != "MX")
#Color coding for patholic TNM stage classification
stage_cols = RColorBrewer::brewer.pal(n = 4,name = 'Spectral')
names(stage_cols) = c("T1","T2","T3","T4")
stage_cols = list(T_stage = stage_cols)
oncoplot(
  maf = maf,
  clinicalFeatures = c('T_stage','N_stage','M_stage'),
  sortByAnnotation = TRUE,
  annotationColor = stage_cols
)

#Cancer Specific stage by stage T
proj
T1_clin
T1_clin = T1_clin[T_stage %in% 'T1', Tumor_Sample_Barcode]
T1_clin <- maf@clinical.data %>%
dplyr::filter(maf@clinical.data$T_stage == "T1")
projects
proj <- projects[1]
qmaf <- substr(proj, 6,10)
maf <- TCGAmutations::tcga_load(study = qmaf, source = "MC3", reassign = TRUE)
maf@clinical.data$T_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "T\\w")
T1_clin <- maf@clinical.data %>%
dplyr::filter(maf@clinical.data$T_stage == "T1")
T1_maf <- maftools::subsetMaf(maf = maf, tsb = T1_clin, mafObj = TRUE)
T1_clin <- maf %>%
dplyr::filter(maf@clinical.data$T_stage == "T1")
maf_clin <- getClinicalData(maf)
maf_clin$T_stage <- str_extract(maf@clinical.data$stage_event_tnm_categories, "T\\w")
View(maf_clin)
T1.samples <- maf_clin[T_stage %in% 'T1', Tumor_Sample_Barcode]
T1.maf <- maftools::subsetMaf(maf = maf, tsb = T1_clin, mafObj = TRUE)
table(maf_clin$T_stage)
T1.samples <- maf_clin[maf_clin$T_stage == 'T1',]
T1.samples <- T1.samples$T_stage
T1.maf <- maftools::subsetMaf(maf = maf, tsb = T1.samples, mafObj = TRUE)
projects

#Cancer Specific stage by stage N
#Cancer Specific stage by stage M


# PanCancer Aggregation of samples
# PanCancer Aggregation of samples T
# PanCancer Aggregation of samples N
# PanCancer Aggregation of samples M

# Pancancer Progression by location

#1:summaryPlot()
#2:OncoPlot()
#3:OncogeneticPathways()
#4:druggableCataogires()


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
