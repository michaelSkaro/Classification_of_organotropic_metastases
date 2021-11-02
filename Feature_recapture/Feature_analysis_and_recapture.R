library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
setwd("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets")
files_in_dir <- list.files()
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$sample_type == "Primary Tumor",]
tumor.samples <- tumor.samples[tumor.samples$project %in% projects,]

tumor.samples$bcr_patient_barcode <- substr(tumor.samples$barcode_short, 0, 12)
clinical$barcode_short <- substr(clinical$barcode, 0,16)

for(proj in projects){
  
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
      # complete data set of IG selected features. The train test splits will be merged 
      val <- str_sub(proj, 6,9)
      dat <- rbind(train, test)
      organ_labels <- dat %>%dplyr::select(org)
      
      dat <- dat %>%
        dplyr::select(-org)
      dat <- as.data.frame(t(log2(dat)))
      names(dat) <- organ_labels[,1]
      
      dat[mapply(is.infinite, dat)] <- 0
      basemean = rowMeans(dat)
      type = gsub("s\\d+_", "", colnames(dat))
      dat <- as.matrix(dat)
      ha = HeatmapAnnotation(type = type, annotation_name_side ="right")
      pdf(str_glue("heatmap_{proj}_{org}_log2_rasteration.pdf"), width = 5, height = 5)
      ht_list = Heatmap(dat, name = " Log2 counts",
                        border = TRUE,
                        cluster_columns =TRUE,
                        cluster_rows = TRUE,
                        row_dend_reorder = TRUE,
                        column_dend_reorder = TRUE,
                        show_heatmap_legend = TRUE,
                        show_column_names = FALSE, show_row_names = FALSE,
                        show_column_dend = TRUE, show_row_dend = FALSE, 
                        use_raster = TRUE, raster_resize = TRUE,
                        top_annotation = ha,
                        column_title = paste0(val," ", "Metastasize to ",str_glue("{org}")),
                        row_title = "Selected features")
      #top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), height = unit(2, "mm"))),
      #width = unit(15, "mm"))
      
      draw(ht_list)
      
      
      dev.off()
      
    }
  }
}

library(UpSetR)
library(stringr)
library(GeneOverlap)


for(proj in projects){
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
      # complete data set of IG selected features. The train test splits will be merged 
      val <- str_sub(proj, 6,9)
      dat <- rbind(train, test)
      organ_labels <- dat %>%dplyr::select(org)
      dat <- dat %>%
        dplyr::select(-org)
      
      # cancer type, metastatic location organ specific gene set enrichment 
      
      ego <-clusterProfiler::enrichGO(gene  = substring(colnames(dat),1,15),
                                      OrgDb         = "org.Hs.eg.db",
                                      keyType       = 'ENSEMBL',
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
      
      #Get only  the significant results from enrichment
      res <- ego@result
      res <- res[res$pvalue<0.05,]
      write.csv2(res, file = str_glue("ego_{proj}_{org}.csv"))
      
    }
  }
}

# semantic analysis of pathways and the overlaps
library(tibble)
library(dplyr)
library(plyr)
compare_list <- as.data.frame(data.table::fread("compare_list.csv", header = TRUE))

i <- 1
out <- as.data.frame(t(c("CancerType1"= NA,
                         "CancerType2"= NA,
                         "Seeding_Location" = NA,
                         "odds.ratio" = NA,
                         "Enriched Processes CT1" = NA,
                         "Enriched Processes CT2" = NA,
                         "intersection"= NA,
                         "p.value" = NA)))
for(i in 1:length(CT1)){
  print(i)
  print(paste0("ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  print(paste0("ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  c1 <- data.table::fread(paste0("ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  c2 <- data.table::fread(paste0("ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  
  # add all of the metric to a dataframe. 
  
  
    go.obj <- newGeneOverlap(listA = c1$ID, listB = c2$ID, genome.size = 23393)
    go.obj <- testGeneOverlap(go.obj)
    df <- as.data.frame(t(c("CancerType1"= compare_list[i,1],
               "CancerType2"= compare_list[i,2],
               "Seeding_Location" = compare_list[i,3],
               "Enriched Processes CT1" = length(c1$ID),
               "Enriched Processes CT2" = length(c2$ID),
               "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
               "intersection"= length(go.obj@intersection),
               "p.value" = go.obj@pval)))
    out <- rbind(out,df)
    
    
}
out <- out[2:length(out$CancerType1),]
out <- out[sort(out$p.value),]  
write.csv2(out, file = "/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/Fisher_exact_test_semantic.csv")
  
sessionInfo()

library(simplifyEnrichment)
library(stringr)
library(tibble)
dat <- data.table::fread("ego_TCGA-HNSC_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")

mat = GO_similarity(dat$ID, ont = "BP")
simplifyGO(mat, method = "binary_cut",column_title = "Biological Processes clustered by Sematic Similarity", plot = TRUE)

dat <- data.table::fread("ego_TCGA-BLCA_Lymph_Node.csv",header = TRUE) %>%
  column_to_rownames("V1")
mat = GO_similarity(dat$ID, ont = "BP")
simplifyGO(mat, method = "binary_cut",column_title = "Biological Processes clustered by Sematic Similarity", plot = TRUE)

# make upsetR plot to finish it off

library(UpSetR)
seed_loc <- unique(c("Lymph_Node","Lung","Prostate","Lung","Lung","Liver","Lung","Liver","Lung","Liver","Liver","Lung","Lung","Lung","Liver","Lung","Liver","Bone","Lung"))

#Bone
BLCA_Bone <- data.table::fread("ego_TCGA-BLCA_Bone.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA_Bone <- data.table::fread("ego_TCGA-BRCA_Bone.csv", header = TRUE) %>%
  column_to_rownames("V1")
listInput <- list('BLCA Bone' = BLCA_Bone$ID ,'BRCA Bone' = BRCA_Bone$ID)




# Liver 
BLCA_Liver <- data.table::fread("ego_TCGA-BLCA_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA_Liver <- data.table::fread("ego_TCGA-BRCA_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
LIHC_Liver <- data.table::fread("ego_TCGA-LIHC_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")
COAD_Liver <- data.table::fread("ego_TCGA-COAD_Liver.csv", header = TRUE) %>%
  column_to_rownames("V1")


library(UpSetR)
listInput <- list('BLCA Liver' = BLCA_Liver$ID ,'BRCA Liver' = BRCA_Liver$ID ,
                  'LIHC Liver' = LIHC_Liver$ID,'COAD Liver' = COAD_Liver$ID)

UpSetR::upset(fromList(listInput), order.by ="freq")





# Lung
BLCA<- data.table::fread("ego_TCGA-BLCA_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
BRCA<- data.table::fread("ego_TCGA-BRCA_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
LIHC<- data.table::fread("ego_TCGA-LIHC_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
HNSC<- data.table::fread("ego_TCGA-HNSC_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")
LUAD<- data.table::fread("ego_TCGA-LUAD_Lung.csv", header = TRUE) %>%
  column_to_rownames("V1")

library(UpSetR)
listInput <- list('BLCA Lung' = BLCA$ID ,'BRCA Lung' = BRCA$ID , "HNSC Lung"=HNSC$ID,
                  'LIHC Lung' = LIHC$ID,'LUAD Lung' = LUAD$ID)

# Lymph Node 

#Lymph_Node
BLCA_Lymph_Node <- data.table::fread("ego_TCGA-BLCA_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")
HNSC_Lymph_Node <- data.table::fread("ego_TCGA-HNSC_Lymph_Node.csv", header = TRUE) %>%
  column_to_rownames("V1")
listInput <- list('BLCA Lymph_Node' = BLCA_Lymph_Node$ID ,'HNSC Lymph_Node' = HNSC_Lymph_Node$ID)

d <- UpSetR::upset(fromList(listInput), order.by ="freq")

dat <- intersect(HNSC_Lymph_Node$ID, BLCA_Lymph_Node$ID)
mat = GO_similarity(dat, ont = "BP")
simplifyGO(mat, method = "binary_cut",column_title = "Biological Processes clustered by Sematic Similarity", plot = TRUE)


# picture time :)

library(tidyverse)
library(tibble)
library(dplyr)
library(tidyr)


setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels")

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]

for(proj in projects){
  print(proj)
  for(org in organs){
    print(org)
    # volcano
    if(file.exists(str_glue("{proj}_met_{org}_RNAseq.csv"))){
      
      # read in DE file, make the filters
      res <- data.table::fread(str_glue("{proj}_met_{org}_RNAseq.csv")) %>%
        tibble::column_to_rownames("V1") %>%
        #dplyr::filter(padj <0.05) %>%
        #dplyr::filter(abs(log2FoldChange) > 0.5) %>%
        tibble::rownames_to_column("ENSEMBL")
      # remove NAs if they exist
      res <- res[complete.cases(res),]
      
      # add a column of NAs
      res$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "DOWN"
      
      # make the 
      res$delabel <- NA
      res$delabel[res$diffexpressed != "NO"] <- res$ENSEMBL[res$diffexpressed != "NO"]
      
      ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
        geom_point() + 
        theme_minimal() +
        geom_text()
      
      
      # Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
      # load library
      library(ggrepel)
      # plot adding up all layers we have seen so far
      p1 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "gray", "red")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      
      p2 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("purple", "black", "orange")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/Volcano/{proj}_{org}"),"_DE_vol_",".pdf"),
            plot = p1, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/Volcano/{proj}_{org}_2"),"_DE_vol_",".pdf"),
             plot = p2, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      
      
      
    }
    
  }
}




# Bar plot where up regs are red and down regs are red: Lung, Bone,Liver,LymphNode
# Lung upset plot
BLCA <- data.table::fread(str_glue("TCGA-BLCA_met_Lung_RNAseq.csv")) %>%
  tibble::column_to_rownames("V1") %>%
  dplyr::filter(padj <0.05) %>%
  dplyr::filter(abs(log2FoldChange) > 0.5) %>%
  tibble::rownames_to_column("ENSEMBL")
# remove NAs if they exist
BLCA <- BLCA[complete.cases(BLCA),]

# add a column of NAs
BLCA$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
BLCA$diffexpressed[BLCA$log2FoldChange > 0.5 & BLCA$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
BLCA$diffexpressed[BLCA$log2FoldChange < -0.5 & BLCA$padj < 0.05] <- "DOWN"

# make the 
BLCA$delabel <- NA
BLCA$delabel[BLCA$diffexpressed != "NO"] <- BLCA$ENSEMBL[BLCA$diffexpressed != "NO"]


BRCA <- data.table::fread(str_glue("TCGA-BRCA_met_Lung_RNAseq.csv")) %>%
  tibble::column_to_rownames("V1") %>%
  dplyr::filter(padj <0.05) %>%
  dplyr::filter(abs(log2FoldChange) > 0.5) %>%
  tibble::rownames_to_column("ENSEMBL")
# remove NAs if they exist
BRCA <- BRCA[complete.cases(BRCA),]

# add a column of NAs
BRCA$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
BRCA$diffexpressed[BRCA$log2FoldChange > 0.5 & BRCA$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
BRCA$diffexpressed[BRCA$log2FoldChange < -0.5 & BRCA$padj < 0.05] <- "DOWN"

# make the 
BRCA$delabel <- NA
BRCA$delabel[BRCA$diffexpressed != "NO"] <- BRCA$ENSEMBL[BRCA$diffexpressed != "NO"]



HNSC <- data.table::fread(str_glue("TCGA-HNSC_met_Lung_RNAseq.csv")) %>%
  tibble::column_to_rownames("V1") %>%
  dplyr::filter(padj <0.05) %>%
  dplyr::filter(abs(log2FoldChange) > 0.5) %>%
  tibble::rownames_to_column("ENSEMBL")
# remove NAs if they exist
HNSC <- HNSC[complete.cases(HNSC),]

# add a column of NAs
HNSC$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
HNSC$diffexpressed[HNSC$log2FoldChange > 0.5 & HNSC$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
HNSC$diffexpressed[HNSC$log2FoldChange < -0.5 & HNSC$padj < 0.05] <- "DOWN"

# make the 
HNSC$delabel <- NA
HNSC$delabel[HNSC$diffexpressed != "NO"] <- HNSC$ENSEMBL[HNSC$diffexpressed != "NO"]

LIHC <- data.table::fread(str_glue("TCGA-LIHC_met_Lung_RNAseq.csv")) %>%
  tibble::column_to_rownames("V1") %>%
  dplyr::filter(padj <0.05) %>%
  dplyr::filter(abs(log2FoldChange) > 0.5) %>%
  tibble::rownames_to_column("ENSEMBL")
# remove NAs if they exist
LIHC <- LIHC[complete.cases(LIHC),]

# add a column of NAs
LIHC$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
LIHC$diffexpressed[LIHC$log2FoldChange > 0.5 & LIHC$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
LIHC$diffexpressed[LIHC$log2FoldChange < -0.5 & LIHC$padj < 0.05] <- "DOWN"

# make the 
LIHC$delabel <- NA
LIHC$delabel[LIHC$diffexpressed != "NO"] <- LIHC$ENSEMBL[LIHC$diffexpressed != "NO"]


LUAD <- data.table::fread(str_glue("TCGA-LUAD_met_Lung_RNAseq.csv")) %>%
  tibble::column_to_rownames("V1") %>%
  dplyr::filter(padj <0.05) %>%
  dplyr::filter(abs(log2FoldChange) > 0.5) %>%
  tibble::rownames_to_column("ENSEMBL")
# remove NAs if they exist
LUAD <- LUAD[complete.cases(LUAD),]

# add a column of NAs
LUAD$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
LUAD$diffexpressed[LUAD$log2FoldChange > 0.5 & LUAD$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
LUAD$diffexpressed[LUAD$log2FoldChange < -0.5 & LUAD$padj < 0.05] <- "DOWN"

# make the 
LUAD$delabel <- NA
LUAD$delabel[LUAD$diffexpressed != "NO"] <- LUAD$ENSEMBL[LUAD$diffexpressed != "NO"]

library(UpSetR)
listInput <- list('BLCA Lung' = BLCA$ENSEMBL ,'BRCA Lung' = BRCA$ENSEMBL , "HNSC Lung"=HNSC$ENSEMBL,
                  'LIHC Lung' = LIHC$ENSEMBL,'LUAD Lung' = LUAD$ENSEMBL)
UpSetR::upset(fromList(listInput), order.by ="freq")



###############################################
###############################################
###############################################

setwd("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets")

files_in_dir <- list.files()
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[1]

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$sample_type == "Primary Tumor",]
tumor.samples <- tumor.samples[tumor.samples$project %in% projects,]

tumor.samples$bcr_patient_barcode <- substr(tumor.samples$barcode_short, 0, 12)
clinical$barcode_short <- substr(clinical$barcode, 0,16)



# Up and down regulated pathways in positive tumors
for(proj in projects){
  for(org in organs){
    if(file.exists(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv"))){
      train <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_train.csv")))
      test <- as.data.frame(data.table::fread(str_glue("{proj}_metastatic_data_RNAseq_{org}_feature_selected_test.csv")))
      # complete data set of IG selected features. The train test splits will be merged 
      val <- str_sub(proj, 6,9)
      dat <- rbind(train, test)
      organ_labels <- dat %>%dplyr::select(org)
      
      dat <- dat %>%
        dplyr::select(-org)
      dat <- as.data.frame(t(dat))
      names(dat) <- organ_labels[,1]
      
      dat[mapply(is.infinite, dat)] <- 0
      dat <- ceiling(dat)
      
      drop <- c("Negative")
      dat = dat[,!(names(dat) %in% drop)]
      
      # replace the period in the colnames for an underscore
      
      newNames <- str_replace_all(colnames(dat), pattern = "\\.", replacement = "_")
      colnames(dat) <- newNames
      colnames(dat)
      # read in the normal data
      
      df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
        as_tibble() %>%
        tibble::column_to_rownames(var = "Ensembl")
      
      n<-dim(df.exp)[1]
      df.exp<-df.exp[1:(n-5),]
      
      # subset to only normal samples
     
      
      coldata.n <- normal.samples[normal.samples$project == proj,]
      
      df.exp <- df.exp[ ,colnames(df.exp) %in% coldata.n$barcode]
      rownames(coldata.n) <- coldata.n$barcode
      coldata.n$sample_type <- gsub(" ", "_", x = coldata.n$sample_type)
      
      # make a coldata for the current dat file
      
      coldata.t <- as.data.frame(colnames(dat))
      coldata.t$sample_type <- "Primary_Tumor"
      colnames(coldata.t)[1] <- "Samples"
      # mthe the rownames match the samples column names 
      rownames(coldata.t) <- coldata.t$Samples
      
      # extract only the rows we need for the selected features
      
      df.exp <- df.exp[rownames(df.exp) %in% rownames(dat),]
      
      # make a column value to left join the two data sets
      
      df.exp <- df.exp %>%
        tibble::rownames_to_column("ENSEMBL")
      dat <- dat %>%
        tibble::rownames_to_column("ENSEMBL")
      
      df.exp <- left_join(df.exp,dat, by ="ENSEMBL")
      
      df.exp <- df.exp %>%
        column_to_rownames("ENSEMBL")
      
      
      # merge the coldata files
      # I think there is something with ordering the columns as well
      
      coldata.n <- coldata.n %>%
        dplyr::select(barcode,sample_type)
      colnames(coldata.n)[1] <- "Samples"
      
      
      coldata <- rbind(coldata.n,coldata.t)
      
      coldata <- coldata[order(coldata$Samples)]
      df.exp <- df.exp[order(colnames(df.exp))]
      
      rownames(coldata) %in% colnames(df.exp)
      
      # Differential expression for just selected features
      # DEseq2
      dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ sample_type)
      
      dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
      
      dds <- DESeq(dds)
      
      save(dds, file = str_glue("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/Rdata/{proj}_DE_met_{org}.RData"))
      
      res <- results(dds)
      resOrdered <- res[order(res$pvalue),]
      resOrdered <- as.data.frame(resOrdered)
      res <- as.data.frame(res)
      
      write.csv(resOrdered, file = str_glue("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/res/{proj}_DE_met_{org}.csv"))
      
      # volcano, enrichment for the DE genes
      
      res <- res %>%
        tibble::rownames_to_column("ENSEMBL")
      
      
      res <- res[complete.cases(res),]
      
      # add a column of NAs
      res$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      res$diffexpressed[res$log2FoldChange > 0.5 & res$padj < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      res$diffexpressed[res$log2FoldChange < -0.5 & res$padj < 0.05] <- "DOWN"
      
      # make the 
      res$delabel <- NA
      res$delabel[res$diffexpressed != "NO"] <- res$ENSEMBL[res$diffexpressed != "NO"]
      
      ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
        geom_point() + 
        theme_minimal() +
        geom_text()
      
      
      # Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
      # load library
      library(ggrepel)
      # plot adding up all layers we have seen so far
      p1 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "gray", "red")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      
      p2 <- ggplot(data=res, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label="")) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("purple", "black", "orange")) +
        geom_vline(xintercept=c(-0.5, 0.5), col="red") +
        geom_hline(yintercept=-log10(0.001), col="red") +
        labs(title = str_glue("{proj} tumors metastasize to {org}"))
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/Volcano/{proj}_{org}"),"_DE_vol_",".pdf"),
             plot = p1, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/Volcano/{proj}_{org}_2"),"_DE_vol_",".pdf"),
             plot = p2, device = "pdf", width = 5, height = 3, units = "in", dpi = "retina")
      
      
      
    }
  }
}

###############################################
###############################################
###############################################
setwd("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/res")
library(tidyverse)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
library(cowplot)
# read in the res values 
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD", "TCGA-HNSC", "TCGA-LUAD", "TCGA-LIHC")
proj <- projects[1]
organs <- c("Bladder", "Liver", "Lung", "Bone", "Lymph_Node", "Pelvis", "Prostate")
org <- organs[3]

#Lung
# gene set enrichment analysis for the DE genes
# hallmark human cancer gene set
all_gene_sets = msigdbr(species = "Homo sapiens")
h_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

for(proj in projects){
  print(str_glue("Processing: {proj}"))
  for(org in organs){
    print(str_glue("Analyzing: {org}"))
    if(file.exists(str_glue("{proj}_DE_met_{org}.csv"))){
      dat <- data.table::fread(str_glue("{proj}_DE_met_{org}.csv", header = TRUE)) %>%
        tibble::column_to_rownames("V1") %>%
        #dplyr::filter(padj <0.05) %>%
        #dplyr::filter(abs(log2FoldChange) > 0.5) %>%
        tibble::rownames_to_column("ENSEMBL")
      # remove NAs if they exist
      dat <- dat[complete.cases(dat),]
      dat$ENSEMBL <- substr(dat$ENSEMBL, 1, 15)
      
      geneList = dat[,1]
      
      ## feature 2: named vector
      names(geneList) = dat[,3]
      
      ## feature 3: decreasing orde
      geneList = sort(geneList, decreasing = TRUE)
      
      
      ego <-clusterProfiler::enrichGO(gene  = geneList,
                                      OrgDb         = "org.Hs.eg.db",
                                      keyType       = 'ENSEMBL',
                                      ont           = "MF",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
      
      #Get only  the significant results from enrichment
      res <- ego@result
      res <- res[res$pvalue<0.05,]
      write.csv(res, file = str_glue("MF_ego_{proj}_{org}.csv"))
      
      
    }
  }
}



# Modified pathways and the overlaps
library(tibble)
library(dplyr)
library(plyr)
compare_list <- as.data.frame(data.table::fread("/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/compare_list.csv", header = TRUE))

i <- 1
out <- as.data.frame(t(c("CancerType1"= NA,
                         "CancerType2"= NA,
                         "Seeding_Location" = NA,
                         "odds.ratio" = NA,
                         "Enriched Processes CT1" = NA,
                         "Enriched Processes CT2" = NA,
                         "intersection"= NA,
                         "p.value" = NA)))

for(i in 1:length(compare_list$CancerType1)){
  print(i)
  print(paste0("BP_ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  print(paste0("BP_ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  c1 <- data.table::fread(paste0("BP_ego","_",compare_list[i,1],"_",compare_list[i,3],".csv"))
  c2 <- data.table::fread(paste0("BP_ego","_",compare_list[i,2],"_",compare_list[i,3],".csv"))
  
  # add all of the metric to a dataframe. 
  
  
  go.obj <- newGeneOverlap(listA = c1$ID, listB = c2$ID, genome.size = 23393)
  go.obj <- testGeneOverlap(go.obj)
  df <- as.data.frame(t(c("CancerType1"= compare_list[i,1],
                          "CancerType2"= compare_list[i,2],
                          "Seeding_Location" = compare_list[i,3],
                          "Enriched Processes CT1" = length(c1$ID),
                          "Enriched Processes CT2" = length(c2$ID),
                          "odds.ratio" = sprintf("%.3f", go.obj@odds.ratio),
                          "intersection"= length(go.obj@intersection),
                          "p.value" = go.obj@pval)))
  out <- rbind(out,df)
  
  
}
out <- out[2:length(out$CancerType1),]
out <- out[order(out$p.value),]  
write.csv(out, file = "/mnt/storage/mskaro1/Machine_Learning/All_MOT_selected_features/feature-selected-datasets/res/Fisher_exact_test_MF_semantic.csv")





















# R version 4.0.2 (2020-06-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04 LTS
# 
# Matrix products: default
# BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-openmp/libopenblasp-r0.3.8.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats4    grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] UpSetR_1.4.0              ComplexHeatmap_2.7.1      plyr_1.8.6                stringr_1.4.0            
# [5] GeneOverlap_1.24.0        org.Hs.eg.db_3.11.4       AnnotationDbi_1.50.3      IRanges_2.22.2           
# [9] S4Vectors_0.26.1          Biobase_2.48.0            simplifyEnrichment_0.99.5 BiocGenerics_0.34.0      
# [13] tibble_3.0.4              dplyr_1.0.2              
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.14.0           colorspace_2.0-0       rjson_0.2.20           ellipsis_0.3.1         ggridges_0.5.2        
# [6] circlize_0.4.11        qvalue_2.20.0          GlobalOptions_0.1.2    clue_0.3-57            rstudioapi_0.13       
# [11] farver_2.0.3           urltools_1.7.3         graphlayouts_0.7.1     ggrepel_0.8.2          bit64_4.0.5           
# [16] fansi_0.4.1            scatterpie_0.1.5       xml2_1.3.2             splines_4.0.2          GOSemSim_2.14.2       
# [21] knitr_1.30             polyclip_1.10-0        jsonlite_1.7.1         Cairo_1.5-12.2         cluster_2.1.0         
# [26] GO.db_3.11.4           png_0.1-7              ggforce_0.3.2          BiocManager_1.30.10    compiler_4.0.2        
# [31] httr_1.4.2             rvcheck_0.1.8          assertthat_0.2.1       Matrix_1.2-18          cli_2.1.0             
# [36] tweenr_1.0.1           htmltools_0.5.0        prettyunits_1.1.1      tools_4.0.2            igraph_1.2.6          
# [41] NLP_0.2-1              gtable_0.3.0           glue_1.4.2             reshape2_1.4.4         DO.db_2.9             
# [46] tinytex_0.27           fastmatch_1.1-0        Rcpp_1.0.5             enrichplot_1.8.1       slam_0.1-47           
# [51] vctrs_0.3.4            ggraph_2.0.3           xfun_0.19              lifecycle_0.2.0        clusterProfiler_3.16.1
# [56] gtools_3.8.2           DOSE_3.14.0            europepmc_0.4          MASS_7.3-53            scales_1.1.1          
# [61] tidygraph_1.2.0        hms_0.5.3              RColorBrewer_1.1-2     yaml_2.2.1             memoise_1.1.0         
# [66] gridExtra_2.3          ggplot2_3.3.2          downloader_0.4         triebeard_0.3.0        stringi_1.5.3         
# [71] RSQLite_2.2.1          caTools_1.18.0         BiocParallel_1.22.0    shape_1.4.5            bitops_1.0-6          
# [76] rlang_0.4.8            pkgconfig_2.0.3        matrixStats_0.57.0     evaluate_0.14          lattice_0.20-41       
# [81] purrr_0.3.4            cowplot_1.1.0          bit_4.0.4              tidyselect_1.1.0       magrittr_1.5          
# [86] R6_2.5.0               gplots_3.1.0           generics_0.1.0         DBI_1.1.0              pillar_1.4.6          
# [91] proxyC_0.1.5           crayon_1.3.4           KernSmooth_2.23-17     rmarkdown_2.5          viridis_0.5.1         
# [96] GetoptLong_1.0.4       progress_1.2.2         data.table_1.13.2      blob_1.2.1             digest_0.6.27         
# [101] tm_0.7-7               tidyr_1.1.2            gridGraphics_0.5-0     RcppParallel_5.0.2     munsell_0.5.0         
# [106] viridisLite_0.3.0      ggplotify_0.0.5 

      
          
