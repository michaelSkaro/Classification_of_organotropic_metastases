# In this investigation we will coduct these analyses:


# Differenital expression between normal tissue and primary tumors 
# Differenital expression between primary tumors that have progressed and those that have not yet. 
# Differenital expression between seeding loci from one primary location
# Differenital expression between primary tumors arising in different locations but seed in the same location

# Construct an adjacency matrix using WGCNA
# Extract the adjacency matricies. Select matricies with significant of coverage of samples and 
# use these as input matricies for learning

# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)



#projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

# add the clinical annotations of progression

# read in the files in one by one and select the correct columns. 
# concatenate the files and end with one uniform progression file. 

proj <- projects[1]

foo <- prog_file

for(proj in projects){
  
  prog_file <- data.table::fread(str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{proj}_metastatic_status.csv"), 
                                 stringsAsFactors = TRUE, header = TRUE) %>%
    dplyr::select(bcr_patient_barcode,met_loc,LymphNodeStatus,Metastatic_status)

  foo <- rbind(foo, prog_file)
  
}


foo <- foo[order(foo$bcr_patient_barcode),]

foo <- foo[!duplicated(foo$bcr_patient_barcode),]

prog_file <- foo




View(prog_file)

normal.samples <- clinical[sample_type == "Solid Tissue Normal"]
tumor.samples <- clinical[sample_type != "Solid Tissue Normal"]
tumor.samples <- tumor.samples[tumor.samples$sample_type == "Primary Tumor",]
tumor.samples <- tumor.samples[tumor.samples$project %in% projects,]

tumor.samples$bcr_patient_barcode <- substr(tumor.samples$barcode_short, 0, 12)
clinical$barcode_short <- substr(clinical$barcode, 0,16)




tumor.samples <- left_join(tumor.samples, prog_file, by = "bcr_patient_barcode")



write.csv(tumor.samples, file = str_glue("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv"))


tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")


proj <- projects[12]

for (proj in projects) {
  
  df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
    as_tibble() %>%
    tibble::column_to_rownames(var = "Ensembl")
  
  n<-dim(df.exp)[1]
  df.exp<-df.exp[1:(n-5),]
  
  coldata<- tumor.samples[tumor.samples$project == proj,]
  #coldata.n <- normal.samples[normal.samples$project == proj,]
  
  df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
  rownames(coldata) <- coldata$barcode
  coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)
  
  rownames(coldata) <- sort(rownames(coldata))
  colnames(df.exp) <- sort(colnames(df.exp))

  coldata$Metastatic_status <- as.factor(coldata$Metastatic_status)
  
  index <- is.na(coldata$Metastatic_status)
  coldata$Metastatic_status[index] <- 0
  
  
  dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ Metastatic_status)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds$Metastatic_status <- relevel(dds$Metastatic_status, ref = "0")

  dds <- DESeq(dds)
  
  save(dds, file = str_glue("{proj}_DE_met.RData"))

  #res <- results(dds)
  #resOrdered <- res[order(res$pvalue),]
  #resOrdered <- as.data.frame(resOrdered)
  #res <- as.data.frame(res)
  
  #write.csv(resOrdered, file = str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/{proj}_DE.csv"))


  res <- data.table::fread(file = str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/{proj}_DE.csv"), header = TRUE) %>%
    column_to_rownames("V1")
  
  
  
  p <- EnhancedVolcano(res, lab = rep(" ", length(rownames(res))),x = "log2FoldChange", y = "padj", xlab = bquote(~Log[2]~ "fold change"),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)), pCutoff = 10e-2, FCcutoff = 2.0, ylim=c(0,7.5), xlim = c(-3,3),
                  transcriptLabSize = 3.0, colAlpha = 0.7,legend=c("NS","Log2 FC","Adjusted p-value", "Adjusted p-value & Log2 FC"),legendPosition = "bottom",
                  legendLabSize = 10, legendIconSize = 3.0, border = "full", borderWidth = 1.5, borderColour = "black", gridlines.major = FALSE,
                  gridlines.minor = FALSE)
  
  
  ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Viz/{proj}"),"_TvsN_Volcano",".png"),
         plot = p, device = "png", width = 8, height =5, units = "in", dpi = "retina")
  
  # pathway enrichment plots
  
  
  dat <- data.table::fread(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/{proj}_DE.csv"))
  dat <- as.data.frame(dat)
  dat <- column_to_rownames(dat, "V1")
  dat$ENSEMBL <- substr(rownames(dat), 1, 15)
  
  dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
  dat <- dat[dat$padj < 1.0e-3,]
  dat <- dat[!is.na(dat$ENSEMBL),]
  
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Viz/{proj}"),"_MF_",".png"),
         plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Pathway_enrichment_flat_files/{proj}"),"_MF_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Viz/{proj}"),"_CC_",".png"),
         plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Pathway_enrichment_flat_files/{proj}"),"_CC_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  ggsave(filename=paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Viz/{proj}"),"_BP_",".png"),
         plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/Pathway_enrichment_flat_files/{proj}"),"_BP_",".txt"))
  
  
  
}
  
# in this analysis we will subset the DEG data into non-metastatic tumor vs normal
# and metastatic tissues vs normal.

# this will allow us to invesitgate if there are genes/gene sets that are specific for the expansion of the disease

# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)



#projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")


normal.samples <- clinical[sample_type == "Solid Tissue Normal"]

tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv", stringsAsFactors = TRUE) %>%
  column_to_rownames("V1")


proj <- projects[3]

for(proj in projects){
  df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"), stringsAsFactors = TRUE) %>%
    as_tibble() %>%
    tibble::column_to_rownames(var = "Ensembl")
  
  n<-dim(df.exp)[1]
  df.exp<-df.exp[1:(n-5),]
  
  coldata.t<- tumor.samples %>%
    dplyr::filter(project == proj) %>%
    dplyr::filter(Metastatic_status==0)
  
  coldata.n <- normal.samples %>%
    dplyr::filter(normal.samples$project == proj)
  
  coldata <- rbind.fill(coldata.t, coldata.n)
  rownames(coldata) <- coldata$barcode
  
  df.exp <- df.exp[ ,colnames(df.exp) %in% coldata$barcode]
  
  rownames(coldata) <- sort(rownames(coldata))
  colnames(df.exp) <- sort(colnames(df.exp))
  
  setdiff(rownames(coldata), colnames(df.exp))
  
  coldata$sample_type <- gsub(" ", "_", x = coldata$sample_type)
  
  dds <- DESeqDataSetFromMatrix(countData = df.exp, colData = coldata, design = ~ sample_type)
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
  
  dds <- DESeq(dds)
  
  save(dds, file = str_glue("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/{proj}_DE_No_met_vs_normal.RData"))
  
  res <- results(dds)
  resOrdered <- res[order(res$pvalue),]
  resOrdered <- as.data.frame(resOrdered)
  res <- as.data.frame(res)
  
  write.csv(resOrdered, file = str_glue("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/{proj}_DE_No_met_vs_normal_DE_res.csv"))
  
  rm(res)
  rm(resOrdered)
  rm(dds)
  rm(keep)
  
}

# in this section we will use the dds objects to produce the intersections of up and down
# regulated gene sets for cancer vs normal. We will then explre the up and down regulated 
# gene sets that set them apart. 
# we will cross referecnce these gene sets with the HRF sets for each cancer type.


# In this section we will conduct gene set enrichment analysis on the DEG sets.
# I like the 4 quadrant pop-art in the TCGAbiolinks viz seciton better than the 
# clusterProfiler look. 



library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)

setwd("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/res")

DE.projects<- list.files()




for(DE.file in DE.projects){
  
  dat <- data.table::fread(str_glue("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/res/{DE.file}"))
  dat <- as.data.frame(dat)
  dat <- column_to_rownames(dat, "V1")
  dat$ENSEMBL <- substr(rownames(dat), 1, 15)
  
  dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
  dat <- dat[dat$padj < 1.0e-1,]
  dat <- dat[!is.na(dat$ENSEMBL),]
  
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
  p1<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,26)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/Pathway_Enrichment/",savR,"_MF_",".png"),
         plot = p1, device = "png", width = 15, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  #write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_MF_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC", pAdjustMethod = "BH")
  p2<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,26)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/Pathway_Enrichment/",savR,"_CC_",".png"),
         plot = p2, device = "png", width = 15, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  #write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_CC_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
  p3<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,26)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/Pathway_Enrichment/",savR,"_BP_",".png"),
         plot = p3, device = "png", width = 15, height = 8, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  #write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_BP_",".txt"))
  
 
  
  
}

library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(org.Hs.eg.db)

setwd("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/res")

DE.projects<- list.files()
DE.file <- DE.projects[1]

for(DE.file in DE.projects){
  
  dat <- data.table::fread(str_glue("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/res/{DE.file}"))
  dat <- as.data.frame(dat)
  dat <- column_to_rownames(dat, "V1")
  dat$ENSEMBL <- substr(rownames(dat), 1, 15)
  
  dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
  dat <- dat[dat$padj < 1.0e-1,]
  dat <- dat[!is.na(dat$ENSEMBL),]

  #library(TCGAbiolinks)
  # Enrichment Analysis EA
  # Gene Ontology (GO) and Pathway enrichment by DEGs list
  #bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
  Genelist <- clusterProfiler::bitr(dat$ENSEMBL, OrgDb = org.Hs.eg.db ,fromType = "ENSEMBL", toType = "SYMBOL", drop = TRUE)
  
  system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist$SYMBOL))
  
  # Enrichment Analysis EA (TCGAVisualize)
  # Gene Ontology (GO) and Pathway enrichment barPlot
  savR<- substr(DE.file, 1,26)
  TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                              GOBPTab = ansEA$ResBP,
                              GOCCTab = ansEA$ResCC,
                              GOMFTab = ansEA$ResMF,
                              PathTab = ansEA$ResPat,
                              nRGTab = Genelist, 
                              nBar = 10,
                              filename=paste0("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/Pathway_Enrichment/",
                                              savR,"TCGA_viz",".pdf"))
  
  #ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/Pathway_Enrichment/",savR,"TCGA_viz",".png"),
   #      plot = p, device = "png", width = 15, height = 8, units = "in", dpi = "retina")
  
}

DE.projects<- list.files()
DE.file <- DE.projects[1]


for(DE.file in DE.projects){
  
  dat <- data.table::fread(str_glue("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/res/{DE.file}"))
  dat <- as.data.frame(dat)
  dat <- column_to_rownames(dat, "V1")
  dat$ENSEMBL <- substr(rownames(dat), 1, 15)
  
  dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
  dat <- dat[dat$padj < 1.0e-1,]
  dat <- dat[!is.na(dat$ENSEMBL),]
  
  #library(TCGAbiolinks)
  # Enrichment Analysis EA
  # Gene Ontology (GO) and Pathway enrichment by DEGs list
  #bitr(geneID, fromType, toType, OrgDb, drop = TRUE)
  Genelist <- clusterProfiler::bitr(dat$ENSEMBL, OrgDb = org.Hs.eg.db ,fromType = "ENSEMBL", toType = "SYMBOL", drop = TRUE)
  
  system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist$SYMBOL))
  
  # Enrichment Analysis EA (TCGAVisualize)
  # Gene Ontology (GO) and Pathway enrichment barPlot
  savR<- substr(DE.file, 1,26)
  TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP), 
                          GOBPTab = ansEA$ResBP,
                          GOCCTab = ansEA$ResCC,
                          GOMFTab = ansEA$ResMF,
                          PathTab = ansEA$ResPat,
                          nRGTab = Genelist, 
                          nBar = 10,
                          filename=paste0("~/storage/Metastatic_Organo_Tropism/NoMetTumor_vs_Normal/Pathway_Enrichment/",
                                          savR,"TCGA_viz",".pdf"))
  
  #ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/MetTumor_vs_Normal/Pathway_Enrichment/",savR,"TCGA_viz",".png"),
  #      plot = p, device = "png", width = 15, height = 8, units = "in", dpi = "retina")
  
}




# conduct PCA analysis on metastatic vs normal and non-metastatic vs normal

# set working directory
setwd("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/dds")

# make dds.files list 
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

proj <- projects[1]


for(proj in projects){
  # read in dds data
  
  load(str_glue("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/dds/{proj}_DE_met.RData"))
  
  # annotate
  
  dds <- estimateSizeFactors(dds)
  
  se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                             colData=colData(dds))
  
  p <- plotPCA( DESeqTransform( se ), intgroup = "Metastatic_status")
  savR<- str_glue("{proj}_Met_vsNormal_PCA")
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/metastatic_vs_non_metastatic_DGE_analysis/PCA/",savR,".pdf"),
         plot = p, device = "pdf", width = 8, height = 8, units = "in", dpi = "retina")
  
  
}  
  
  
  
  
# trying some new things :)

# install and load the packages necessary for learning

#install.packages("randomForest")
library(randomForest)





# activate the use the libraries that we will need.
library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)
library(stringr)
library(tidyverse)
library(dplyr)

library(caret)

#projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")

# add the clinical annotations of progression

# read in the files in one by one and select the correct columns. 
# concatenate the files and end with one uniform progression file. 

tumor.samples <- data.table::fread("~/storage/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv") %>% 
  tibble::column_to_rownames(var ="V1")


proj <- projects[13]
setwd("~/storage/PanCancerAnalysis/ML_2019")
for(proj in projects){

  # read in data   
  dat <- data.table::fread(str_glue("{proj}_learning.csv", header = TRUE))
  colnames(dat)[1] <- "barcode"
  
  # remove the normal samples, in this section we are only classifying the metastatic cases from 
  # here on out
  
  dat <- dat[dat$sample_type == "Primary_Tumor",]
  
  
  # add column annotation for the metastatic staus of the sample's metastatic status
  
  dat<- left_join(dat, tumor.samples, by="barcode")
  
  # get rid of some of the unnecessary columns that were added, we want the majority
  # of the data to be numeric with only a few feactorial variables.
  
  dat <- dat %>%
    dplyr::select(-barcode,-sample_type.x,-sample_type.y, -project, -tumor_stage, -caseID, -fileID, -filename, 
                  -barcode_short, -bcr_patient_barcode,-met_loc, -LymphNodeStatus)
  
  
  
  # split the datas so that we get some training and some testing shoooorrdd
  
  set.seed(2)
  
  id <- sample(2, nrow(dat), prob = c(0.8,0.2), replace =TRUE)
  
  dat.train <- dat[id ==1,]
  dat.test <- dat[id ==2,]
  
  # tune RF finds the optimal mtry for the learning. We do not want to use all of the 
  # predictor variables in data set or it may actually learn nothing. We want each tree to be different
  # and therefore more powerful

  set.seed(1234)
  
  # Define the control
  trControl <- trainControl(method = "cv",
                            number = 10,
                            search = "grid")
  
  
  
  
  
  # Search best mtry
  
  set.seed(1234)
  tuneGrid <- expand.grid(.mtry = c(1: 10))
  rf_mtry <- train(Metastatic_status~.,
                   data = data.train,
                   method = "rf",
                   metric = "Accuracy",
                   tuneGrid = tuneGrid,
                   trControl = trControl,
                   importance = TRUE,
                   nodesize = 14,
                   ntree = 1000)
  
  best_mtry <- rf_mtry$bestTune$mtry 
  
  
  # Run the model

  dat.forest <- randomForest::randomForest(Metastatic_status~., data= dat.train)

  dat.pred <- stats::predict(dat.forest, newdata = dat.test, type = "class")
  
  efficiency <- caret::confusionMatrix(table(dat.pred, dat.test$Metastatic_status))
  
  save(efficiency, file = str_glue("~/storage/Metastatic_Organo_Tropism/prediction/{proj}_random_forest_prediction_DLvsNDL.RData"))
  
  
  
}


####### add some random forests.

# random forest training sessions

#https://www.guru99.com/r-random-forest-tutorial.html



library(dplyr)
library(randomForest)
library(caret)
library(e1071)

# import the data i.e get the training sets
data_train <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/train.csv")
glimpse(data_train)
data_test <- read.csv("https://raw.githubusercontent.com/guru99-edu/R-Programming/master/test.csv") 
glimpse(data_test)

data_train <- na.exclude(data_train)
data_test <- na.exclude(data_test)

data_train$Survived <- as.factor(data_train$Survived)

# Define the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")


set.seed(1234)
# Run the model
rf_default <- train(Survived~.,
                    data = data_train,
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)
# Print the results
print(rf_default)


# Search best mtry

set.seed(1234)
tuneGrid <- expand.grid(.mtry = c(1: 10))
rf_mtry <- train(Survived~.,
                 data = data_train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
print(rf_mtry)


best_mtry <- rf_mtry$bestTune$mtry 


# Search the best maxnodes

Code explanation:
  
  # store_maxnode <- list(): The results of the model will be stored in this list
  # expand.grid(.mtry=best_mtry): Use the best value of mtry
  # for (maxnodes in c(15:25)) { ... }: Compute the model with values of maxnodes starting from 15 to 25.
  # maxnodes=maxnodes: For each iteration, maxnodes is equal to the current value of maxnodes. i.e 15, 16, 17, ...
  # key <- toString(maxnodes): Store as a string variable the value of maxnode.
  # store_maxnode[[key]] <- rf_maxnode: Save the result of the model in the list.
  # resamples(store_maxnode): Arrange the results of the model
  # summary(results_mtry): Print the summary of all the combination
  
  


store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5: 15)) {
  set.seed(1234)
  rf_maxnode <- train(Survived~.,
                      data = data_train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

# using the maxnodes

store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(20: 30)) {
  set.seed(1234)
  rf_maxnode <- train(survived~.,
                      data = data_train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  key <- toString(maxnodes)
  store_maxnode[[key]] <- rf_maxnode
}
results_node <- resamples(store_maxnode)
summary(results_node)


#Step 4) Search the best ntrees

store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  rf_maxtrees <- train(Survived~.,
                       data = data_train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 14,
                       maxnodes = 24,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)


# fit with fine tuned model

fit_rf <- train(Survived~.,
                data_train,
                method = "rf",
                metric = "Accuracy",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 14,
                ntree = 800,
                maxnodes = 24)


# make your predictions

prediction <-predict(fit_rf, data_test)


# assess the model

efficiency<- confusionMatrix(prediction, data_test$survived)















# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
# [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] candisc_0.8-0               heplots_1.3-5               car_3.0-2                   carData_3.0-2               corrplot_0.84              
# [6] wesanderson_0.3.6           limma_3.38.3                psych_1.8.12                scales_1.0.0                GEOquery_2.50.5            
# [11] hexbin_1.27.2               vsn_3.50.0                  pheatmap_1.0.12             EnhancedVolcano_1.0.1       ggrepel_0.8.0              
# [16] org.Hs.eg.db_3.7.0          AnnotationDbi_1.44.0        DESeq2_1.22.2               SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
# [21] BiocParallel_1.16.6         matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
# [26] IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0         clusterProfiler_3.10.1      dendsort_0.3.3             
# [31] ggpubr_0.2                  magrittr_1.5                WGCNA_1.66                  fastcluster_1.1.25          dynamicTreeCut_1.63-1      
# [36] forcats_0.4.0               stringr_1.4.0               dplyr_0.8.0.1               purrr_0.3.2                 readr_1.3.1                
# [41] tidyr_0.8.3                 tibble_2.1.1                ggplot2_3.1.0               tidyverse_1.2.1            
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_0.2.5       robust_0.4-18          RSQLite_2.1.1          htmlwidgets_1.3        grid_3.5.1             munsell_0.5.0         
# [7] codetools_0.2-15       preprocessCore_1.44.0  withr_2.1.2            colorspace_1.4-1       GOSemSim_2.8.0         knitr_1.22            
# [13] rstudioapi_0.10        robustbase_0.93-4      DOSE_3.8.2             urltools_1.7.2         GenomeInfoDbData_1.2.0 mnormt_1.5-5          
# [19] polyclip_1.10-0        bit64_0.9-7            farver_1.1.0           generics_0.0.2         xfun_0.5               R6_2.4.0              
# [25] doParallel_1.0.14      locfit_1.5-9.1         bitops_1.0-6           fgsea_1.8.0            gridGraphics_0.3-0     assertthat_0.2.1      
# [31] ggraph_1.0.2           nnet_7.3-12            enrichplot_1.2.0       gtable_0.3.0           affy_1.60.0            rlang_0.3.3           
# [37] genefilter_1.64.0      splines_3.5.1          lazyeval_0.2.2         acepack_1.4.1          impute_1.56.0          broom_0.5.1           
# [43] europepmc_0.3          checkmate_1.9.1        yaml_2.2.0             BiocManager_1.30.4     reshape2_1.4.3         abind_1.4-5           
# [49] modelr_0.1.4           backports_1.1.3        qvalue_2.14.1          Hmisc_4.2-0            tools_3.5.1            ggplotify_0.0.3       
# [55] affyio_1.52.0          RColorBrewer_1.1-2     ggridges_0.5.1         Rcpp_1.0.1             plyr_1.8.4             base64enc_0.1-3       
# [61] progress_1.2.0         zlibbioc_1.28.0        RCurl_1.95-4.12        prettyunits_1.0.2      rpart_4.1-13           viridis_0.5.1         
# [67] cowplot_0.9.4          haven_2.1.0            cluster_2.0.7-1        data.table_1.12.0      DO.db_2.9              openxlsx_4.1.0        
# [73] triebeard_0.3.0        mvtnorm_1.0-10         hms_0.4.2              xtable_1.8-3           XML_3.98-1.19          rio_0.5.16            
# [79] readxl_1.3.1           gridExtra_2.3          compiler_3.5.1         crayon_1.3.4           htmltools_0.3.6        pcaPP_1.9-73          
# [85] Formula_1.2-3          geneplotter_1.60.0     rrcov_1.4-7            lubridate_1.7.4        DBI_1.0.0              tweenr_1.0.1          
# [91] MASS_7.3-51.1          Matrix_1.2-15          cli_1.1.0              igraph_1.2.4           pkgconfig_2.0.2        fit.models_0.5-14     
# [97] rvcheck_0.1.3          foreign_0.8-71         xml2_1.2.0             foreach_1.4.4          annotate_1.60.1        XVector_0.22.0        
# [103] rvest_0.3.2            digest_0.6.18          cellranger_1.1.0       fastmatch_1.1-0        htmlTable_1.13.1       curl_3.3              
# [109] nlme_3.1-137           jsonlite_1.6           viridisLite_0.3.0      pillar_1.3.1           lattice_0.20-38        httr_1.4.0            
# [115] DEoptimR_1.0-8         survival_2.43-3        GO.db_3.7.0            glue_1.3.1             zip_2.0.1              UpSetR_1.3.3          
# [121] iterators_1.0.10       bit_1.1-14             ggforce_0.2.1          stringi_1.4.3          blob_1.1.1             latticeExtra_0.6-28   
# [127] memoise_1.1.0
# 
# 
