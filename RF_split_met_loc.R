library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
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

for(proj in projects){
  # iterate over the projects to extract the counts files from each of the 
  # respective cancer types. Only pull people with a clinically identified 
  # distant lesion.
  
  # read in project and match th samples
  
  dat <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv"))
  transcript_ids <- dat$Ensembl
  
  proj_met_samples <- met_samples[met_samples$project == proj,]
  proj_met_samples <- proj_met_samples[proj_met_samples$met_loc != "NA",]
  dat <- dat %>% select(one_of(proj_met_samples$barcode))
  dat <- cbind(transcript_ids,dat)
  
  n<-dim(dat)[1]
  dat<-dat[1:(n-4),]
  
  #transpose the data so it can be learned on properly.
  #annotate the data frame with the correct seeding locations.
  dat <- as.data.frame(t(dat))
  colnames(dat) <- dat[1,]
  n<-dim(dat)[1]
  dat<-dat[2:n,]
  dat$barcode <- rownames(dat)
  dat<- left_join(dat, proj_met_samples, by="barcode")
  
  
  
  
  #expand the comma separated column met_loc and add samples for each 
  #one of the unique locations
  
  foo <- dat$barcode
  bar <- dat$met_loc
  
  foobar <- as.data.frame(cbind(foo,bar))
  s <- strsplit(foobar$bar, split = ",")
  foobar <- data.frame(foo = rep(foobar$foo, sapply(s, length)), bar = unlist(s))
  
  colnames(foobar) <- c("barcode","exp.locs")
  
  dat2 <- left_join(foobar,dat, by ="barcode")
  
  write.csv2(x = dat2, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/RF_input_Metastatic/{proj}_metastatic_data_RNAseq.csv"))
  
  
}

# make reduced data sets with some of the major seeding areas. 
# This will allow fo feature selection for each of the areas 
# possible options for better classifications with reduced noise.

setwd('/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/RF_input_Metastatic')
projects <- unique(met_samples$project)
projects <- sort(projects)
proj <- projects[2]
projects <- c(projects[1:3], projects[7:9])
for(proj in projects){
  # read data in to the file. reverse sort the exp.locs column, reduce to the
  # data plus the top 5, 8 labels. Get ready for classifications. Feature selection and 
  # actual bioinformatic analysis.
  dat <- data.table::fread(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/RF_input_Metastatic/{proj}_metastatic_data_RNAseq.csv", header = TRUE)) %>%
    column_to_rownames("V1")
  
  # Reduce the samples to the top 5 seeding locations

  common_seeding_locations <- as.data.frame(sort(table(dat$exp.locs), decreasing = TRUE))
  colnames(common_seeding_locations) <- c("Location", "Freq")
  names(dat)[names(dat) == 'exp.locs'] <- 'labels'
  n<-dim(dat)[2]
  dat<-dat[,1:(n-12)]
  
  index <- dat$labels == " Brain" 
  dat$labels[index] <- "Brain"
  
  index <- dat$labels == " Bone" 
  dat$labels[index] <- "Bone"
  
  index <- dat$labels == " Brain" 
  dat$labels[index] <- "Brain"
  
  index <- dat$labels == " Lung" 
  dat$labels[index] <- "Lung"
  
  index <- dat$labels == " Liver" 
  dat$labels[index] <- "Liver"
  
  index <- dat$labels == " Lymph Node" 
  dat$labels[index] <- "Lymph Node"
  
  index <- dat$labels == " Lymph  Node" 
  dat$labels[index] <- "Lymph Node"
  
  index <- dat$labels == " Breast" 
  dat$labels[index] <- "Breast"
  
  index <- dat$labels == " Skin" 
  dat$labels[index] <- "Skin"
  
  index <- dat$labels == " Adrenal Gland" 
  dat$labels[index] <- "Adrenal Gland"
  
  index <- dat$labels == " Chest Wall" 
  dat$labels[index] <- "Chest"
  
  index <- dat$labels == " Colon" 
  dat$labels[index] <- "Colon"
  
  index <- dat$labels == " Omentum" 
  dat$labels[index] <- "Omentum"
  
  common_seeding_locations <- common_seeding_locations[common_seeding_locations$Freq >=8,]
  
  dat <- as.data.frame(dat[dat$labels %in% common_seeding_locations$Location,])
  
  # now we need to one hot encode the samples so we can do proper classification
  
  dat$seen =1
  
  index <- dat$labels == "" 
  dat$labels[index] <- "Unknown"
  
  
  dat <- dat %>%
    pivot_wider(names_from = labels, values_from = seen, values_fill =0)
  
  # count how many times the unique onehot encoding happens. We need to eliminate 
  # ones that do not occur atleast once
  y=dat %>% select(Breast:Lung) %>% group_by_all() %>% summarise(COUNT = n())
  
  res <- dat[duplicated(dat %>% select(Colon:Liver)),]
  d <- duplicated(dat %>% select(Colon:Liver)) | duplicated(dat %>% select(Colon:Liver), fromLast = TRUE)
  
  dat2 <- dat[d,]
  
  write.csv(x = dat2, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/{proj}_metastatic_data_RNAseq.csv"))
  
  
  # transpose to complete DEA from normal
  
  # extract the samples that are from each loc, DEA for the tumor vs normal
  
  # feature extracts from the random forest. Compare the DEA vs. selected features
  
  # Enriched pathways for each of the seeding locations
  
} 
setwd("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels")
for(proj in projects){
  # read in file that contains the metastatic cases
  if(proj =="TCGA-BLCA"){
    df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
      as_tibble() %>%
      tibble::column_to_rownames(var = "Ensembl")
    dat <- data.table::fread("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/one_hot_encoded_labels/TCGA-BLCA_metastatic_data_RNAseq.csv", header = TRUE) %>%
      column_to_rownames("V1")
    # cut the columns off to make a coldata draft
    coldata <- dat %>% dplyr::select(c(barcode,Bone:Pelvis))
    coldata.all.mets <- coldata
    dat <- dat %>% dplyr::select(-c(Bone:Pelvis))
    dat <- as.data.frame(t(dat))
    names(dat) <- dat[1,]
    dat <- as.data.frame(dat[2:length(dat[,1]),])
  
    coldata.n <- normal.samples[normal.samples$project == proj,]
    coldata.n <- coldata.n %>% dplyr::select(barcode,sample_type)
    
    
    
    # Bone
    
    coldata.bone <- dplyr::select(coldata,  c(barcode, Bone)) %>%
      dplyr::filter(Bone ==1)
    coldata.bone$Bone <-  "Advanced_to_Bone"
    coldata.n$sample_type <-  "Solid_Tissue_Normal"
    colnames(coldata.bone) <- c("barcode", "sample_type")
    
    coldata.exp <- rbind(coldata.n, coldata.bone) 
    
    
    rownames(coldata.exp) <- coldata.exp$barcode
    
    df.exp <- df.exp %>%
      dplyr::select(coldata.exp$barcode)
    
    # feed coldata.exp and df.exp to the differential experssion analysis with DEseq2
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = df.exp,colData = coldata.exp, design = ~sample_type)
     
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")

    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$pvalue),]
    resOrdered <- as.data.frame(resOrdered)
    res <- as.data.frame(res)
   
    write.csv(resOrdered, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_Bone_DE.csv"))
    
    save(dds, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_Bone_DE_met.RData"))
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    res <- res[abs(res$log2FoldChange) > 1.0,] 
    res <- res[res$padj < 1.0e-3,]
    res <- res[!is.na(res$ENSEMBL),]
    
    p <- EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlim = c(-5, 8))
    
    p <- EnhancedVolcano(res, lab = rep(" ", length(rownames(res))),x = "log2FoldChange", 
                         y = "padj", xlab = bquote(~Log[2]~ "fold change"), ylab = bquote(~-Log[10]~adjusted~italic(P)), 
                         pCutoff = 10e-2, FCcutoff = 2.0, ylim=c(0,7.5), xlim = c(-3,3), colAlpha = 0.7,
                         legend=c("NS","Log2 FC","Adjusted p-value", "Adjusted p-value & Log2 FC"),legendPosition = "bottom",
                         legendLabSize = 10, legendIconSize = 3.0, border = "full", borderWidth = 1.5, borderColour = "black", 
                         gridlines.major = FALSE,gridlines.minor = FALSE)
    
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    
    ego<- enrichGO(gene = res$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
    p<- clusterProfiler::dotplot(ego, showCategory =20)
    p
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Bone_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Bone_flat_",".txt"))
    
    
    
    
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_CC_Bone_Flat_",".txt"))
    
    
    ego<- enrichGO(gene = res$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
    p<- clusterProfiler::dotplot(ego, showCategory =20)
    p
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_BP_Bone_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_BP_Bone_Flat_",".txt"))
    
    
    rm(res,dds,ego,df.exp,resOrdered,coldata.exp, coldata.bone)
    
    # Prostate
    
    df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
      as_tibble() %>%
      tibble::column_to_rownames(var = "Ensembl")
    
    coldata.prostate <- coldata.all.mets %>% dplyr::select(c(barcode, Prostate))
    coldata.prostate$Prostate <-  "Advanced_to_Prostate"
    coldata.n$sample_type <-  "Solid_Tissue_Normal"
    colnames(coldata.prostate) <- c("barcode", "sample_type")
    coldata.exp <- rbind(coldata.n, coldata.prostate) 
    rownames(coldata.exp) <- coldata.exp$barcode
    
    df.exp <- df.exp %>%
      dplyr::select(coldata.exp$barcode)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = df.exp,colData = coldata.exp, design = ~sample_type)
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
    
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$pvalue),]
    resOrdered <- as.data.frame(resOrdered)
    res <- as.data.frame(res)
    
    write.csv(resOrdered, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_prostate_DE.csv"))
    
    save(dds, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_prostate_DE_met.RData"))
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    res <- res[abs(res$log2FoldChange) > 1.5,] 
    res <- res[res$padj < 1.0e-4,]
    res <- res[!is.na(res$ENSEMBL),]
    
    p <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         xlim = c(-5, 8))
    
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Prostate_volcano_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    
    ego<- enrichGO(gene = res$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
    p<- clusterProfiler::dotplot(ego, showCategory =20)
    p
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_prostate_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_prostate_flat_",".txt"))
    
    
    
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_CC_prostate_Flat_",".txt"))
    
    
    ego<- enrichGO(gene = res$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
    p<- clusterProfiler::dotplot(ego, showCategory =20)
    p
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_BP_prostate_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_BP_prostate_Flat_",".txt"))
    
    
    rm(res,dds,ego,df.exp,resOrdered,coldata.exp, coldata.prostate,ego_vals)
    
    #Lung
    
    df.exp <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/counts/{proj}.counts.csv")) %>%
      as_tibble() %>%
      tibble::column_to_rownames(var = "Ensembl")
    
    coldata.Lung <- coldata.all.mets %>% dplyr::select(c(barcode, Lung))
    coldata.Lung$Lung <-  "Advanced_to_Lung"
    coldata.n$sample_type <-  "Solid_Tissue_Normal"
    colnames(coldata.Lung) <- c("barcode", "sample_type")
    coldata.exp <- rbind(coldata.n, coldata.Lung) 
    rownames(coldata.exp) <- coldata.exp$barcode
    
    df.exp <- df.exp %>%
      dplyr::select(coldata.exp$barcode)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = df.exp,colData = coldata.exp, design = ~sample_type)
    
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    
    dds$sample_type <- relevel(dds$sample_type, ref = "Solid_Tissue_Normal")
    
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$pvalue),]
    resOrdered <- as.data.frame(resOrdered)
    res <- as.data.frame(res)
    
    write.csv(resOrdered, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_Lung_DE.csv"))
    
    save(dds, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}_Lung_DE_met.RData"))
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    res <- res[abs(res$log2FoldChange) > 1.5,] 
    res <- res[res$padj < 1.0e-4,]
    res <- res[!is.na(res$ENSEMBL),]
    
    p <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         xlim = c(-5, 8))
    
    
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Lung_volcano_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    
    res$ENSEMBL <- substr(rownames(res), 1, 15)
    
    ego<- enrichGO(gene = res$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
    p<- clusterProfiler::dotplot(ego, showCategory =20)
    p
    ggsave(filename=paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Lung_",".png"),
           plot = p, device = "png", width = 8, height = 8, units = "in", dpi = "retina")
    ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
    write.table(ego_vals, paste0(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/Metastatic_cancers_by_location_vs_Normal/{proj}"),"_MF_Lung_flat_",".txt"))
    
    
    
    
  
  }
  # get healthy data
  
  
  

  
  
  
  
  # transpose
  
  # use clinical files to pull healthy
  
  # attach healthy
  
  # DeSeq2 each seeding location against the metastatic loci
  
  # GSEA / Cluster profiler for each location
  
  # output them in a style that will allow for them to be publishable
  
  # possibly compare features in the One.vs.Rest RF feature selection (sorted by vals)
  # enrichment to the enrichment of the DEseq2 DEGs
  
}


  
  #For:
  
  # Bone 
    # Make pretty confusion plot or combine data from early classification to make all
    # one figure?
    # Bone seeding vs. Normal DEseq2
    # Upregulated vs. Down Regulated pathways
    # Draw features from random forest

  # Brain
  
  # Lung
  
  # Liver 
  
  # Lymph
  
  
  
  




