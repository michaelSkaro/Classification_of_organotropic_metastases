# metastatic tropism: This work will be on the Gene expression branch
# The leaf will be DE for tumor vs. Normal. 
# I will use the DE we have for every cancer type. I will make figures for supplement.

# what is the objective here?
# establish roles of TS genes in each cancer type. We know the drivers but 
# we can define the the roles of TS for each cancer type in carcinogenesis


# identify pathways enriched in TvsN for each primary tissue.

# we have already completed the RNA-seq DE
library(clusterProfiler)
library(ggplot2)
library(tidyverse)

setwd("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal")

DE.projects<- list.files()
DE.projects <- DE.projects[2:15]

for(DE.file in DE.projects){
  
  dat <- data.table::fread(str_glue("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal/{DE.file}"))
  dat <- as.data.frame(dat)
  dat <- column_to_rownames(dat, "V1")
  dat$ENSEMBL <- substr(rownames(dat), 1, 15)
  
  dat <- dat[abs(dat$log2FoldChange) > 1.0,] 
  dat <- dat[dat$padj < 1.0e-3,]
  dat <- dat[!is.na(dat$ENSEMBL),]
  
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,9)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/visualization/",savR,"_MF_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_MF_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,9)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/visualization/",savR,"_CC_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_CC_",".txt"))
  
  
  ego<- enrichGO(gene = dat$ENSEMBL, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH")
  p<- clusterProfiler::dotplot(ego, showCategory =20)
  savR<- substr(DE.file, 1,9)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/visualization/",savR,"_BP_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  ego_vals <- cbind(ego$ID, ego$Description, ego$GeneRatio, ego$BgRatio, ego$qvalue, ego$geneID)
  write.table(ego_vals, paste0("~/storage/Metastatic_Organo_Tropism/Pathway_Enrichment/flats/",savR,"_BP_",".txt"))
  
  
}
