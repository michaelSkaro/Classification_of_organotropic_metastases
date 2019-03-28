# generate plots(heatmap,PCA, Volcano) for all of the sequencings results for each cancer type.
library(EnhancedVolcano)
library(pheatmap)
library(DESeq2)

setwd("~/CSBL_shared/RNASeq/TCGA/DEA/tumor_vs_normal/dds")

DE.projects<- list.files()

for(DE.file in DE.projects){

  load(DE.file)

  res <- results(dds)
  
  # volcano
  
  p <- EnhancedVolcano(res,
                  
                  #lab = rownames(res),
                  
                  lab = rep(" ", length(rownames(res))),
                  
                  x = "log2FoldChange",
                  
                  y = "padj",
                  
                  xlab = bquote(~Log[2]~ "fold change"),
                  
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  
                  pCutoff = 10e-10,
                  
                  FCcutoff = 2.0,
                  
                  ylim=c(0,60),
                  
                  transcriptLabSize = 3.0,
                  
                  colAlpha = 0.7,
                  
                  legend=c("NS","Log2 FC","Adjusted p-value",
                           "Adjusted p-value & Log2 FC"),
                  
                  legendPosition = "bottom",
                  
                  legendLabSize = 10,
                  
                  legendIconSize = 3.0)
  
  savR<- substr(DE.file, 1,9)
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/Volcano/",savR,"_TvsN_Volcano_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  
  
  # heatmap
  
  vsd <- vst(dds, blind=FALSE)
  library("pheatmap")
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)
  df <- as.data.frame(colData(dds)[,c("condition")])

  p <- pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
           cluster_cols=TRUE, annotation_col=df)
   
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/Heatmaps/",savR,"_TvsN_HeatMap_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  
  # PCA plot
  p <- plotPCA(vsd, intgroup=c("condition"))
  
  ggsave(filename=paste0("~/storage/Metastatic_Organo_Tropism/PCA/",savR,"_TvsN_PCA_",".png"),
         plot = p, device = "png", width = 25, height = 25, units = "in", dpi = "retina")
  
  
  rm(dds)
  
}
