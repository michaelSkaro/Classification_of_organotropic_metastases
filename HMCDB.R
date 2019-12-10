# change of plans

# not all data is being normalized in the same way on the previous note. This will have to be streamlined.

# I will downloadand analyze raw CEL files. 



# set the gse object


library(GEOquery)

setwd("~/storage/Metastatic_Organo_Tropism")
#Metastatic Disease datasheet

refDat <- data.table::fread("Metastatic_database_project_information.csv")

# view the appropriate headers we need to iterate over

colnames(refDat)

# [1] "Experiment_id"     "Pubmed_id"         "Dataset_id"        "Platform_id"       "Standard"          "Class_id"          "Class_name"        "Sample_id"        
# [9] "Cancer_type"       "Cancer_subtype"    "Metastasis_status" "Primary_site"      "Metastasis_site"   "Sample_label" 

# use the experiment id in our i for loop, then use the datasetID and the platformID for j and k

exps <- unique(refDat$Experiment_id)
gseIDs <- unique(refDat$Dataset_id)
GSMs <- unique(refDat$Sample_id)
platforms <- unique(refDat$Platform_id)

gseEntry <- gseIDs[1] 


for(gseEntry in gseIDS){
  my.gse <- gseEntry
  
  if(!file.exists("geo_downloads")) dir.create("geo_downloads")
  if(!file.exists("results"))  dir.create("results", recursive=TRUE)
  ##get data from GEO
  my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
  
  class(my.geo.gse)
  
  length(my.geo.gse)
  
  names(my.geo.gse)
  
  # will have to change for the specific GPL we need
  
  ##get rid of list structure
  my.geo.gse <- my.geo.gse[[1]]
  
  # look at the structure 
  
  ##you can see the structure of the object
  str(my.geo.gse)
  
  # check the names 
  colnames(pData(my.geo.gse))
  
  # check the normalization methods
  
  pData(my.geo.gse)$data_processing[1]
  
  # [1] Data for both channels were Lowess-normalized and then log(2) ratio was taken
  # Levels: Data for both channels were Lowess-normalized and then log(2) ratio was taken
  
  # check the expression set
  
  head(exprs(my.geo.gse))
  
  # raw cel file
  
  if(!file.exists(paste0("./geo_downloads/",my.gse)))
    getGEOSuppFiles(my.gse, makeDirectory=T, baseDir="geo_downloads")
  
  list.files("geo_downloads")
  
  list.files(paste0("geo_downloads/",my.gse))
  
  #file.list <- read.delim(paste0("geo_downloads/",my.gse,"/filelist.txt"), as.is=T)
  
  untar(paste0("geo_downloads/",my.gse,"/",my.gse,"_RAW.tar"), exdir=paste0("geo_downloads/",my.gse,"/CEL"))
  list.files(paste0("geo_downloads/",my.gse,"/CEL"))
  
  ##make data frame of phenoData
  my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
  head(my.pdata)
  
  dim(my.pdata)
  
  colnames(my.pdata)
  
  head(my.pdata[, c("title", "geo_accession", "description")], 10)
  
  my.pdata.clean <- refDat[refDat$Dataset_id ==my.gse,]
  
  
  # make the distringuishing column the factor and make a reference
  
  design.mat <- my.pdata.clean
  
  design.mat$Sample_label <- as.character(design.mat$Sample_label)
  
  
  design.mat$sample.levels <- as.factor(design.mat$Sample_label)
  design.mat$sample.levels <- relevel(design.mat$sample.levels, ref="Primary Normal")
  
  
  table(design.mat$sample.levels)
  
  # Visualize the density plots of the matrix. Show data has been normalized
  level.pal <- brewer.pal(4, "Dark2")
  level.cols <- level.pal[unname(design.mat$sample.levels)]
  
  plotDensities(exprs(my.geo.gse), legend=F, col=level.cols , main="Arrays Normalized")
  legend("topright", legend=levels(design.mat$sample.levels), fill=level.pal)
  
  
  # hierarchical clustering
  
  
  
  cluster.dat <- exprs(my.geo.gse)
  gene.mean <- apply(cluster.dat, 1, mean)
  gene.sd <- apply(cluster.dat, 1, sd)
  cluster.dat <- sweep(cluster.dat, 1, gene.mean, "-")
  cluster.dat <- sweep(cluster.dat, 1, gene.sd, "/")
  
  my.dist <- dist(t(cluster.dat), method="euclidean")
  my.hclust <- hclust(my.dist, method="average")
  my.hclust$labels <- design.mat$sample.levels
  plot(my.hclust$merge, cex=0.75, main="Comparison of Biological Replicates", xlab="Euclidean Distance")
  
  #Correlation matrix
  dat <- as.matrix(exprs(my.geo.gse)) 
  my.cor <- cor(dat)
  my.dist <- as.dist(1 - my.cor)
  my.hclust <- hclust(my.dist, method="average")
  my.hclust$labels <- design.mat$sample.levels
  plot(my.hclust, cex=0.75, main="Comparison of Biological Replicates")
  
  
  # heatmaps
  
  library(ComplexHeatmap)
  
  ComplexHeatmap::Heatmap(my.cor, name = "Correlation_of_expresssion", split =2, show_column_dend = FALSE)
  
  ComplexHeatmap::Heatmap(dat, name = "normalized expression", show_row_names = FALSE)
  
  
  
}





