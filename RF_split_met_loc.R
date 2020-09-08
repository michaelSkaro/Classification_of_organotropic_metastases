library(tidyverse)
met_samples <- data.table::fread("/mnt/storage/mskaro1/Metastatic_Organo_Tropism/tumor_samples_annotated_progression.csv")
met_samples <- met_samples[,2:13]

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")


 
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
for(proj in projects){
  # read data in to the file. reverse sort the exp.locs column, reduce to the
  # data plus the top 5, 8 labels. Get ready for classifications. Feature selection and 
  # actual bioinformatic analysis.
  dat <- data.table::fread(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/RF_input_Metastatic/{proj}_metastatic_data_RNAseq.csv", header = TRUE)) %>%
    column_to_rownames("V1")
  
  # Reduce the samples to the top 5 seeding locations
  
  dat <- data.table::fread(str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/RF_input_Metastatic/TCGA-LUSC_metastatic_data_RNAseq.csv", header = TRUE)) %>%
    column_to_rownames("V1")
  
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
  
  write.csv2(x = dat, file = str_glue("/mnt/storage/mskaro1/PanCancerAnalysis/ML_2019/Metastatic_loci_consolidated/{proj}_metastatic_data_RNAseq.csv"))
  
  
  
  # transpose to complete DEA from normal
  
  # extract the samples that are from each loc, DEA for the tumor vs normal
  
  # feature extracts from the 
  
  
  # Enriched pathways for each of the seeding locations
  
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
  
  
  
  
