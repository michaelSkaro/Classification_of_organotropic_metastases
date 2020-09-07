library(tidyverse)
met_samples <- data.table::fread("tumor_samples_annotated_progression.csv")
met_samples <- met_samples[,2:13]

# gene annotaiton, may need, may not
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")

clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")


 
# objective: pull the samples from the files in the CSBL shared folder and
# combine them into a file for ML classification

projects <- unique(dat$project)
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



dat <- data.table::fread("TCGA-BLCA_metastatic_data_RNAseq.csv", header = TRUE) %>%
  column_to_rownames("V1")
