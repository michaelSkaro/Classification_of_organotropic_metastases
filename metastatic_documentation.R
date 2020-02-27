# I would liek to start by saying before anyone moves forward with this code:

# I am neither proud nor am I happy with 99% of this code. However,
# I really have exhausted many of my options in the forms of REGEX to 
# complete this work. 
# The files I was given from the database were very terribly organized and for that
# I hope no one ever tries to use this code. 

#Thanks. 










# Clear your workspace
rm(list=ls(all=TRUE))
# Load some things
library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(gsubfn)
library(dplyr)

# make a working directory to read all the files from
setwd("~/CSBL_shared/clinical/TCGA_xml") # most comprehensive


# lets decide which data sets we are going to work on
projects <- c("TCGA-BLCA","TCGA-BRCA","TCGA-COAD","TCGA-ESCA","TCGA-HNSC","TCGA-KICH","TCGA-KIRC","TCGA-KIRP","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC","TCGA-PRAD","TCGA-STAD","TCGA-THCA")

#projects <- c("BLCA","BRCA","COAD","HNSC","KICH","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","STAD","THCA")
# annotate weird gene symbols and get some clinical data
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv")
refDat <- data.table::fread("~/storage/Metastatic_Organo_Tropism/Metastatic_database_project_information.csv")
refDat <- refDat[order(refDat$Sample_id),]
refDat<- refDat[match(unique(refDat$Sample_id), refDat$Sample_id),]
clinical <- data.table::fread(
  "~/CSBL_shared/RNASeq/TCGA/annotation/counts_annotation.csv")
clinical$Sample_id <- substr(clinical$barcode, 0,16)





## define a helper function
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) ## since ifelse wont work with factors
  ifelse(as.character(x)!="", x, NA)
}



i <- projects[1]

# do the things

for(i in projects){
  
  if(i== "TCGA-BLCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    BLCA_met <- as.data.frame(dat %>% dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site, metastatic_site,`metastatic_site[1]`,
                                                    `metastatic_site[2]`,`metastatic_site[3]`,new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text,
                                                    number_of_lymphnodes_positive_by_he)) %>%
      mutate(BLCA_met, LymphNodeStatus = ifelse(is.na(number_of_lymphnodes_positive_by_he), 0,
                                                ifelse(number_of_lymphnodes_positive_by_he == 0, 0, 1))) %>%
      tidyr::unite(met_loc, other_malignancy_anatomic_site:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE, sep = ",") %>%
      mutate_each(funs(empty_as_na))
      
      
      index <- BLCA_met$met_loc == ",,,,,," 
      BLCA_met$met_loc[index] <- NA
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,,,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,,", ",")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      
      index <- BLCA_met$met_loc == ",None," 
      BLCA_met$met_loc[index] <- NA
      BLCA_met$met_loc <- stringr::str_replace(BLCA_met$met_loc, ",,", "")
      BLCA_met$met_loc <- stringr::str_replace(BLCA_met$met_loc, ",", "")
      
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Other specify", "")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Other, specify", "")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Other,", "")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "[|]", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, ",,", ",")
      BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "None,", "")
      
      BLCA_met <- BLCA_met %>%
        mutate_each(funs(empty_as_na)) %>%
        mutate(BLCA_met, Metastatic_status = ifelse(is.na(met_loc) | malignancy_type == "Prior Malignancy" & LymphNodeStatus ==0, 0, 1))

    
    index <- is.na(BLCA_met$number_of_lymphnodes_positive_by_he)
    BLCA_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    index <- BLCA_met$met_loc == "Lymph node only," 
    BLCA_met$Metastatic_status[index] <- 0
    BLCA_met$met_loc[index] <- "Lymph node"
    
    index <- BLCA_met$met_loc == "Lymph Node Only," 
    BLCA_met$Metastatic_status[index] <- 0
    BLCA_met$met_loc[index] <- "Lymph node"
    
    index <- BLCA_met$met_loc == "Bone," 
    BLCA_met$met_loc[index] <- "Bone"
    BLCA_met$Metastatic_status[index] <- 1
    
    
    index <- BLCA_met$met_loc == "Bladder," 
    BLCA_met$met_loc[index] <- "Bladder"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == "Liver," 
    BLCA_met$met_loc[index] <- "Liver"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == "Lung," 
    BLCA_met$met_loc[index] <- "Lung"
    BLCA_met$Metastatic_status[index] <- 1
    
    BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Lymph Node Only", "Lymph Node")
    BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Lymph node only", "Lymph Node")
    BLCA_met$met_loc <- stringr::str_replace_all(BLCA_met$met_loc, "Lymph node", "Lymph Node")
    
    
    
    
    index <- is.na(BLCA_met$Metastatic_status)
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Abdominal wall" 
    BLCA_met$met_loc[index] <- "Abdominal wall"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Adrenal" 
    BLCA_met$met_loc[index] <- "Adrenal"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",brain and spinal cord" 
    BLCA_met$met_loc[index] <- "brain,spinal cord"
    BLCA_met$Metastatic_status[index] <- 1
                                                      
    index <- BLCA_met$met_loc == ",brain and spinal cord" 
    BLCA_met$met_loc[index] <- "brain,spinal cord"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Groin" 
    BLCA_met$met_loc[index] <- "Groin"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Left posterior pelvic wall" 
    BLCA_met$met_loc[index] <- "Pelvis"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Lung,Brain" 
    BLCA_met$met_loc[index] <- "Lung,Brain"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Lymph Node,Bone,neck node,neck lymph node" 
    BLCA_met$met_loc[index] <- "Lymph Node,Bone,Head Neck"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Pelvic" 
    BLCA_met$met_loc[index] <- "Pelvis"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",penile urethra" 
    BLCA_met$met_loc[index] <- "Penis"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",peritoneum" 
    BLCA_met$met_loc[index] <- "Peritoneum"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",peritoneum" 
    BLCA_met$met_loc[index] <- "Peritoneum"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",psoas" 
    BLCA_met$met_loc[index] <- "Psoas"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",Rectum" 
    BLCA_met$met_loc[index] <- "Rectum"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",right hemipelvis" 
    BLCA_met$met_loc[index] <- "Pelvis"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",sacrum" 
    BLCA_met$met_loc[index] <- "Sacrum"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",sigmoid colon" 
    BLCA_met$met_loc[index] <- "Colon"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == ",soft tissue pelvic mass" 
    BLCA_met$met_loc[index] <- "Pelvis"
    BLCA_met$Metastatic_status[index] <- 1
    
    BLCA_met$met_loc <- gsub(",$", "", BLCA_met$met_loc)
    BLCA_met$met_loc <- gsub("$,", "", BLCA_met$met_loc)
    
    index <- BLCA_met$met_loc == ",Pelvis" 
    BLCA_met$met_loc[index] <- "Pelvis"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == "Abdominal wall" 
    BLCA_met$met_loc[index] <- "Abdomen"
    BLCA_met$Metastatic_status[index] <- 1
    
    index <- BLCA_met$met_loc == "Ankle" 
    BLCA_met$met_loc[index] <- "Bone"
    
    index <- BLCA_met$met_loc == "Bladder,Bone,Lung,Liver,Lymph Node,Adrenal Metastases,Pelvic Sidewall Soft Tissue Recurrence" 
    BLCA_met$met_loc[index] <- "Bladder,Bone,Lung,Liver,Lymph Node,Adrenal, Pevis, Soft Tissue"
    
    index <- BLCA_met$met_loc == "Bladder,Lung,suprapubic area" 
    BLCA_met$met_loc[index] <- "Bladder,Lung,Abdomen"
   
    index <- BLCA_met$met_loc == "Bone,Lung,Liver,pancreas" 
    BLCA_met$met_loc[index] <- "Bone,Lung,Liver,Pancreas"
    
    index <- BLCA_met$met_loc == "Fallopian tube" 
    BLCA_met$met_loc[index] <- "Ovary"
    
    index <- BLCA_met$met_loc == "Forehead,Prostate" 
    BLCA_met$met_loc[index] <- "Head and Neck, Skin"
    
    index <- BLCA_met$met_loc == "Head - face or neck NOS" 
    BLCA_met$met_loc[index] <- "Head and Neck, Skin"
    
    index <- BLCA_met$met_loc == "Liver,Lung,Bone,adrenal" 
    BLCA_met$met_loc[index] <- "Liver,Lung,Bone,Adrenal"
   
    index <- BLCA_met$met_loc == "Liver,Pelvic Soft Tissue Mass" 
    BLCA_met$met_loc[index] <- "Liver,Pelvis, Soft Tissue"
    
    index <- BLCA_met$met_loc == "lower abdominal wall and right pelvis" 
    BLCA_met$met_loc[index] <- "Abdomen, Pelvis"
    
    index <- BLCA_met$met_loc == "Lung,adjacent to the right iliac crest" 
    BLCA_met$met_loc[index] <- "Lung, Illiac crest"
    
    index <- BLCA_met$met_loc == "Lung,local & retroperitonium" 
    BLCA_met$met_loc[index] <- "Lung, Abdomen"
    
    index <- BLCA_met$met_loc == "Lung,Lung" 
    BLCA_met$met_loc[index] <- "Lung"
    
    index <- BLCA_met$met_loc == "Lung,Lung,Bladder" 
    BLCA_met$met_loc[index] <- "Lung, Bladder"
    
    index <- BLCA_met$met_loc == "Lung,Lymph Node,Bone,Soft tissue,Soft Tissue" 
    BLCA_met$met_loc[index] <- "Lung,Lymph Node,Bone,Soft tissue"
    
    index <- BLCA_met$met_loc == "Lung,Lymph Node,Lung" 
    BLCA_met$met_loc[index] <- "Lung, Lymph Node"
    
    index <- BLCA_met$met_loc == "Lung,Peritoneum- Pelvic Lymph Nodes" 
    BLCA_met$met_loc[index] <- "Lung, Lymph Node"
    
    index <- BLCA_met$met_loc == "Lung,Rectal wall,Rectal Wall" 
    BLCA_met$met_loc[index] <- "Lung, Rectum"
    
    index <- BLCA_met$met_loc == "Lung,soft tissue pelvic mass" 
    BLCA_met$met_loc[index] <- "Lung, Pelvis, Soft Tissue"
    
    index <- BLCA_met$met_loc == "Lymph Node,Bladder,Lymph Node" 
    BLCA_met$met_loc[index] <- "Bladder, Lymph Node"
    
    index <- BLCA_met$met_loc == "Lymph Node,Bone,Head Neck" 
    BLCA_met$met_loc[index] <- "Lymph Node,Bone,Head and Neck"
    
    index <- BLCA_met$met_loc == "Lymph Node,Bone,Liver,Lymph Nodes" 
    BLCA_met$met_loc[index] <- "Lymph Node,Bone,Liver"
    
    index <- BLCA_met$met_loc == "Lymph Node,Liver,Bone,Lymph Node" 
    BLCA_met$met_loc[index] <- "Lymph Node,Liver,Bone"
    
    index <- BLCA_met$met_loc == "Lymph Node,lumbar spine and right hemipelvis" 
    BLCA_met$met_loc[index] <- "Lymph Node, Spine Bone"
    
    index <- BLCA_met$met_loc == "Lymph Node,Lung,Bone,Lymph Node" 
    BLCA_met$met_loc[index] <- "Lymph  Node, Lung"
    
    index <- BLCA_met$met_loc == "Lymph Node,Lung,Liver,Lymph Node,Lung,Liver" 
    BLCA_met$met_loc[index] <- "Lymph Node,Lung,Liver"
    
    index <- BLCA_met$met_loc == "Lymph Node,Lung,Lymph Node" 
    BLCA_met$met_loc[index] <- "Lymph  Node, Lung"
    
    index <- BLCA_met$met_loc == "Lymph Node,Lymph Node" 
    BLCA_met$met_loc[index] <- "Lymph Node"
    
    index <- BLCA_met$met_loc == "Lymph Node,pelvic mass with adherence to pelvic sidewall and sigmoid colon,sigmoid colon,pelvic side wall" 
    BLCA_met$met_loc[index] <- "Lymph Node,Pelvis, Colon"
    
    index <- BLCA_met$met_loc == "Lymph Node,Scrotum and Perineum" 
    BLCA_met$met_loc[index] <- "Lymph Node,Scrotum, Perineum"
    
    index <- BLCA_met$met_loc == "Lymph Node(s)" 
    BLCA_met$met_loc[index] <- "Lymph Node"
    
    index <- BLCA_met$met_loc == "Other" 
    BLCA_met$met_loc[index] <- NA
    
    index <- BLCA_met$met_loc == "Other,Prostate" 
    BLCA_met$met_loc[index] <- "Prostate"
    
    index <- BLCA_met$met_loc == "pelvis" 
    BLCA_met$met_loc[index] <- "Pelvis"
    
    index <- BLCA_met$met_loc == "Pelvis and vagina" 
    BLCA_met$met_loc[index] <- "Pelvis, Vagina"
    
    index <- BLCA_met$met_loc == "Prostate,Other,Bladder" 
    BLCA_met$met_loc[index] <- "Prostate, Bladder"
    
    index <- BLCA_met$met_loc == "Renal Pelvis,Lymph Node,pelvis" 
    BLCA_met$met_loc[index] <- "Renal Pelvis,Lymph Node,Pelvis"
    
    index <- BLCA_met$met_loc == "small bowel" 
    BLCA_met$met_loc[index] <- "Small Bowel"
    
    index <- BLCA_met$met_loc == "Thyroid gland" 
    BLCA_met$met_loc[index] <- "Thyroid"
    
    index <- BLCA_met$met_loc == "Tongue Base of tongue" 
    BLCA_met$met_loc[index] <- "Oral Cavity"
    
    index <- BLCA_met$met_loc == "Urethra,Lymph Node,Blood-Large B Cell Lymphoma" 
    BLCA_met$met_loc[index] <- "Urethra,Lymph Node,Lymphoma"
    
    
    
    write.csv(BLCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_status.csv"))
    
  
    print("BLCA Done")
    rm(BLCA_met)
    rm(dat)
  }
  
    if(i = "TCGA-BRCA"){
      
      dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv")))
    
      BRCA_met <- as.data.frame(dat %>%
                                  dplyr::select(bcr_patient_barcode, metastatic_site_at_diagnosis, number_of_lymphnodes_positive_by_he,
                                               `metastatic_site_at_diagnosis[1]`,`metastatic_site_at_diagnosis[2]`,
                                               `metastatic_site_at_diagnosis[3]`,`metastatic_site_at_diagnosis[4]`,
                                                new_neoplasm_event_occurrence_anatomic_site, 
                                                new_neoplasm_occurrence_anatomic_site_text))%>% 
        mutate_each(funs(empty_as_na)) %>%
        mutate(BRCA_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
        tidyr::unite(met_site_merge, `metastatic_site_at_diagnosis[1]`:`metastatic_site_at_diagnosis[4]`, na.rm =TRUE, sep = ",") %>%
        tidyr::unite(neoplasm_and_distant_met, met_site_merge:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
        dplyr::select(bcr_patient_barcode,metastatic_site_at_diagnosis,neoplasm_and_distant_met,number_of_lymphnodes_positive_by_he,LymphNodeStatus) %>%
        tidyr::unite(met_loc, metastatic_site_at_diagnosis:neoplasm_and_distant_met,na.rm =TRUE,sep = ",") %>%
        mutate_each(funs(empty_as_na))%>%
        mutate(BRCA_met, Metastatic_status = ifelse(is.na(met_loc), 0, 1))
        
      index <- is.na(BRCA_met$LymphNodeStatus)
      BRCA_met$LymphNodeStatus[index] <- 0
      
      #lymph node status
      
      index <- is.na(BRCA_met$number_of_lymphnodes_positive_by_he)
      BRCA_met$number_of_lymphnodes_positive_by_he[index] <- 0
      
      
      
      # consolidate columns
      
      
      index	<-	BRCA_met$met_loc	==	",Bone"
      BRCA_met$met_loc[index]	<-	"Bone"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Brain"
      BRCA_met$met_loc[index]	<-	"Bone, Brain"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Liver"
      BRCA_met$met_loc[index]	<-	"Bone, Liver"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Liver,Liver"
      BRCA_met$met_loc[index]	<-	"Bone, Liver"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Liver|Other, specify,Breast"
      BRCA_met$met_loc[index]	<-	"Bone, Liver, Breast"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Liver|Other, specify|Brain,Intrathoracic Lymph Node"
      BRCA_met$met_loc[index]	<-	"Bone, Liver, Breast, Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Bone|Other, specify,omentum"
      BRCA_met$met_loc[index]	<-	"Bone, Omentum"
      
      index	<-	BRCA_met$met_loc	==	",Brain"
      BRCA_met$met_loc[index]	<-	"Brain"
      
      index	<-	BRCA_met$met_loc	==	",Liver"
      BRCA_met$met_loc[index]	<-	"Liver"
      
      index	<-	BRCA_met$met_loc	==	",Lung"
      BRCA_met$met_loc[index]	<-	"Lung"
      
      index	<-	BRCA_met$met_loc	==	",Lung,Liver"
      BRCA_met$met_loc[index]	<-	"Lung, Liver"
      
      index	<-	BRCA_met$met_loc	==	",Lung|Other, specify|Brain,Mediastinal and Supraclavicular Lymph Nodes"
      BRCA_met$met_loc[index]	<-	"Lung, Brain, Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Bone Marrow"
      BRCA_met$met_loc[index]	<-	"Bone"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,bone, brain"
      BRCA_met$met_loc[index]	<-	"Bone, Brain"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Breast"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,breast|Lymph Node"
      BRCA_met$met_loc[index]	<-	"Breast, Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Chest Wall"
      BRCA_met$met_loc[index]	<-	"Chest Wall"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,colon"
      BRCA_met$met_loc[index]	<-	"Colon"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Contralateral Breast"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Dermis and epidermis"
      BRCA_met$met_loc[index]	<-	"Skin"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Endometrial"
      BRCA_met$met_loc[index]	<-	"Uterus"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Endometrium"
      BRCA_met$met_loc[index]	<-	"Uterus"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Left axilla"
      BRCA_met$met_loc[index]	<-	"Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,left breast"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,LEFT BREAST"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Left Cervical Lymph Node"
      BRCA_met$met_loc[index]	<-	"Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Left Chest Wall"
      BRCA_met$met_loc[index]	<-	"Chest Wall"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Liver and Pleura and Bone"
      BRCA_met$met_loc[index]	<-	"Liver, Bone, Lung"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Lung"
      BRCA_met$met_loc[index]	<-	"Lung"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Lung, Bone, Liver, Brain and Skin Nodules"
      BRCA_met$met_loc[index]	<-	"Lung, Liver, Brain, Bone, Skin"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,malignant melanoma|Malignant melanoma"
      BRCA_met$met_loc[index]	<-	"Skin"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,mediastinal lymph node"
      BRCA_met$met_loc[index]	<-	"Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,mediastinal lymph nodes"
      BRCA_met$met_loc[index]	<-	"Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Pectoral muscle"
      BRCA_met$met_loc[index]	<-	"Muscle"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Rectum"
      BRCA_met$met_loc[index]	<-	"Rectum"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Renal"
      BRCA_met$met_loc[index]	<-	"Renal"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Right Breast"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,right breast cancer contralateral"
      BRCA_met$met_loc[index]	<-	"Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Skin and bone"
      BRCA_met$met_loc[index]	<-	"Bone, Colon"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Skin left chest wall"
      BRCA_met$met_loc[index]	<-	"Skin, Chest Wall"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Skin lesion-Basal Cell Left Lower Lateral Back"
      BRCA_met$met_loc[index]	<-	"Skin"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify,Skin, right leg"
      BRCA_met$met_loc[index]	<-	"Skin"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify|Bone,Breast recurrence|Chest wall|Breast|Chest wall, Breast Recurrence"
      BRCA_met$met_loc[index]	<-	"Bone, Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify|Bone,Left Breast"
      BRCA_met$met_loc[index]	<-	"Bone, Breast"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify|Bone|Lung|Liver,Axilla"
      BRCA_met$met_loc[index]	<-	"Bone, Liver, Lung"
      
      index	<-	BRCA_met$met_loc	==	",Other, specify|Lung|Liver|Bone,chestwall|Adrenal glands|Liver, bone, chest wall, adrenal glands"
      BRCA_met$met_loc[index]	<-	"Lung, Liver, Bone, Adrenal Gland"
      
      index	<-	BRCA_met$met_loc	==	"Bone,"
      BRCA_met$met_loc[index]	<-	"Bone"
      
      index	<-	BRCA_met$met_loc	==	"Bone,,Bone"
      BRCA_met$met_loc[index]	<-	"Bone"
      
      index	<-	BRCA_met$met_loc	==	"Bone,,Bone|Other, specify|Lung,Intrathoracic lymph node"
      BRCA_met$met_loc[index]	<-	"Bone, Lung, Lymph Node"
      
      index	<-	BRCA_met$met_loc	==	"Bone,Liver,Bone,Liver"
      BRCA_met$met_loc[index]	<-	"Bone, Liver"
      
      index	<-	BRCA_met$met_loc	==	"Bone,Liver,Bone|Liver"
      BRCA_met$met_loc[index]	<-	"Bone, Liver"
      
      index	<-	BRCA_met$met_loc	==	"Liver,"
      BRCA_met$met_loc[index]	<-	"Liver"
      
      index	<-	BRCA_met$met_loc	==	"Lung,,Lung"
      BRCA_met$met_loc[index]	<-	"Lung"
      
      index	<-	BRCA_met$met_loc	==	"Lung,Bone,Liver,Other, specify,Lung"
      BRCA_met$met_loc[index]	<-	"Lung, Liver, Bone"
      
      index	<-	BRCA_met$met_loc	==	"Lung,Bone,Lung"
      BRCA_met$met_loc[index]	<-	"Lung,Bone"
      
      index	<-	BRCA_met$met_loc	==	"Lung,Other, specify"
      BRCA_met$met_loc[index]	<-	"Lung"
      
      index	<-	BRCA_met$met_loc	==	"Other, specify,"
      BRCA_met$met_loc[index]	<-	NA
      
      index	<-	BRCA_met$met_loc	==	"Other, specify,,Brain"
      BRCA_met$met_loc[index]	<-	"Brain"
      
      index	<-	BRCA_met$met_loc	==	"Other, specify,,Other, specify,lung, bone, liver"
      BRCA_met$met_loc[index]	<-	"Lung, Bone, Liver"
      
      
      
      
      write.csv(BRCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
      
      print("BRCA Done")
      rm(BRCA_met)
    }
  
  if(i== "TCGA-COAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    COAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,new_neoplasm_event_type, other_malignancy_anatomic_site, 
                                              site_of_additional_surgery_new_tumor_event_mets,number_of_lymphnodes_positive_by_he)) %>%
                                tidyr::unite(met_loc,new_neoplasm_event_type:site_of_additional_surgery_new_tumor_event_mets, na.rm =TRUE,sep = ",") %>%
                                mutate_each(funs(empty_as_na)) %>%
                                mutate(COAD_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                mutate(COAD_met, Metastatic_status = ifelse(met_loc ==",," & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 
                                
      
                            
                              
    index <- is.na(COAD_met$LymphNodeStatus)
    COAD_met$LymphNodeStatus[index] <- 0
    
    
    index <- COAD_met$malignancy_type == "Prior Malignancy"
    COAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(COAD_met$Metastatic_status)
    COAD_met$Metastatic_status[index] <- 0
    
    
    index <- COAD_met$met_loc == ",,"
    COAD_met$met_loc[index] <- NA
    
    index <- is.na(COAD_met$number_of_lymphnodes_positive_by_he)
    COAD_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    index <- is.na(COAD_met$met_loc)
    COAD_met$Metastatic_status[index] <- 0
    
    index	<-	COAD_met$met_loc	==	",,Liver"
    COAD_met$met_loc[index]	<-	"Liver"
    
    index	<-	COAD_met$met_loc	==	",,Other"
    COAD_met$met_loc[index]	<-	NA
    
    index	<-	COAD_met$met_loc	==	",Ascending colon|Lymphoma,"
    COAD_met$met_loc[index]	<-	"Colon, Lymphoma"
    
    index	<-	COAD_met$met_loc	==	",Back,"
    COAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	COAD_met$met_loc	==	",Back|Breast,"
    COAD_met$met_loc[index]	<-	"Breast, Skin"
    
    index	<-	COAD_met$met_loc	==	",Bladder,"
    COAD_met$met_loc[index]	<-	"Bladder"
    
    index	<-	COAD_met$met_loc	==	",Blood,"
    COAD_met$met_loc[index]	<-	"Blood"
    
    index	<-	COAD_met$met_loc	==	",Breast,"
    COAD_met$met_loc[index]	<-	"Breast"
    
    index	<-	COAD_met$met_loc	==	",Breast|Other,"
    COAD_met$met_loc[index]	<-	"Breast"
    
    index	<-	COAD_met$met_loc	==	",Bronchus and Lung,"
    COAD_met$met_loc[index]	<-	"Lung"
    
    index	<-	COAD_met$met_loc	==	",Cecum,"
    COAD_met$met_loc[index]	<-	"Cecum"
    
    index	<-	COAD_met$met_loc	==	",Colon,"
    COAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	COAD_met$met_loc	==	",Colon|Bladder|Prostate,"
    COAD_met$met_loc[index]	<-	"Colon, Bladder, Prostate"
    
    index	<-	COAD_met$met_loc	==	",Ear,"
    COAD_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	COAD_met$met_loc	==	",Gum,"
    COAD_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	COAD_met$met_loc	==	",Kidney,"
    COAD_met$met_loc[index]	<-	"Kidney"
    
    index	<-	COAD_met$met_loc	==	",Kidney|Bladder,"
    COAD_met$met_loc[index]	<-	"Kidney, Bladder"
    
    index	<-	COAD_met$met_loc	==	",Larynx,"
    COAD_met$met_loc[index]	<-	"Larynx"
    
    index	<-	COAD_met$met_loc	==	",Lung,"
    COAD_met$met_loc[index]	<-	"Lung"
    
    index	<-	COAD_met$met_loc	==	",Lymphoma,"
    COAD_met$met_loc[index]	<-	"Lymphoma"
    
    index	<-	COAD_met$met_loc	==	",Other,"
    COAD_met$met_loc[index]	<-	NA
    
    index	<-	COAD_met$met_loc	==	",Other|Rectum,"
    COAD_met$met_loc[index]	<-	"Rectum"
    
    index	<-	COAD_met$met_loc	==	",Parotid gland,"
    COAD_met$met_loc[index]	<-	"Parotid"
    
    index	<-	COAD_met$met_loc	==	",Prostate,"
    COAD_met$met_loc[index]	<-	"Prostate"
    
    index	<-	COAD_met$met_loc	==	",Prostate|Bone marrow,"
    COAD_met$met_loc[index]	<-	"Prostate, Bone"
    
    index	<-	COAD_met$met_loc	==	",Prostate|Cecum,"
    COAD_met$met_loc[index]	<-	"Prostate, Cecum"
    
    index	<-	COAD_met$met_loc	==	",Rectum,"
    COAD_met$met_loc[index]	<-	"Rectum"
    
    index	<-	COAD_met$met_loc	==	",Scalp,"
    COAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	COAD_met$met_loc	==	",Sigmoid colon,"
    COAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	COAD_met$met_loc	==	",Sigmoid colon|Other,"
    COAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	COAD_met$met_loc	==	",Stomach,"
    COAD_met$met_loc[index]	<-	"Stomach"
    
    index	<-	COAD_met$met_loc	==	",Testicle,"
    COAD_met$met_loc[index]	<-	"Testicle"
    
    index	<-	COAD_met$met_loc	==	",Thyroid gland,"
    COAD_met$met_loc[index]	<-	"Thyroid"
    
    index	<-	COAD_met$met_loc	==	",Uterus,"
    COAD_met$met_loc[index]	<-	"Uterus"
    
    index	<-	COAD_met$met_loc	==	"Locoregional Disease,,"
    COAD_met$met_loc[index]	<-	NA
    
    index	<-	COAD_met$met_loc	==	"Locoregional Disease|Metastatic,,Liver"
    COAD_met$met_loc[index]	<-	"Liver"
    
    index	<-	COAD_met$met_loc	==	"Metastatic,,Liver"
    COAD_met$met_loc[index]	<-	"Liver"
    
    index	<-	COAD_met$met_loc	==	"Metastatic,,Liver|Lung"
    COAD_met$met_loc[index]	<-	"Liver, Lung"
    
    index	<-	COAD_met$met_loc	==	"Metastatic,,Other"
    COAD_met$met_loc[index]	<-	NA
    COAD_met$Metastatic_status[index] <- 1
    
    
  
    
    write.csv(COAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("COAD Done")
    rm(COAD_met)
  }
  
  if(i== "TCGA-ESCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    ESCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive_by_he,new_neoplasm_occurrence_anatomic_site_text,other_malignancy_anatomic_site)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                tidyr::unite(met_loc,new_neoplasm_occurrence_anatomic_site_text:other_malignancy_anatomic_site, na.rm =TRUE,sep = ",") %>%
                                mutate(ESCA_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(ESCA_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 
      
    index <- is.na(ESCA_met$LymphNodeStatus)
    ESCA_met$LymphNodeStatus[index] <- 0
    
    
    index <- ESCA_met$malignancy_type == "Prior Malignancy"
    ESCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(ESCA_met$Metastatic_status)
    ESCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(ESCA_met$number_of_lymphnodes_positive_by_he)
    ESCA_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    index <- is.na(ESCA_met$met_loc)
    ESCA_met$Metastatic_status[index] <- 0
    
    index <- ESCA_met$met_loc == "Adrenal" 
    ESCA_met$met_loc[index] <- "Adrenal Gland"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Adrenal Gland|Colon|Thigh" 
    ESCA_met$met_loc[index] <- "Adrenal Gland,Colon,Thigh"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "bone" 
    ESCA_met$met_loc[index] <- "Bone"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "brain" 
    ESCA_met$met_loc[index] <- "Brain"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Brain - right occipital lobe and cerebellum" 
    ESCA_met$met_loc[index] <- "Brain"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Ear" 
    ESCA_met$met_loc[index] <- "Head and Neck"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Head and neck lymph node" 
    ESCA_met$met_loc[index] <- "Head and Neck,Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Neck" 
    ESCA_met$met_loc[index] <- "Head and Neck"
    ESCA_met$Metastatic_status[index] <- 1 
    
    index <- ESCA_met$met_loc == "Pleural Mets" 
    ESCA_met$met_loc[index] <- "Lung"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "regional lymph nodes" 
    ESCA_met$met_loc[index] <- "Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Right diaphragm-HG spindle cell sarcoma" 
    ESCA_met$met_loc[index] <- "Lung"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "stomach" 
    ESCA_met$met_loc[index] <- "Stomach"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "left supraclavicular nodal metastasis" 
    ESCA_met$met_loc[index] <- "Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "left gastric node" 
    ESCA_met$met_loc[index] <- "Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "lymph node|Lymph Node" 
    ESCA_met$met_loc[index] <- "Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "peritoneum" 
    ESCA_met$met_loc[index] <- "Peritoneum"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$met_loc == "Mediastinal lymph node" 
    ESCA_met$met_loc[index] <- "Lymph Node"
    ESCA_met$Metastatic_status[index] <- 1
    
    index <- ESCA_met$malignancy_type == "Prior Malignancy"
    ESCA_met$Metastatic_status[index] <- 0
  
    
    
    write.csv(ESCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("ESCA Done")
    rm(ESCA_met)
  } 
  
  if(i== "TCGA-HNSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    HNSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive_by_he,other_malignancy_anatomic_site, other_malignancy_anatomic_site_text,
                                      new_neoplasm_event_occurrence_anatomic_site,new_neoplasm_occurrence_anatomic_site_text)) %>%
                                      mutate_each(funs(empty_as_na))%>%
                                      tidyr::unite(met_loc,other_malignancy_anatomic_site:new_neoplasm_occurrence_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
                                      mutate(HNSC_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive_by_he) & number_of_lymphnodes_positive_by_he >0, 1,0)) %>%
                                      mutate_each(funs(empty_as_na))%>%
                                      mutate(HNSC_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0 & malignancy_type != "Prior Malignancy", 0, 1)) 
    
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, "Distant Metastasis|", "")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, "Distant Metastasis,", "")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, ",lung", "Lung")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, ",Lung", "Lung")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, "LUNG", "Lung")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, "Lungs", "Lung")
    HNSC_met$met_loc <- stringr::str_replace_all(HNSC_met$met_loc, "Lung (distant mets or possible new primary)", "Lung")
    
    
    
    
    index	<-	HNSC_met$met_loc	==	",Brain Metastasis"
    HNSC_met$met_loc[index]	<-	"Brain"
    index	<-	HNSC_met$met_loc	==	",Brain|Cerebellum"
    HNSC_met$met_loc[index]	<-	"Brain"
    index	<-	HNSC_met$met_loc	==	",esophagus|Esophagus"
    HNSC_met$met_loc[index]	<-	"Esophagus"
    index	<-	HNSC_met$met_loc	==	",L neck subdermal mets"
    HNSC_met$met_loc[index]	<-	"Skin"
    index	<-	HNSC_met$met_loc	==	",left axilla"
    HNSC_met$met_loc[index]	<-	"Brain"
    index	<-	HNSC_met$met_loc	==	",suprasternal"
    HNSC_met$met_loc[index]	<-	"Suprasternal"
    index	<-	HNSC_met$met_loc	==	",trachea|Lung, Right Lower Lobe"
    HNSC_met$met_loc[index]	<-	"Lung"
    index	<-	HNSC_met$met_loc	==	",tracheostoma"
    HNSC_met$met_loc[index]	<-	"Lung"
    index	<-	HNSC_met$met_loc	==	"|Cervical Lymph NodesLung"
    HNSC_met$met_loc[index]	<-	"Lymph Node"
    index	<-	HNSC_met$met_loc	==	"|Cervical Lymph NodesLung|parapharyngeal space"
    HNSC_met$met_loc[index]	<-	"Lymph Node"
    index	<-	HNSC_met$met_loc	==	"Acute myeloid leukemia"
    HNSC_met$met_loc[index]	<-	"Bone"
    index	<-	HNSC_met$met_loc	==	"Oropharynx,OROPHARYNX, OVERLAPPING|oropharynx"
    HNSC_met$met_loc[index]	<-	"Oropharynx"
    index	<-	HNSC_met$met_loc	==	"Cervical Lymph Nodes,bilateral"
    HNSC_met$met_loc[index]	<-	"Lymph Node"
    index	<-	HNSC_met$met_loc	==	"Cervical Lymph Nodes"
    HNSC_met$met_loc[index]	<-	"Lymph Node"
    index	<-	HNSC_met$met_loc	==	"R. Lung"
    HNSC_met$met_loc[index]	<-	"Lung"
    index	<-	HNSC_met$met_loc	==	"Oropharynx,Oropharynx and cervical lymph node"
    HNSC_met$met_loc[index]	<-	"Oropharynx, Lymph Node"
    
    index	<-	HNSC_met$met_loc	==	"Prostate|Kidney"
    HNSC_met$met_loc[index]	<-	"Prostate, Kidney"
    
    index	<-	HNSC_met$met_loc	==	"Prostate|Hypopharynx"
    HNSC_met$met_loc[index]	<-	"Prostate, Hypopharynx"
    
    index	<-	HNSC_met$met_loc	==	"Other,skin of right cheek"
    HNSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	HNSC_met$met_loc	==	"Other,skin NOS"
    HNSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	HNSC_met$met_loc	==	"Other,skin"
    HNSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	HNSC_met$met_loc	==	"Oropharynx|Oral Cavity,Submandibular gland|tongue and floor of mouth"
    HNSC_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	HNSC_met$met_loc	==	"Lung, mediastinal lymph node|lung, mediastinal and hilar lymph nodes|liver, spleen, kidney, adrenal, bone|brain"
    HNSC_met$met_loc[index]	<-	"Lung,Liver,Brain,Bone,Adrenal Gland,Lymph Node"
    
    index	<-	HNSC_met$met_loc	==	"Oropharynx,L tonsil"
    HNSC_met$met_loc[index]	<-	"Oropharynx, Bone"
    
    index	<-	HNSC_met$met_loc	==	"mediastinum"
    HNSC_met$met_loc[index]	<-	"Mediastinum"
    
    
    index	<-	HNSC_met$met_loc	==	"LUL Lung"
    HNSC_met$met_loc[index]	<-	"Lung"
    
    index	<-	HNSC_met$met_loc	==	"Cervical Lymph Nodes,L NECK"
    HNSC_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	HNSC_met$met_loc	==	"Lung (distant mets or possible new primary)"
    HNSC_met$met_loc[index]	<-	"Lung"
    
    index	<-	HNSC_met$met_loc	==	"Oropharynx,tonsil"
    HNSC_met$met_loc[index]	<-	"Oropharynx,Bone"
    
    index <- HNSC_met$met_loc == "Distant Metastasis|" 
    HNSC_met$met_loc[index] <- ""
    HNSC_met$Metastatic_status[index] <- 1
    
    index <- HNSC_met$met_loc == "Distant Metastasis," 
    HNSC_met$met_loc[index] <- ""
    HNSC_met$Metastatic_status[index] <- 1
    
    
    index <- is.na(HNSC_met$LymphNodeStatus)
    HNSC_met$LymphNodeStatus[index] <- 0
    
    
    index <- HNSC_met$malignancy_type == "Prior Malignancy"
    HNSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(HNSC_met$Metastatic_status)
    HNSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(HNSC_met$number_of_lymphnodes_positive_by_he)
    HNSC_met$number_of_lymphnodes_positive_by_he[index] <- 0
    
    index <- HNSC_met$met_loc == "Lymph Node"
    HNSC_met$Metastatic_status[index] <- 0
    HNSC_met$LymphNodeStatus[index] <- 1
    
    write.csv(HNSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus.csv"))
    
    print("HNSC Done")
    rm(HNSC_met)
  } 
  
  if(i== "TCGA-KICH"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    KICH_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,number_of_lymphnodes_positive, 
                                              other_malignancy_anatomic_site)) %>%
                                mutate(KICH_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(KICH_met, met_loc = ifelse(malignancy_type == "Synchronous Malignancy", malignancy_type, NA)) %>%
                                mutate(KICH_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus ==0, 0, 1)) 
    
    write.csv(KICH_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KICH is bulshit but it's Done")
    rm(KICH_met)
  } 
  
  if(i== "TCGA-KIRC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive,malignancy_type,
                                              other_malignancy_anatomic_site, other_malignancy_anatomic_site_text)) %>% 
                              mutate_each(funs(empty_as_na))%>%
                              tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
                              mutate(KIRC_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
                              mutate(KIRC_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus >0 | is.na(LymphNodeStatus) | malignancy_type =="Prior Malignnacy" , 0, 1))
    
    index <- is.na(KIRC_met$LymphNodeStatus)
    KIRC_met$LymphNodeStatus[index] <- 0
    
    index <- KIRC_met$LymphNodeStatus >0
    KIRC_met$Metastatic_status[index] <- 1
    
    index <- KIRC_met$malignancy_type == "Prior Malignancy"
    KIRC_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRC_met$Metastatic_status)
    KIRC_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRC_met$number_of_lymphnodes)
    KIRC_met$number_of_lymphnodes_positive[index] <- 0
    
    
    index	<-	KIRC_met$met_loc	==	"Anus|Kidney"
    KIRC_met$met_loc[index]	<-	"Anus, Kidney"
    
    index	<-	KIRC_met$met_loc	==	"Back|Uterus"
    KIRC_met$met_loc[index]	<-	"Skin, Uterus"
    
    index	<-	KIRC_met$met_loc	==	"Bladder|Prostate"
    KIRC_met$met_loc[index]	<-	"Bladder, Prostate"
    
    index	<-	KIRC_met$met_loc	==	"Bladder|Uterus"
    KIRC_met$met_loc[index]	<-	"Bladder, Uterus"
    
    index	<-	KIRC_met$met_loc	==	"Bone marrow"
    KIRC_met$met_loc[index]	<-	"Bone"
    
    index	<-	KIRC_met$met_loc	==	"Colon|Prostate"
    KIRC_met$met_loc[index]	<-	"Colon,Prostate"
    
    index	<-	KIRC_met$met_loc	==	"Ear"
    KIRC_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	KIRC_met$met_loc	==	"Forehead"
    KIRC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	KIRC_met$met_loc	==	"Head/Neck"
    KIRC_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	KIRC_met$met_loc	==	"Kidney|Eye"
    KIRC_met$met_loc[index]	<-	"Kidney, Eye"
    
    index	<-	KIRC_met$met_loc	==	"Kidney|Pancreas"
    KIRC_met$met_loc[index]	<-	"Kidney, Pancreas"
    
    index	<-	KIRC_met$met_loc	==	"Kidney|Prostate"
    KIRC_met$met_loc[index]	<-	"Kidney, Prostate"
    
    index	<-	KIRC_met$met_loc	==	"Lymph node(s)"
    KIRC_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	KIRC_met$met_loc	==	"Neck"
    KIRC_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	KIRC_met$met_loc	==	"Other"
    KIRC_met$met_loc[index]	<-	NA
    
    index	<-	KIRC_met$met_loc	==	"other"
    KIRC_met$met_loc[index]	<-	NA
    
    index	<-	KIRC_met$met_loc	==	"Other,site not stated"
    KIRC_met$met_loc[index]	<-	NA
    
    index	<-	KIRC_met$met_loc	==	"Other,Skin"
    KIRC_met$met_loc[index]	<-	"Skin"
    
    index	<-	KIRC_met$met_loc	==	"Pancreas|Kidney"
    KIRC_met$met_loc[index]	<-	"Pancreas, Kidney"
    
    index	<-	KIRC_met$met_loc	==	"Prostate|Kidney"
    KIRC_met$met_loc[index]	<-	"Prostate, Kidney"
    
    index	<-	KIRC_met$met_loc	==	"Spine"
    KIRC_met$met_loc[index]	<-	"Bone"
    
    index	<-	KIRC_met$met_loc	==	"Thyroid gland|Kidney"
    KIRC_met$met_loc[index]	<-	"Thyroid gland, Kidney"
    
    index	<-	KIRC_met$met_loc	==	"Other,Skin Cancer"
    KIRC_met$met_loc[index]	<-	"Skin"
    
    
    write.csv(KIRC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRC Done")
    rm(KIRC_met)
  } 
  
  if(i== "TCGA-KIRP"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    KIRP_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, number_of_lymphnodes_positive,malignancy_type,
                                              other_malignancy_anatomic_site, other_malignancy_anatomic_site_text)) %>% 
      mutate_each(funs(empty_as_na))%>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate(KIRP_met, LymphNodeStatus = ifelse(!is.na(number_of_lymphnodes_positive) & number_of_lymphnodes_positive >0, 1,0)) %>%
      mutate(KIRP_met, Metastatic_status = ifelse(is.na(met_loc) & LymphNodeStatus >0 | is.na(LymphNodeStatus) | malignancy_type =="Prior Malignnacy" , 0, 1))
    
    index <- is.na(KIRP_met$LymphNodeStatus)
    KIRP_met$LymphNodeStatus[index] <- 0
    
    index <- KIRP_met$LymphNodeStatus >0
    KIRP_met$Metastatic_status[index] <- 1
    
    index <- KIRP_met$malignancy_type == "Prior Malignancy"
    KIRP_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRP_met$Metastatic_status)
    KIRP_met$Metastatic_status[index] <- 0
    
    index <- is.na(KIRP_met$number_of_lymphnodes)
    KIRP_met$number_of_lymphnodes_positive[index] <- 0
    
    
    index	<-	KIRP_met$met_loc	==	"Back"
    KIRP_met$met_loc[index]	<-	"Skin"
    
    index	<-	KIRP_met$met_loc	==	"Back|Eye"
    KIRP_met$met_loc[index]	<-	"Back,Eye"
    
    index	<-	KIRP_met$met_loc	==	"Breast|Ovary"
    KIRP_met$met_loc[index]	<-	"Breast,Ovary"
    
    index	<-	KIRP_met$met_loc	==	"Colon and rectum|Colon"
    KIRP_met$met_loc[index]	<-	"Colon, Rectum"
    
    index	<-	KIRP_met$met_loc	==	"Ear|Prostate"
    KIRP_met$met_loc[index]	<-	"Head and Neck, Prostate"
    
    index	<-	KIRP_met$met_loc	==	"Kidney|Prostate"
    KIRP_met$met_loc[index]	<-	"Kidney, Prostate"
    
    index	<-	KIRP_met$met_loc	==	"Lymph node(s)"
    KIRP_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	KIRP_met$met_loc	==	"Other,Face"
    KIRP_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	KIRP_met$met_loc	==	"Other,Skin, nothing more specific provided"
    KIRP_met$met_loc[index]	<-	"Skin"
    
    index	<-	KIRP_met$met_loc	==	"Other|Kidney,Skin"
    KIRP_met$met_loc[index]	<-	"Kidney, Skin"
    
    index	<-	KIRP_met$met_loc	==	"Thigh"
    KIRP_met$met_loc[index]	<-	"Skin"
    
    index	<-	KIRP_met$met_loc	==	"Thyroid gland|Prostate"
    KIRP_met$met_loc[index]	<-	"Thyroid gland, Prostate"
    
    index	<-	KIRP_met$met_loc	==	"Transverse Colon"
    KIRP_met$met_loc[index]	<-	"Colon"
    
    
    
    write.csv(KIRP_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("KIRP Done")
    rm(KIRP_met)
  
  
  
  } 
  
  if(i== "TCGA-LIHC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LIHC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode, malignancy_type, other_malignancy_anatomic_site, other_malignancy_anatomic_site_text,new_neoplasm_event_occurrence_anatomic_site)) %>% 
      tidyr::unite(met_loc,other_malignancy_anatomic_site:new_neoplasm_event_occurrence_anatomic_site, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(LIHC_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" , 1, 0))
    
    LIHC_met[LIHC_met == ",,"] <- NA
    
    
    index <- !is.na(LIHC_met$met_loc) & is.na(LIHC_met$Metastatic_status)
    LIHC_met$Metastatic_status[index] <- 1
    
    index <- LIHC_met$met_loc
    LIHC_met$Metastatic_status[index] <- 1
    
    index <- LIHC_met$malignancy_type == "Prior Malignancy"
    LIHC_met$Metastatic_status[index] <- 0
    
    index <- is.na(LIHC_met$Metastatic_status)
    LIHC_met$Metastatic_status[index] <- 0
    
    LIHC_met$met_loc <- str_replace_all(LIHC_met$met_loc, ",,","")
    
    index	<-	LIHC_met$met_loc	==	"Arm"
    LIHC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LIHC_met$met_loc	==	"Back"
    LIHC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LIHC_met$met_loc	==	"Back|Breast"
    LIHC_met$met_loc[index]	<-	"Skin, Breast"
    
    index	<-	LIHC_met$met_loc	==	"Bladder|Prostate"
    LIHC_met$met_loc[index]	<-	"Bladder, Prostate"
    
    index	<-	LIHC_met$met_loc	==	"Bone|Liver"
    LIHC_met$met_loc[index]	<-	"Bone, Liver"
    
    index	<-	LIHC_met$met_loc	==	"Brain|Liver"
    LIHC_met$met_loc[index]	<-	"Brain, Liver"
    
    index	<-	LIHC_met$met_loc	==	"Cecum"
    LIHC_met$met_loc[index]	<-	"Intestine"
    
    index	<-	LIHC_met$met_loc	==	"Colon and rectum"
    LIHC_met$met_loc[index]	<-	"Colon, Rectum"
    
    index	<-	LIHC_met$met_loc	==	"Ear"
    LIHC_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	LIHC_met$met_loc	==	"Forehead"
    LIHC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LIHC_met$met_loc	==	"Head - face or neck, NOS"
    LIHC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LIHC_met$met_loc	==	"Leg"
    LIHC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Bone"
    LIHC_met$met_loc[index]	<-	"Liver, Bone"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Brain"
    LIHC_met$met_loc[index]	<-	"Liver, Brain"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Lung"
    LIHC_met$met_loc[index]	<-	"Liver, Lung"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Lung|Other, specify"
    LIHC_met$met_loc[index]	<-	"Liver, Lung"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Other, specify"
    LIHC_met$met_loc[index]	<-	"Liver"
    
    index	<-	LIHC_met$met_loc	==	"Liver|Other, specify|Lung"
    LIHC_met$met_loc[index]	<-	"Liver, Lung"
    
    index	<-	LIHC_met$met_loc	==	"Lung|Bone"
    LIHC_met$met_loc[index]	<-	"Lung, Bone"
    
    index	<-	LIHC_met$met_loc	==	"Lung|Brain"
    LIHC_met$met_loc[index]	<-	"Lung, Brain"
    
    index	<-	LIHC_met$met_loc	==	"Lymph Node(s)"
    LIHC_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	LIHC_met$met_loc	==	"Other, specify"
    LIHC_met$met_loc[index]	<-	NA
    
    index	<-	LIHC_met$met_loc	==	"Other, specify|Bone"
    LIHC_met$met_loc[index]	<-	"Bone"
    
    index	<-	LIHC_met$met_loc	==	"Other,Endometrium,"
    LIHC_met$met_loc[index]	<-	"Uterus"
    
    index	<-	LIHC_met$met_loc	==	"Other,SKIN,"
    LIHC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LIHC_met$met_loc	==	"Other|Scalp,Cheek,"
    LIHC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LIHC_met$met_loc	==	"Prostate|Laryngopharynx"
    LIHC_met$met_loc[index]	<-	"Prostate, Laryngopharynx"
    
    index	<-	LIHC_met$met_loc	==	"Pulmonary"
    LIHC_met$met_loc[index]	<-	"Lung"
    
    index	<-	LIHC_met$met_loc	==	"Sigmoid colon|Prostate"
    LIHC_met$met_loc[index]	<-	"Colon, Prosate"
    
    
    View(table(LIHC_met$met_loc))
    
    
    write.csv(LIHC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LIHC Done")
    rm(LIHC_met)
    
    
  } 
  
  
  if(i== "TCGA-LUAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 

    LUAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,
                                              other_malignancy_anatomic_site, pathologic_N,pathologic_M))%>% 
      tidyr::unite(met_loc,other_malignancy_anatomic_site, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(LUAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(LUAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    LUAD_met[LUAD_met == ",,"] <- NA
    
    LUAD_met[LUAD_met == "MX"] <- "M0"
    
    index <- is.na(LUAD_met$pathologic_M)
    LUAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(LUAD_met$pathologic_N)
    LUAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(LUAD_met$met_loc) & is.na(LUAD_met$Metastatic_status)
    LUAD_met$Metastatic_status[index] <- 0
    
    index <- LUAD_met$malignancy_type == "Prior Malignancy"
    LUAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(LUAD_met$Metastatic_status)
    LUAD_met$Metastatic_status[index] <- 0
    
    LUAD_met$met_loc <- str_replace_all(LUAD_met$met_loc, ",,","")
    
    
    View(table(LUAD_met$met_loc))
    
    
    index	<-	LUAD_met$met_loc	==	"Arm"
    LUAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUAD_met$met_loc	==	"Arm|Lung"
    LUAD_met$met_loc[index]	<-	"Arm, Lung"
    
    index	<-	LUAD_met$met_loc	==	"Breast|Uterus|Cervix"
    LUAD_met$met_loc[index]	<-	"Breast, Uterus, Cervix"
    
    index	<-	LUAD_met$met_loc	==	"Chest wall"
    LUAD_met$met_loc[index]	<-	"Bone"
    
    index	<-	LUAD_met$met_loc	==	"Descending colon"
    LUAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	LUAD_met$met_loc	==	"Ear|Other"
    LUAD_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	LUAD_met$met_loc	==	"Forehead"
    LUAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUAD_met$met_loc	==	"Hand"
    LUAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUAD_met$met_loc	==	"Head - face or neck, NOS"
    LUAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUAD_met$met_loc	==	"Kidney|Bowel"
    LUAD_met$met_loc[index]	<-	"Kidney, Intestine"
    
    index	<-	LUAD_met$met_loc	==	"Lymph Node(s) Cervical"
    LUAD_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	LUAD_met$met_loc	==	"Mouth"
    LUAD_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	LUAD_met$met_loc	==	"Neck"
    LUAD_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	LUAD_met$met_loc	==	"Other"
    LUAD_met$met_loc[index]	<-	NA
    
    index	<-	LUAD_met$met_loc	==	"Other|Blood"
    LUAD_met$met_loc[index]	<-	"Blood"
    
    index	<-	LUAD_met$met_loc	==	"Rectum|Bladder"
    LUAD_met$met_loc[index]	<-	"Rectum, Bladder"
    
    index	<-	LUAD_met$met_loc	==	"Rectum|Prostate"
    LUAD_met$met_loc[index]	<-	"Rectum, Prostate"
    
    index	<-	LUAD_met$met_loc	==	"Spleen|Breast"
    LUAD_met$met_loc[index]	<-	"Spleen, Breast"
    
    index	<-	LUAD_met$met_loc	==	"Throat|Jaw"
    LUAD_met$met_loc[index]	<-	"Esophagus, Bone"
    
    
  
    write.csv(LUAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUAD Done")
    rm(LUAD_met)
    
    
    
  } 
  
  if(i== "TCGA-LUSC"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    LUSC_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
                                tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
                                mutate_each(funs(empty_as_na))%>%
                                mutate(LUSC_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                                            | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
                                mutate(LUSC_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    LUSC_met[LUSC_met == ","] <- NA
    
    LUSC_met[LUSC_met == "MX"] <- "M0"
    
    index <- is.na(LUSC_met$pathologic_M)
    LUSC_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(LUSC_met$pathologic_N)
    LUSC_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(LUSC_met$met_loc) & is.na(LUSC_met$Metastatic_status)
    LUSC_met$Metastatic_status[index] <- 0
    
    index <- LUSC_met$malignancy_type == "Prior Malignancy"
    LUSC_met$Metastatic_status[index] <- 0
    
    index <- is.na(LUSC_met$Metastatic_status)
    LUSC_met$Metastatic_status[index] <- 0
    
    LUSC_met$met_loc <- str_replace_all(LUSC_met$met_loc, ",,","")
    
    index	<-	LUSC_met$met_loc	==	"Bladder,"
    LUSC_met$met_loc[index]	<-	"Bladder"
    
    index	<-	LUSC_met$met_loc	==	"Blood,"
    LUSC_met$met_loc[index]	<-	"Blood"
    
    index	<-	LUSC_met$met_loc	==	"Breast,"
    LUSC_met$met_loc[index]	<-	"Breast"
    
    index	<-	LUSC_met$met_loc	==	"Breast|Thymus,"
    LUSC_met$met_loc[index]	<-	"Breast, Thymus"
    
    index	<-	LUSC_met$met_loc	==	"Cervix,"
    LUSC_met$met_loc[index]	<-	"Cervix"
    
    index	<-	LUSC_met$met_loc	==	"Colon,"
    LUSC_met$met_loc[index]	<-	"Colon"
    
    index	<-	LUSC_met$met_loc	==	"Eye,"
    LUSC_met$met_loc[index]	<-	"Eye"
    
    index	<-	LUSC_met$met_loc	==	"Foot,"
    LUSC_met$met_loc[index]	<-	"Foot, Bone"
    
    index	<-	LUSC_met$met_loc	==	"Forehead|Lung,"
    LUSC_met$met_loc[index]	<-	"Head and Neck, Skin, Lung"
    
    index	<-	LUSC_met$met_loc	==	"Forehead|Other,Face"
    LUSC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUSC_met$met_loc	==	"Head - face or neck, NOS,"
    LUSC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUSC_met$met_loc	==	"Kidney,"
    LUSC_met$met_loc[index]	<-	"Kidney"
    
    index	<-	LUSC_met$met_loc	==	"Larynx,"
    LUSC_met$met_loc[index]	<-	"Larynx"
    
    index	<-	LUSC_met$met_loc	==	"Lip|Other|Shoulder,NOSE"
    LUSC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUSC_met$met_loc	==	"Liver,"
    LUSC_met$met_loc[index]	<-	"Liver"
    
    index	<-	LUSC_met$met_loc	==	"Lung,"
    LUSC_met$met_loc[index]	<-	"Lung"
    
    index	<-	LUSC_met$met_loc	==	"Lymphoma,"
    LUSC_met$met_loc[index]	<-	"Lymphoma"
    
    index	<-	LUSC_met$met_loc	==	"Nasopharynx,"
    LUSC_met$met_loc[index]	<-	"Nasopharynx"
    
    index	<-	LUSC_met$met_loc	==	"Neck"
    LUSC_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	LUSC_met$met_loc	==	"Other,Basel Cell Skin Cancer"
    LUSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUSC_met$met_loc	==	"Other,gynecological (not otherwise specified)"
    LUSC_met$met_loc[index]	<-	"Cervix"
    
    index	<-	LUSC_met$met_loc	==	"Other,Right Tongue"
    LUSC_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	LUSC_met$met_loc	==	"Other,skin"
    LUSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUSC_met$met_loc	==	"Other,Skin"
    LUSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUSC_met$met_loc	==	"Other,Skin- Nose - Nasal tip"
    LUSC_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	LUSC_met$met_loc	==	"Other,Skin, NOS"
    LUSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUSC_met$met_loc	==	"Other|Arm,Left Chest"
    LUSC_met$met_loc[index]	<-	"Skin"
    
    index	<-	LUSC_met$met_loc	==	"Pelvis|Bladder,"
    LUSC_met$met_loc[index]	<-	"Pelvis, Bladder"
    
    index	<-	LUSC_met$met_loc	==	"Pharynx,"
    LUSC_met$met_loc[index]	<-	"Pharynx"
    
    index	<-	LUSC_met$met_loc	==	"Prostate,"
    LUSC_met$met_loc[index]	<-	"Prostate"
    
    index	<-	LUSC_met$met_loc	==	"Prostate|Back,"
    LUSC_met$met_loc[index]	<-	"Prostate, Skin"
    
    index	<-	LUSC_met$met_loc	==	"Prostate|Lung,"
    LUSC_met$met_loc[index]	<-	"Prostate, Lung"
    
    index	<-	LUSC_met$met_loc	==	"Thymus,"
    LUSC_met$met_loc[index]	<-	"Thymus"
    
    index	<-	LUSC_met$met_loc	==	"Tongue,"
    LUSC_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	LUSC_met$met_loc	==	"Tonsil,"
    LUSC_met$met_loc[index]	<-	"Oral Cavity, Bone"
    
    index	<-	LUSC_met$met_loc	==	"Uterus,"
    LUSC_met$met_loc[index]	<-	"Uterus"
    
    
    
    
    write.csv(LUSC_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("LUSC Done")
    rm(LUSC_met)
    
    
    
  } 
  
  if(i== "TCGA-PRAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    PRAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(PRAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(PRAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    PRAD_met[PRAD_met == ","] <- NA
    
    PRAD_met[PRAD_met == "MX"] <- "M0"
    
    index <- is.na(PRAD_met$pathologic_M)
    PRAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(PRAD_met$pathologic_N)
    PRAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(PRAD_met$met_loc) & is.na(PRAD_met$Metastatic_status)
    PRAD_met$Metastatic_status[index] <- 0
    
    index <- PRAD_met$malignancy_type == "Prior Malignancy"
    PRAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(PRAD_met$Metastatic_status)
    PRAD_met$Metastatic_status[index] <- 0
    
    PRAD_met$met_loc <- str_replace_all(PRAD_met$met_loc, ",,","")
    
    
    index <- PRAD_met$LymphNodeStatus > 0
    PRAD_met$Metastatic_status[index] <- 1
    
    index	<-	PRAD_met$met_loc	==	"Arm,"
    PRAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	PRAD_met$met_loc	==	"Bladder,"
    PRAD_met$met_loc[index]	<-	"Bladder"
    
    index	<-	PRAD_met$met_loc	==	"Chest wall,"
    PRAD_met$met_loc[index]	<-	"Bone"
    
    index	<-	PRAD_met$met_loc	==	"Colon,"
    PRAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	PRAD_met$met_loc	==	"Esophagus,"
    PRAD_met$met_loc[index]	<-	"Esophagus"
    
    index	<-	PRAD_met$met_loc	==	"Femur,"
    PRAD_met$met_loc[index]	<-	"Bone"
    
    index	<-	PRAD_met$met_loc	==	"Head - face or neck, NOS,"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Kidney,"
    PRAD_met$met_loc[index]	<-	"Kidney"
    
    index	<-	PRAD_met$met_loc	==	"Lung,"
    PRAD_met$met_loc[index]	<-	"Lung"
    
    index	<-	PRAD_met$met_loc	==	"Lymph node(s),"
    PRAD_met$met_loc[index]	<-	"Lymph Node"
    
    index	<-	PRAD_met$met_loc	==	"Lymphoma,"
    PRAD_met$met_loc[index]	<-	"Lymphoma"
    
    index	<-	PRAD_met$met_loc	==	"Mammary Gland,"
    PRAD_met$met_loc[index]	<-	"Breast"
    
    index	<-	PRAD_met$met_loc	==	"Neck,"
    PRAD_met$met_loc[index]	<-	"Head and Neck"
    
    index	<-	PRAD_met$met_loc	==	"Other,Eye lid"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Other,Nose"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Other,skin"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Other,Skin"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Other,Skin of Arms and Back"
    PRAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	PRAD_met$met_loc	==	"Other,Skin right front parietal scalp"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Pancreas,"
    PRAD_met$met_loc[index]	<-	"Pancreas"
    PRAD_met$Metastatic_status[index] <- 1
    
    index	<-	PRAD_met$met_loc	==	"Scalp,"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Scalp|Other,Infraorbital region"
    PRAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	PRAD_met$met_loc	==	"Thyroid gland,"
    PRAD_met$met_loc[index]	<-	"Thyroid"
    
    index <- PRAD_met$pathologic_N =="N1"
    PRAD_met$Metastatic_status[index] <- 0
    
    
    write.csv(PRAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("PRAD Done")
    rm(PRAD_met)
    
  } 
  
  if(i== "TCGA-STAD"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    STAD_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(STAD_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(STAD_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    STAD_met[STAD_met == ","] <- NA
    
    STAD_met[STAD_met == "MX"] <- "M0"
    
    index <- is.na(STAD_met$pathologic_M)
    STAD_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(STAD_met$pathologic_N)
    STAD_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(STAD_met$met_loc) & is.na(STAD_met$Metastatic_status)
    STAD_met$Metastatic_status[index] <- 0
    
    index <- STAD_met$malignancy_type == "Prior Malignancy"
    STAD_met$Metastatic_status[index] <- 0
    
    index <- is.na(STAD_met$Metastatic_status)
    STAD_met$Metastatic_status[index] <- 0
    
    STAD_met$met_loc <- str_replace_all(STAD_met$met_loc, ",,","")
    
    
    index <- STAD_met$LymphNodeStatus > 0
    STAD_met$Metastatic_status[index] <- 1
    
    
    index	<-	STAD_met$met_loc	==	"Bladder,"
    STAD_met$met_loc[index]	<-	"Bladder"
    
    index	<-	STAD_met$met_loc	==	"Breast,"
    STAD_met$met_loc[index]	<-	"Breast"
    
    index	<-	STAD_met$met_loc	==	"Colon and rectum,"
    STAD_met$met_loc[index]	<-	"Colon, Rectum"
    
    index	<-	STAD_met$met_loc	==	"Kidney,"
    STAD_met$met_loc[index]	<-	"Kidney"
    
    index	<-	STAD_met$met_loc	==	"Kidney|Prostate,"
    STAD_met$met_loc[index]	<-	"Kidney, Prostate"
    
    index	<-	STAD_met$met_loc	==	"Other,Skin"
    STAD_met$met_loc[index]	<-	"Skin"
    
    index	<-	STAD_met$met_loc	==	"Other,Skin - face"
    STAD_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	STAD_met$met_loc	==	"Ovary,"
    STAD_met$met_loc[index]	<-	"Ovary"
    
    index	<-	STAD_met$met_loc	==	"Prostate,"
    STAD_met$met_loc[index]	<-	"Prostate"
    
    index	<-	STAD_met$met_loc	==	"Prostate|Testicle|Bladder,"
    STAD_met$met_loc[index]	<-	"Prostate, Testicle, Bladder"
    
    index	<-	STAD_met$met_loc	==	"Sigmoid colon,"
    STAD_met$met_loc[index]	<-	"Colon"
    
    index	<-	STAD_met$met_loc	==	"Testicle,"
    STAD_met$met_loc[index]	<-	"Testicle"
    
    index <- STAD_met$pathologic_N =="N1"
    STAD_met$Metastatic_status[index] <- 0
    
    

    write.csv(STAD_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("STAD Done")
    rm(STAD_met)
    
    
    
  } 
  
  if(i== "TCGA-THCA"){
    
    dat <- as.data.frame(data.table::fread(str_glue("~/CSBL_shared/clinical/TCGA_xml/{i}.csv"))) 
    
    
    THCA_met <- as.data.frame(dat %>%
                                dplyr::select(bcr_patient_barcode,malignancy_type,other_malignancy_anatomic_site,other_malignancy_anatomic_site_text,
                                              other_malignancy_anatomic_site,pathologic_N,pathologic_M)) %>%
      tidyr::unite(met_loc,other_malignancy_anatomic_site:other_malignancy_anatomic_site_text, na.rm =TRUE,sep = ",") %>%
      mutate_each(funs(empty_as_na))%>%
      mutate(THCA_met, Metastatic_status = ifelse(!is.na(met_loc) & malignancy_type !="Prior Malignnacy" 
                                                  | pathologic_M == "M1"| pathologic_M == "M1a" | pathologic_M == "M1b" , 1, 0)) %>%
      mutate(THCA_met, LymphNodeStatus = ifelse(is.na(pathologic_N) | pathologic_N == "NX" | pathologic_N == "N0" , 0,1)) 
    
    THCA_met[THCA_met == ","] <- NA
    
    THCA_met[THCA_met == "MX"] <- "M0"
    
    index <- is.na(THCA_met$pathologic_M)
    THCA_met$pathologic_M[index] <- "M0"
    
    
    index <- is.na(THCA_met$pathologic_N)
    THCA_met$pathologic_N[index] <- "N0"
    
    
    index <- !is.na(THCA_met$met_loc) & is.na(THCA_met$Metastatic_status)
    THCA_met$Metastatic_status[index] <- 0
    
    index <- THCA_met$malignancy_type == "Prior Malignancy"
    THCA_met$Metastatic_status[index] <- 0
    
    index <- is.na(THCA_met$Metastatic_status)
    THCA_met$Metastatic_status[index] <- 0
    
    THCA_met$met_loc <- str_replace_all(THCA_met$met_loc, ",,","")
    
    
    index <- THCA_met$LymphNodeStatus > 0
    THCA_met$Metastatic_status[index] <- 1
    
    index	<-	THCA_met$met_loc	==	"Back,"
    THCA_met$met_loc[index]	<-	"Skin"
    
    index	<-	THCA_met$met_loc	==	"Brain,"
    THCA_met$met_loc[index]	<-	"Brain"
    
    index	<-	THCA_met$met_loc	==	"Breast,"
    THCA_met$met_loc[index]	<-	"Breast"
    
    index	<-	THCA_met$met_loc	==	"Buccal mucosa,"
    THCA_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	THCA_met$met_loc	==	"Calf,"
    THCA_met$met_loc[index]	<-	"Muscle"
    
    index	<-	THCA_met$met_loc	==	"Colon,"
    THCA_met$met_loc[index]	<-	"Colon"
    
    index	<-	THCA_met$met_loc	==	"Head - face or neck, NOS,"
    THCA_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	THCA_met$met_loc	==	"Head/Neck/Chest/Abd/Pelvis,"
    THCA_met$met_loc[index]	<-	"Head and Neck, Skin, Bone,Abdomen,Pelvis"
    
    index	<-	THCA_met$met_loc	==	"Lung,"
    THCA_met$met_loc[index]	<-	"Lung"
    
    index	<-	THCA_met$met_loc	==	"Lymphoma,"
    THCA_met$met_loc[index]	<-	"Lymphoma"
    
    index	<-	THCA_met$met_loc	==	"Mandible,"
    THCA_met$met_loc[index]	<-	"Head and Neck, Bone"
    
    index	<-	THCA_met$met_loc	==	"Oral cavity,"
    THCA_met$met_loc[index]	<-	"Oral Cavity"
    
    index	<-	THCA_met$met_loc	==	"Other,Nose"
    THCA_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	THCA_met$met_loc	==	"Other,skin of face back and neck"
    THCA_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	THCA_met$met_loc	==	"Prostate,"
    THCA_met$met_loc[index]	<-	"Prostate"
    
    index	<-	THCA_met$met_loc	==	"Rectum,"
    THCA_met$met_loc[index]	<-	"Rectum"
    
    index	<-	THCA_met$met_loc	==	"Scalp,"
    THCA_met$met_loc[index]	<-	"Head and Neck, Skin"
    
    index	<-	THCA_met$met_loc	==	"Thyroid gland,"
    THCA_met$met_loc[index]	<-	"Thyroid"
    
    
    write.csv(THCA_met, file = str_glue("~/storage/PanCancerAnalysis/TCGABiolinks/metastatic_clin_info/{i}_metastatic_staus_.csv"))
    
    print("THCA Done")
    rm(THCA_met)
    
    
    
  } 
  
}

# R version 3.5.1 (2018-07-02)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.19.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C              LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pbapply_1.4-2                            EBImage_4.24.0                           remotes_2.1.0                            gsubfn_0.7                              
# [5] proto_1.0.0                              forcats_0.4.0                            purrr_0.3.3                              readr_1.3.1                             
# [9] tidyr_1.0.2                              tibble_2.1.3                             tidyverse_1.3.0                          targetscan.Mm.eg.db_0.6.1               
# [13] data.table_1.12.8                        RColorBrewer_1.1-2                       gplots_3.0.1.2                           Glimma_1.10.1                           
# [17] edgeR_3.24.3                             limma_3.38.3                             vsn_3.50.0                               pheatmap_1.0.12                         
# [21] apeglm_1.4.2                             EnsDb.Mmusculus.v79_2.99.0               ensembldb_2.6.8                          AnnotationFilter_1.6.0                  
# [25] TxDb.Mmusculus.UCSC.mm10.knownGene_3.4.4 GenomicFeatures_1.34.8                   tximport_1.10.1                          org.Mm.eg.db_3.7.0                      
# [29] rpart_4.1-13                             randomForest_4.6-14                      keras_2.2.5.0                            kernlab_0.9-29                          
# [33] caret_6.0-85                             lattice_0.20-38                          plyr_1.8.5                               stringr_1.4.0                           
# [37] dplyr_0.8.4                              WGCNA_1.68                               fastcluster_1.1.25                       dynamicTreeCut_1.63-1                   
# [41] org.Hs.eg.db_3.7.0                       AnnotationDbi_1.44.0                     clusterProfiler_3.10.1                   EnhancedVolcano_1.0.1                   
# [45] ggrepel_0.8.1                            ggplot2_3.2.1                            DESeq2_1.22.2                            SummarizedExperiment_1.12.0             
# [49] DelayedArray_0.8.0                       BiocParallel_1.16.6                      matrixStats_0.55.0                       Biobase_2.42.0                          
# [53] GenomicRanges_1.34.0                     GenomeInfoDb_1.18.2                      IRanges_2.16.0                           S4Vectors_0.20.1                        
# [57] BiocGenerics_0.28.0                     
# 
# loaded via a namespace (and not attached):
# [1] rtracklayer_1.42.2       ModelMetrics_1.2.2.1     coda_0.19-3              acepack_1.4.1            bit64_0.9-7              knitr_1.27               RCurl_1.98-1.1           doParallel_1.0.15       
# [9] generics_0.0.2           preprocessCore_1.44.0    cowplot_1.0.0            RSQLite_2.2.0            europepmc_0.3            tensorflow_2.0.0         bit_1.1-15.1             enrichplot_1.2.0        
# [17] xml2_1.2.2               lubridate_1.7.4          assertthat_0.2.1         viridis_0.5.1            gower_0.2.1              xfun_0.12                hms_0.5.3                DEoptimR_1.0-8          
# [25] fansi_0.4.1              progress_1.2.2           readxl_1.3.1             dbplyr_1.4.2             caTools_1.17.1.3         igraph_1.2.4.2           DBI_1.1.0                geneplotter_1.60.0      
# [33] htmlwidgets_1.5.1        ellipsis_0.3.0           backports_1.1.5          annotate_1.60.1          biomaRt_2.38.0           vctrs_0.2.2              abind_1.4-5              withr_2.1.2             
# [41] ggforce_0.3.1            triebeard_0.3.0          robustbase_0.93-5        bdsmatrix_1.3-4          checkmate_1.9.4          GenomicAlignments_1.18.1 prettyunits_1.1.1        cluster_2.0.7-1         
# [49] DOSE_3.8.2               lazyeval_0.2.2           crayon_1.3.4             genefilter_1.64.0        recipes_0.1.9            pkgconfig_2.0.3          tweenr_1.0.1             nlme_3.1-137            
# [57] ProtGenerics_1.14.0      nnet_7.3-12              rlang_0.4.4              lifecycle_0.1.0          affyio_1.52.0            modelr_0.1.5             cellranger_1.1.0         polyclip_1.10-0         
# [65] tiff_0.1-5               Matrix_1.2-15            urltools_1.7.3           reprex_0.3.0             base64enc_0.1-3          whisker_0.4              ggridges_0.5.2           png_0.1-7               
# [73] viridisLite_0.3.0        bitops_1.0-6             KernSmooth_2.23-15       pROC_1.16.1              Biostrings_2.50.2        blob_1.2.1               qvalue_2.14.1            robust_0.4-18.2         
# [81] jpeg_0.1-8.1             gridGraphics_0.4-1       scales_1.1.0             memoise_1.1.0            magrittr_1.5             gdata_2.18.0             zlibbioc_1.28.0          compiler_3.5.1          
# [89] bbmle_1.0.23.1           rrcov_1.5-2              Rsamtools_1.34.1         cli_2.0.1                affy_1.60.0              XVector_0.22.0           htmlTable_1.13.3         Formula_1.2-3           
# [97] MASS_7.3-51.1            tidyselect_1.0.0         stringi_1.4.5            yaml_2.2.1               emdbook_1.3.11           GOSemSim_2.8.0           locfit_1.5-9.1           latticeExtra_0.6-28     
# [105] grid_3.5.1               fastmatch_1.1-0          tools_3.5.1              rstudioapi_0.10          foreach_1.4.7            foreign_0.8-71           gridExtra_2.3            prodlim_2019.11.13      
# [113] farver_2.0.3             ggraph_2.0.0             digest_0.6.23            rvcheck_0.1.7            BiocManager_1.30.10      lava_1.6.6               Rcpp_1.0.3               broom_0.5.4             
# [121] httr_1.4.1               colorspace_1.4-1         rvest_0.3.5              fs_1.3.1                 XML_3.99-0.3             reticulate_1.14          splines_3.5.1            graphlayouts_0.5.0      
# [129] ggplotify_0.0.4          fit.models_0.5-14        xtable_1.8-4             jsonlite_1.6.1           tidygraph_1.1.2          timeDate_3043.102        UpSetR_1.4.0             zeallot_0.1.0           
# [137] ipred_0.9-9              R6_2.4.1                 Hmisc_4.3-0              pillar_1.4.3             htmltools_0.4.0          glue_1.3.1               fftwtools_0.9-8          class_7.3-14            
# [145] codetools_0.2-15         fgsea_1.8.0              pcaPP_1.9-73             mvtnorm_1.0-12           numDeriv_2016.8-1.1      curl_4.3                 tfruns_1.4               gtools_3.8.1            
# [153] GO.db_3.7.0              survival_3.1-8           munsell_0.5.0            DO.db_2.9                GenomeInfoDbData_1.2.0   iterators_1.0.12         impute_1.56.0            haven_2.2.0             
# [161] reshape2_1.4.3           gtable_0.3.
# 
