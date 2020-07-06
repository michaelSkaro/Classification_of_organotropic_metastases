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

plat <- platforms[3]
for(plat in platforms){
  # only use this block if the GPL == GPL1390
  if(plat == "GPL1390"){
    setwd("~/storage/Metastatic_Organo_Tropism")
    # subset by platforms
    coldata <- refDat %>%
      filter(Platform_id == plat)
    # reduce download to only desired GSM list
    GSMs <- unique(coldata$Sample_id)
    ##download gsms and concatenate
    genesSets <- GSMs[1] 
    # download samples
    gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL1390.annot")
    var <- as.data.frame(gsm@dataTable@table)%>%
      dplyr::select(ID_REF, VALUE)
    var$Platform_SPOTID <- gpl@dataTable@table$Platform_SPOTID
    var <- var[c("ID_REF" , "Platform_SPOTID", "VALUE")]
    colnames(var) <- c("ID_REF" , "Platform_SPOTID", genesSets)  
    dat <- var
    # incriment the genesSets
    
  # iterate over the  rest of the GSM list, download each one, bind each one to make a 
  # values matrix
    #incriment the genesSets because we started at the first one for initial download
    for(genesSets in GSMs[2:188]){
      
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      
      #clear the temp_dircache
      tmp_dir <- tempdir()
      list.files(tmp_dir)
      files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
      file.exists(files)
      file.remove(files)
      file.exists(files)
      list.files(tmp_dir)
      
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      colnames(var) <- c("ID_REF" , genesSets)
      dat <- left_join(dat, var, by= "ID_REF")
    }
  }
    # write dat to downloads folder
    
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    # incriment to the next platform because the stupid annot gpl doesnt work
    
    if(plat == "GPL887"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL887.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      var$Platform_SPOTID <- gpl@dataTable@table$Platform_SPOTID
      var <- var[c("ID_REF" , "Platform_SPOTID", "VALUE")]
      colnames(var) <- c("ID_REF" , "Platform_SPOTID", genesSets)  
      dat <- var
    
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL570"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL570.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      var$Platform_Probe_ID <- gpl@dataTable@table$ID
      var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL96"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL96.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      var$Platform_Probe_ID <- gpl@dataTable@table$ID
      var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    # if(plat == "GPL6370"){
    #   setwd("~/storage/Metastatic_Organo_Tropism")
    #   # subset by platforms
    #   coldata <- refDat %>%
    #     filter(Platform_id == plat)
    #   # reduce download to only desired GSM list
    #   GSMs <- unique(coldata$Sample_id)
    #   ##download gsms and concatenate
    #   genesSets <- GSMs[1] 
    #   # download samples
    #   gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #   gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6370.txt")
    #   var <- as.data.frame(gsm@dataTable@table)%>%
    #     dplyr::select(ID_REF, VALUE)
    #   var$Platform_Probe_ID <- gpl@dataTable@table$ID
    #   var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
    #   colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
    #   dat <- var
    #   
    #   for(genesSets in GSMs[2:length(GSMs)]){
    #     
    #     gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #     
    #     #clear the temp_dircache
    #     tmp_dir <- tempdir()
    #     list.files(tmp_dir)
    #     files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
    #     file.exists(files)
    #     file.remove(files)
    #     file.exists(files)
    #     list.files(tmp_dir)
    #     
    #     var <- as.data.frame(gsm@dataTable@table)%>%
    #       dplyr::select(ID_REF, VALUE)
    #     colnames(var) <- c("ID_REF" , genesSets)
    #     dat <- left_join(dat, var, by= "ID_REF")
    #   }
    # }
    # 
    # write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    # write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    # rm(dat)
    
    #Skip bead chip GPL6370
    
    # if(plat == "GPL8128"){
    #   setwd("~/storage/Metastatic_Organo_Tropism")
    #   # subset by platforms
    #   coldata <- refDat %>%
    #     filter(Platform_id == plat)
    #   # reduce download to only desired GSM list
    #   GSMs <- unique(coldata$Sample_id)
    #   ##download gsms and concatenate
    #   genesSets <- GSMs[1] 
    #   # download samples
    #   gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #   gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL8128.annot")
    #   var <- as.data.frame(gsm@dataTable@table)%>%
    #     dplyr::select(ID_REF, VALUE)
    #   var$Platform_Probe_ID <- gpl@dataTable@table$ID
    #   var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
    #   colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
    #   dat <- var
    #   
    #   for(genesSets in GSMs[2:length(GSMs)]){
    #     
    #     gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #     
    #     #clear the temp_dircache
    #     tmp_dir <- tempdir()
    #     list.files(tmp_dir)
    #     files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
    #     file.exists(files)
    #     file.remove(files)
    #     file.exists(files)
    #     list.files(tmp_dir)
    #     
    #     var <- as.data.frame(gsm@dataTable@table)%>%
    #       dplyr::select(ID_REF, VALUE)
    #     colnames(var) <- c("ID_REF" , genesSets)
    #     dat <- left_join(dat, var, by= "ID_REF")
    #   }
    # }
    
    # write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    # write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    # rm(dat)
    
    
    if(plat == "GPL97"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL97.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL91"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL91.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    
    if(plat == "GPL6102"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6102.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL80"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL80.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL8469"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL8469.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL13140"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL13140.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL5858"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL5858.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL334"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL334.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL6480"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6480.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL6244"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6244.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL8432"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL8432.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    if(plat == "GPL11350"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6480.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL3307"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6480.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL4133"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL4133.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL15236"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL15236.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL13703"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL13703.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL5826"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL5826.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL3282"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL3282.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL13497"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL13497.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    if(plat == "GPL3985"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL3985.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL7722"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL7722.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL5918"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL5918.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL10379"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL10379.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL6884"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL6884.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL10558"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL10558.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL8178"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL8178.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL10262"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL10379.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL17077"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL17077.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL16744"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL16744.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL17537"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL17537.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    # if(plat == "GPL11154"){
    #   setwd("~/storage/Metastatic_Organo_Tropism")
    #   # subset by platforms
    #   coldata <- refDat %>%
    #     filter(Platform_id == plat)
    #   # reduce download to only desired GSM list
    #   GSMs <- unique(coldata$Sample_id)
    #   ##download gsms and concatenate
    #   genesSets <- GSMs[1] 
    #   # download samples
    #   gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #   #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL10379.annot")
    #   var <- as.data.frame(gsm@dataTable@table)%>%
    #     dplyr::select(ID_REF, VALUE)
    #   #var$Platform_Probe_ID <- gpl@dataTable@table$ID
    #   #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
    #   #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
    #   colnames(var) <- c("ID_REF", genesSets)  
    #   dat <- var
    #   
    #   for(genesSets in GSMs[2:length(GSMs)]){
    #     
    #     gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
    #     
    #     #clear the temp_dircache
    #     tmp_dir <- tempdir()
    #     list.files(tmp_dir)
    #     files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
    #     file.exists(files)
    #     file.remove(files)
    #     file.exists(files)
    #     list.files(tmp_dir)
    #     
    #     var <- as.data.frame(gsm@dataTable@table)%>%
    #       dplyr::select(ID_REF, VALUE)
    #     colnames(var) <- c("ID_REF" , genesSets)
    #     dat <- left_join(dat, var, by= "ID_REF")
    #   }
    # }
    # 
    # write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    # write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    # rm(dat)
    
    if(plat == "GPL5175"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL5175.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL14951"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL14951.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL16686"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL16686.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL8300"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL8300.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL92"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      #gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL92.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL93"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL93.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL2891"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL2891.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL15018"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL15018.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL17537-2"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL17537.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL15018-2"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL15018.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL10332"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL10332.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL20712"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL20712.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL20945"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL20945.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL11010"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL11010.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    if(plat == "GPL1708-2"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL1708.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL15659"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL15659.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL16699"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL16699.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL14550"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL14550.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL18044"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL18044.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    if(plat == "GPL22330"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL22330.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL201"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL201.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    if(plat == "GPL5049"){
      setwd("~/storage/Metastatic_Organo_Tropism")
      # subset by platforms
      coldata <- refDat %>%
        filter(Platform_id == plat)
      # reduce download to only desired GSM list
      GSMs <- unique(coldata$Sample_id)
      ##download gsms and concatenate
      genesSets <- GSMs[1] 
      # download samples
      gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
      gpl <- getGEO(filename="~/storage/Metastatic_Organo_Tropism/annot_files/GPL5049.annot")
      var <- as.data.frame(gsm@dataTable@table)%>%
        dplyr::select(ID_REF, VALUE)
      #var$Platform_Probe_ID <- gpl@dataTable@table$ID
      #var <- var[c("ID_REF" , "Platform_Probe_ID", "VALUE")]
      #colnames(var) <- c("ID_REF" , "Platform_Probe_ID", genesSets)  
      colnames(var) <- c("ID_REF", genesSets)  
      dat <- var
      
      for(genesSets in GSMs[2:length(GSMs)]){
        
        gsm <- GEOquery::getGEO(genesSets, AnnotGPL = TRUE,  getGPL = TRUE)
        
        #clear the temp_dircache
        tmp_dir <- tempdir()
        list.files(tmp_dir)
        files <- list.files(tmp_dir, full.names = T, pattern = ".soft")
        file.exists(files)
        file.remove(files)
        file.exists(files)
        list.files(tmp_dir)
        
        var <- as.data.frame(gsm@dataTable@table)%>%
          dplyr::select(ID_REF, VALUE)
        colnames(var) <- c("ID_REF" , genesSets)
        dat <- left_join(dat, var, by= "ID_REF")
      }
    }
    
    write.csv(dat, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_raw_values.csv"))  
    write.csv(coldata, file = str_glue("~/storage/Metastatic_Organo_Tropism/results/{plat}_design_matrix_raw.csv"))  
    rm(dat)
    
    
    # get anno data  done :)
    
    # intput DESeq2 for the TCGA data
    
    
    
    
    
    
    
    # add a few methods today
    

    
    
    
    
    
    
    # analyze the data for the expression levels of the 
    
    
}









# 
# eset <- GDS2eSet(gsm,do.log2=TRUE)
# var <- as.data.frame(eset@assayData$exprs) %>%
#   tibble::rownames_to_column("Probes")
# dat <- var
# annotGPL <- pData(eset)



# get each model GSE
##################################################################
#GSE10893
#GPL1390

#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE10893", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("03000000000030000000300000011000000330000000004000",
               "10000000000000000000000000000302222222220000000000",
               "0004440000000000000000000000000000000000000000X100",
               "00000000000000000000000001110000011")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G4-G0, G1-G0, G2-G1, G3-G2, G4-G3, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE10893", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1390", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- paste0("03000000000030000000300000011000000330000000004000",
               "10000000000000000000000000000302222222220000000000",
               "0004440000000000000000000000000000000000000000X100",
               "00000000000000000000000001110000011")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
#set group names
sml <- paste("G", sml, sep="")  

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("Breast+Cancer","Normal+Breast","Treated+Breast+Cancer","Breast+Cancer+Brain+Met","Breast+Cancer+Lymph+Met")

# set parameters and draw the plot
palette(c("#f4dfdf","#dfeaf4","#dfeaf4","#f4dfdf","#f2cb98", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE10893", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")








# get the libraries
library(Biobase)
library(oligoClasses)
library(knitr)
library(BiocStyle)
library(oligo)
library(geneplotter)
library(ggplot2)
library(dplyr)
library(LSD)
library(gplots)
library(RColorBrewer)
library(ArrayExpress)
library(arrayQualityMetrics)
library(stringr)
library(matrixStats)
library(topGO)
library(genefilter)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(pheatmap)
library(mvtnorm)
library(DAAG)
library(multcomp)
library(limma)
library(ReactomePA)
library(clusterProfiler)
library(openxlsx)
library(devtools)
library(biomaRt)
library(ArrayExpress)


