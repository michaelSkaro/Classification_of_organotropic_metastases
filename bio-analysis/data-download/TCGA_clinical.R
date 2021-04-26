library(TCGAbiolinks)
library(xml2)
library(tidyverse)

projects <- c("TCGA-ACC","TCGA-BLCA","TCGA-BRCA","TCGA-CESC","TCGA-CHOL","TCGA-COAD",
              "TCGA-DLBC","TCGA-ESCA","TCGA-GBM","TCGA-HNSC","TCGA-KICH","TCGA-KIRC",
              "TCGA-KIRP","TCGA-LAML","TCGA-LGG","TCGA-LIHC","TCGA-LUAD","TCGA-LUSC",
              "TCGA-MESO","TCGA-OV","TCGA-PAAD","TCGA-PCPG","TCGA-PRAD","TCGA-READ",
              "TCGA-SARC","TCGA-SKCM","TCGA-STAD","TCGA-TGCT","TCGA-THCA","TCGA-THYM",
              "TCGA-UCEC","TCGA-UCS","TCGA-UVM")


for (proj in projects) {
  query <- GDCquery(project = proj,
                    data.category = "Clinical",
                    file.type = "xml")
  GDCdownload(query)
}

# Parse xml files for each cancer type ---------------------------------------
projects.downloaded <- list.files("./GDCdata/")
proj <- projects.downloaded[1]
for (proj in projects.downloaded) {
  
  patients <- list.files(
    stringr::str_glue("./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/"),
    pattern = "xml$", recursive = T)
  filenames <- stringr::str_glue(
    "./GDCdata/{proj}/harmonized/Clinical/Clinical_Supplement/{patients}"
  )
  patients <- map_chr(str_split(patients, "/"), function(x) x[[1]])
  
  res <- NULL
  for (i in seq_along(filenames)) {
    df <- read_xml(filenames[i])
    barcode <- df %>%
      xml_find_first("//shared:bcr_patient_barcode") %>%
      xml_text()
    
    # In project:patient, look for all fields with text
    df <- xml_find_all(df, "//*[local-name(.) = 'patient']//*[text()]")
    df <- tibble(
      fileID = patients[i],
      field = map_chr(df, xml_path),
      value = map_chr(df, xml_text)
    ) %>%
      mutate(field = str_extract(field, "[^:]+$"))
    res <- bind_rows(res, df)
  }
  
  # The arguments to spread():
  # - data: Data object
  # - key: Name of column containing the new column names
  # - value: Name of column containing values
  
  fields <- c("fileID", unique(res$field))
  res <- res %>%
    distinct() %>%
    group_by(fileID,field) %>% 
    dplyr::mutate(group_row = 1:n()) %>%
    spread(key = field, value = value) %>%
    select(-group_row)
  
  res <- res[match(unique(res$fileID), res$fileID),]
  data.table::fwrite(res, str_glue("/mnt/storage/mskaro1/Clinical_annotation_{proj}.csv"))
}
