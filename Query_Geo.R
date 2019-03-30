# download metastatic tissues exoression data

BiocManager::install("GEOquery", version = "3.8")
library(GEOquery)

# gpl97 <- getGEO('GPL97')
# Meta(gpl97)$title

MetDB_annot <- data.table::fread("/home/mskaro1/storage/Metastatic_Organo_Tropism/Metastatic_database_project_information.csv")

Geo_platform_list<- unique(MetDB_annot$Platform_id)
platform <- Geo_platform_list[1]

res<- data.frame()
for (platform in Geo_platform_list) {
  tmp <- getGEO(platform)
  
  foo<- Meta(tmp)$title
  
  res <- rbind.data.frame(res, foo)
}

write.csv(res, "probably_did_not_work.csv")
