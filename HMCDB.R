# data analysis from data download 

library(GEOquery)

setwd("~/storage/Metastatic_Organo_Tropism")
#Metastatic Disease datasheet

refDat <- data.table::fread("Metastatic_database_project_information.csv")

setwd("~/storage/Metastatic_Organo_Tropism/downloads")


files.list <- list.files()

coldata <- data.table::fread("GPL1390_design_matrix_raw.csv")

coldata <- coldata %>%
  dplyr::select(-"V1")

dat <- data.table::fread("GPL1390_raw_values.csv")

id.map <- dat %>%
  dplyr::select(c("ID_REF", "Platform_SPOTID"))

dat <- dat %>%
  dplyr::select(-c("V1", "ID_REF", "Platform_SPOTID"))

coldata <- coldata[order(coldata$Sample_id),]


coldata<- coldata[match(unique(coldata$Sample_id), coldata$Sample_id),]

design.mat <- coldata

design.mat$Sample_label <- as.character(design.mat$Sample_label)


design.mat$sample.levels <- as.factor(design.mat$Sample_label)
design.mat$sample.levels <- relevel(design.mat$sample.levels, ref="Primary Normal")


table(design.mat$sample.levels)

# Visualize the density plots of the matrix. Show data has been normalized
level.pal <- brewer.pal(4, "Dark2")
level.cols <- level.pal[unname(design.mat$sample.levels)]

plotDensities(dat, legend=F, col=level.cols , main="Arrays Normalized")
legend("topright", legend=levels(design.mat$sample.levels), fill=level.pal)

# hierarchical clustering

annot.dat <- data.table::fread("GPL1390.annot")



dat <- as.data.frame(dat)


# heatmaps

v <- colnames(dat) %in% coldata$Sample_id
w <- coldata$Sample_id %in% colnames(dat)
colnames(dat)[v] <- coldata$Sample_label[w]

dat[is.na(dat)] <- 0

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(dat, name = "Expression", show_row_names = FALSE, cluster_columns = FALSE, column_order = c("Metastasis Tumor","Primary Normal","Primary Tumor"))

design.mat <- design.mat[design.mat$Sample_id %in% colnames(dat),]


# new designmatrix with binary operators

design.mat$sample.levels <- str_replace_all(design.mat$sample.levels,"[[:punct:]]","_")
design.mat$sample.levels <- str_replace_all(design.mat$sample.levels," ","_")
design.mat$sample.levels <- str_replace_all(design.mat$sample.levels,"__","_")


my.design <- model.matrix(~0 + sample.levels, design.mat)

colnames(my.design) <- c("Metastasis_Tumor", "Primary_Normal", "Primary_Tumor")
colnames(my.design)
rownames(my.design) <- design.mat$sample.levels

# linear model fitting

my.fit <- lmFit(dat, my.design)

fit2 <- eBayes(my.fit, 0.01)
x <- c("Metastasis_Tumor-Primary_Normal","Primary_Tumor-Primary_Normal","Primary_Tumor-Metastasis_Tumor")
cont.matrix <- makeContrasts(contrasts=x,levels=c("Metastasis_Tumor", "Primary_Normal", "Primary_Tumor"))
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number = 22576)

View(tT)
