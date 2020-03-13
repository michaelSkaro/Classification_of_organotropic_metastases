# metastatic classification
library(ggplot2)
library(ggrepel)
dat <- data.table::fread("Results_summary_metastatic_info.txt", header = TRUE, stringsAsFactors = TRUE)

dat.pos <- dat[dat$Class == "positive",]
dat.neg <- dat[dat$Class == "negative",]

p <- ggplot2::ggplot(dat, aes(`FP Rate`, `TP Rate`, label = Project, color = Class)) +
  geom_point(color ="black") + geom_text_repel(size =3) + geom_abline(intercept = 0, slope = 1) + theme_classic() + xlim(0, 1) + ylim(0, 1) 


dat <- data.table::fread("classification_accuracy_mets.txt", header = TRUE, stringsAsFactors = TRUE)

p <- ggplot2::ggplot(dat, aes(Percent, Cases, label = Project, color = Class)) +
  geom_point(color ="black") + geom_text_repel(size =3) + geom_abline(intercept = 0, slope = 1) 
+ theme_classic() + xlim(0, 1) + ylim(0, 1) 




