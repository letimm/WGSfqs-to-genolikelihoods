#!/usr/bin/env Rscript

packages_needed <- c("stringr")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

args=commandArgs(trailingOnly = TRUE)
FILENAME <- args[1]
ALPHA <- args[2]
lr <- read.table(FILENAME) # read in ngsParalog calcLR output
lr$pval <- 0.5*pchisq(lr$V5,df=1,lower.tail=FALSE) # append column of p-values
lr$pval.adj <- p.adjust(lr$pval, method="bonferroni") # p-values adjusted for number of tested sites

#full_report_filename <- str_replace(FILENAME, ".lr", "_sig.out")
#write.table(lr, file = full_report_filename, row.names=FALSE, sep="\t")

homologous_sites <- lr[which(lr$pval.adj > ALPHA),1:2]
retain_sites_filename <- str_replace(FILENAME, ".lr", "_retain.sites")
write.table(homologous_sites, file = retain_sites_filename, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
