#NA's to 0
d[is.na(d)] <- 0

#Pull function from other R script
source("~/Desktop/R_scripts/universal_names_function.R")

#Substitute exact match
gsub("\\<AC093012.1\\>","IRAK4",external_gene_name)
