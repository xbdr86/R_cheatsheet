#NA's to 0
df[is.na(df)] <- 0

#Delete NA rows
df <- df[complete.cases(df), ]

#Pull function from other R script
source("~/Desktop/R_scripts/universal_names_function.R")

#Substitute exact match
gsub("\\<AC093012.1\\>","IRAK4",external_gene_name)
