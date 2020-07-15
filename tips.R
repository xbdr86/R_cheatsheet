fileNames <- Sys.glob("*.csv")

#NA's to 0
df[is.na(df)] <- 0

#Delete NA rows
df <- df[complete.cases(df), ]

#Pull function from other R script
source("~/Desktop/R_scripts/universal_names_function.R")

#Substitute exact match
gsub("\\<AC093012.1\\>","IRAK4",external_gene_name)

#All combinations
card_numbers <- c(c(2:10),"J","Q","K","A")
suits <- c("Hearts","Tiles","Clovers","Pikes")
deck <- expand.grid(card_numbers,suits)

#not in function
`%notin%` <- Negate(`%in%`)


#read a text
library(readtext)
gene <- readtext("DROSHA_cDNA.txt")

#order
newdata <- mtcars[order(mpg, -cyl),]

#downsampling
set.seed(123)
index <- sample(x = 1:nrow(tsne), size = 10000)
small <- tsne[index, ]
plot(small$V2, small$V3)
