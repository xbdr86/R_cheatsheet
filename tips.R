# R cheatsheets
https://www.rstudio.com/resources/cheatsheets/?s=03
# Create, Move and delete files and folders
https://www.r-bloggers.com/2018/07/r-how-to-create-delete-move-and-more-with-files/?s=03 
directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE) #list of directories underneath
file.exists("matrix_exactMatches_RPM_all.tsv") #check if file exists
fileNames <- Sys.glob("*.csv")

#Print cat output
cat(paste("Line 1","Line 2",sep="\n"))

#Excel
#install.packages("tidyverse")
#install.packages("readxl")
library(readxl)
WT <- read_excel("journal.ppat.1010464.s009.xlsx", sheet = "WT")
##older
require(gdata)
excel_file <- read.xls(file_name, sheet = "Sheet1", header = TRUE)

#fread
install.packages('data.table')
require(data.table)
total2 <- unique(fread(file = "data.tsv", header = T))

#Count instances
read_dist <- as.data.frame(table(master_file$reads))

#RNA logos in R
# https://omarwagih.github.io/ggseqlogo/
#install.packages("ggseqlogo")
require(ggseqlogo)
ggseqlogo(ap1_st$variant )
ggplot() + geom_logo(ap1_st$variant) + theme_logo()

#Anchors start and end positions in string 
^start
send$
apic1_5 <- "^\\Q((((((..(((.((((((((((((.((((((((.(.((.\\E"
apic1_3 <- "\\Q.)).).)))))))).)).)))))))))).)))..))))))\\E$"

#reverse complementary
test <- "AAUGC"
library(seqinr)
toupper(c2s(rev(comp(s2c(gsub("U","T",test))))))

#reverse sequence
splits <- strsplit(dna, "")[[1]]
reversed <- rev(splits)
final_result <- paste(reversed, collapse = "")


#mode for factors
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
name_selected <- as.character(calculate_mode(as.factor(WT$miR_name)))

#string distance calculation
require('stringdist')
potential$dist5 <- stringdist(potential$part5, part5ref, method="dl")

#NA's to 0
df[is.na(df)] <- 0
df$var1[is.na(df$var1)] <- 0

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
mini <- mini[order(-mini$WT_GFP),]

#downsampling
set.seed(123)
index <- sample(x = 1:nrow(tsne), size = 10000)
small <- tsne[index, ]
plot(small$V2, small$V3)

#compare two sequences nucleotide by nucleotide
found_seq <- unlist(strsplit(found_seq, split = ""))
ref_seq <- unlist(strsplit(ref_seq, split = ""))
test_table <- data.frame(found_seq,ref_seq)
test_table$equal <- test_table$found_seq==test_table$ref_seq

#time keeping
start.time <- Sys.time()
end.time <- Sys.time()
gene_names$time[a] <- difftime(end.time, start.time, units = "mins")

# http://r-statistics.co/Strategies-To-Improve-And-Speedup-R-Code.html
# https://www.data-to-viz.com/caveat/overplotting.html?s=03 overplotting
# https://youtu.be/B880KgTpG0Q 
