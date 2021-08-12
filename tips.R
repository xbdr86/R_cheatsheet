# R cheatsheets
https://www.rstudio.com/resources/cheatsheets/?s=03

#Excel
require(gdata)
excel_file <- read.xls(file_name, sheet = "Sheet1", header = TRUE) 

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

#mode for factors
calculate_mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}
name_selected <- as.character(calculate_mode(as.factor(WT$miR_name)))

#string distance calculation
require('stringdist')
potential$dist5 <- stringdist(potential$part5, part5ref, method="dl")


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

#time keeping
start.time <- Sys.time()
end.time <- Sys.time()
gene_names$time[a] <- difftime(end.time, start.time, units = "mins")

# http://r-statistics.co/Strategies-To-Improve-And-Speedup-R-Code.html
# https://www.data-to-viz.com/caveat/overplotting.html?s=03 overplotting
# https://youtu.be/B880KgTpG0Q 
