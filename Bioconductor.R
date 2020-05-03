
library(BiocInstaller) # Load the BiocInstaller packag
version # Check R version

BiocInstaller::biocVersion() # Check Bioconductor version (For versions <= 3.7)
# or biocVersion()

library(BSgenome) #Load the BSgenome package
packageVersion("BSgenome") #Check the version of the BSgenome package
available.genomes()
# Load the yeast genome
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Assign data to the yeastGenome object
yeastGenome <- BSgenome.Scerevisiae.UCSC.sacCer3
# Get the head of seqnames and tail of seqlengths for yeastGenome
head(seqnames(yeastGenome))
tail(seqlengths(yeastGenome))

# Print chromosome M, alias chrM
print(yeastGenome$chrM)

# Count characters of the chrM sequence
nchar(yeastGenome$chrM)

getSeq(yeastGenome, names = "chrI", start = 100, end = 150)
getSeq(yeastGenome, end=30) # Get the first 30 bases of each chromosome
