
library(BiocInstaller) # Load the BiocInstaller packag
library(BiocManager) 
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

# Check the alphabet of the zikaVirus
alphabet(zikaVirus)

# Check the alphabetFrequency of the zikaVirus
alphabetFrequency(zikaVirus)
# Check alphabet of the zikaVirus using baseOnly = TRUE
alphabet(zikaVirus, baseOnly = TRUE)
readDNAStringSet(zikaVirus)

# Unlist the set, select the first 21 letters, and assign to dna_seq
dna_seq <- subseq(unlist(zikaVirus), end = 21)
dna_seq

# Transcribe dna_seq into an RNAString object and print it
rna_seq <- RNAString(dna_seq) 

# Translate rna_seq into an AAString object and print it
aa_seq <- translate(rna_seq)

# read the sequence as a set
zikaVirus <- readDNAStringSet("data/zika.fa")
# Create zikv with one collated sequence using `zikaVirus`
zikv <- unlist(zikaVirus)

# Check the length of zikaVirus and zikv
length(zikaVirus)
length(zikv)

# Reverse the zikv sequence
reverse(zikv)

# Complement the zikv sequence
complement(zikv)

# Reverse complement the zikv sequence
reverseComplement(zikv)

# Translate the zikv sequence
translate(zikv)

# For Sets
vmatchPattern(pattern = "ACATGGGCCTACCATGGGAG", 
              subject = zikaVirus, max.mismatch = 1)
# For single sequences
matchPattern(pattern = "ACATGGGCCTACCATGGGAG", 
              subject = zikv, max.mismatch = 1)


findPalindromes(zikv) # Find palindromes in zikv


rnaframesZikaSet # print the rnaframesZikaSet
AAzika6F <- translate(rnaframesZikaSet) # translate all 6 reading frames 
vcountPattern(pattern = NS5, subject = AAzika6F, max.mismatch = 15) # Count the matches allowing 15 mistmatches
selectedSet <- AAzika6F[3] # Select the frame that contains the match
selectedSeq <- unlist(selectedSet) #Convert this frame into a single sequence
matchPattern(ns5, selectedSeq, max.mismatch = 15)
# load package IRanges

library(IRanges)
IRnum1 <- IRanges(start = c(1:5), end = 100) # start vector 1 through 5, end 100 
IRnum2 <- IRanges(end = 100, width = c(89,10)) # end 100 and width 89 and 10

library(GenomicRanges)
myGR <- as(seq_intervals,"GRanges")
seqnames(gr)
ranges(gr)
mcols(gr)
genome(gr)
seqinfo(gr)
# Load human reference genome hg38
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assign hg38 to hg, then print it
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg
# Extract all positive stranded genes in chromosome X, assign to hg_chrXgp, then sort it
hg_chrXgp <- genes(hg, filter = list(tx_chrom = "chrX", tx_strand = "+"))
sort(hg_chrXgp)

rangefound <- subsetByOverlaps(hg_chrX, ABCD1) # Store the overlapping range in rangefound
names(rangefound) # Check names of rangefound
ABCD1 # Check the gene of interest 
rangefound # Check rangefound

library(TxDb.Hsapiens.UCSC.hg38.knownGene) # Load the human transcripts DB to hg
hg <- TxDb.Hsapiens.UCSC.hg38.knownGene
seqlevels(hg) <- c("chrX") # Prefilter chromosome X "chrX" using seqlevels()
hg_chrXt <- transcriptsBy(hg, by = "gene") # Get all transcripts by gene and print it
hg_chrXt$`215` # Select gene `215` from the transcripts

library(ShortRead)

fqsample <- readFastq(dirPath = "data/", pattern = "fastq") # read fastq

writeFastq(fqsample, file = "data/sample.fastq.gz")

# Subsample of 500 bases
sampler <- FastqSampler(con="data/SRR1971253.fastq", n=500)
# save the yield of 500 read sequences
sample_small <- yield(sampler)

sread(fqsample)[1]
quality(fqsample)[1] # Quality is represented with ASCII characters 
pq <- PhredQuality(quality(fqsample)) ## PhredQuality instance
qs <- as(pq, "IntegerList") # transform encoding into scores 

qaSummary <- qa(fqsample, lane = 1)
browseURL(report(qaSummary))

# Check quality
quality(fqsample)

# Check encoding of quality
encoding(quality(fqsample))

# Check baseQuality
qaSummary[["baseQuality"]]


# Glimpse nucByCycle
glimpse(nucByCycle)

# Create a line plot of cycle vs count
nucByCycle %>% 
  # Gather the nucleotide letters in alphabet and get a new count column
  gather(key = alphabet, value = count , -cycle) %>% 
  ggplot(aes(x = cycle, y =  count, color = alphabet)) +
  geom_line(size = 0.5 ) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank())

# Sample with duplicates of class: ShortReadQ
dfqsample

# Get the reads from dfqsample
mydReads <- sread(dfqsample)

# Create myFil using polynFilter
myFil <- polynFilter(threshold = 3, nuc = c("A"))

# Apply your filter to fqsample
filterCondition <- myFil(fqsample)

# Use myFil with fqsample
filteredSequences <- fqsample[filterCondition==TRUE]

# Check reads of filteredSequences
sread(filteredSequences)


library(Rqc)
files <- # get the full path of the files you want to assess
qaRqc <- rqcQA(files)
qaRqc <- rqcQA(files, workers = 4))
# sample of sequences
set.seed(1111)

qaRqc_sample <- rqcQA(files, workers = 4, sample = TRUE, n = 500))
# create a report
reportFile <- rqcReport(qaRqc, templateFile = "myReport.Rmd")
browseURL(reportFile)
#The class of qaRqc is rqcResultSet
methods(class = "RqcResultSet")
