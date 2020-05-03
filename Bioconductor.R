
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
