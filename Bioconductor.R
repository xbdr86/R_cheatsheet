
library(BiocInstaller) # Load the BiocInstaller packag
version # Check R version

BiocInstaller::biocVersion() # Check Bioconductor version (For versions <= 3.7)
# or biocVersion()

library(BSgenome) #Load the BSgenome package
packageVersion("BSgenome") #Check the version of the BSgenome package
