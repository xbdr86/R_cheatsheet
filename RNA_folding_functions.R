getRNADuplex<-function(ref_seq,target_seq){
  ref_seq <- gsub("T","U",ref_seq)
  target_seq <- gsub("T","U",target_seq)
  input_data<-c(paste('>ref', sep=''), as.character(ref_seq), paste('>test', sep=''), as.character(target_seq))
  input_file<-'temp-miRNA-seq.fa'
  unlink(input_file)
  file<-file(input_file)
  writeLines(input_data, file)
  close(file)
  energy<-system('RNAduplex < temp-miRNA-seq.fa', intern=T)[3]
  #energy <- substr(energy,nchar(energy)-10,nchar(energy))
  #energy<-gsub("\\(", "\\*", energy)
  #energy<-gsub("\\)", "\\*", energy)
  #energy<-strsplit(energy, "*", fixed=T)[[1]][2]
  #energy<-as.numeric(energy)
  return(energy)
}

getEnergy<-function(duplex){
  as.numeric(gsub("\\)","",substr(duplex,regexpr("-",duplex),nchar(duplex))))
}

get_ref_seq_start <- function(duplex){
  start <- substr(duplex,1,regexpr(":",duplex)) 
  start <- gsub("\\(","",start) 
  start <- gsub("\\)","",start)
  start <- gsub("&","",start)
  start <- gsub("\\.","",start)
  start <- gsub(" ","",start)
  start <- substr(start,1,regexpr(",",start)-1)
  start <- as.numeric(start)
  return(start)
}

get_ref_seq_end <- function(duplex){
  end <- substr(duplex,1,regexpr(":",duplex)) 
  end <- gsub("\\(","",end) 
  end <- gsub("\\)","",end)
  end <- gsub("&","",end)
  end <- gsub("\\.","",end)
  end <- gsub(" ","",end)
  end <- gsub(":","",end)
  end <- substr(end,regexpr(",",end)+1,nchar(end))
  end <- as.numeric(end)
  return(end)
}

get_footprint_seq <- function(duplex,target_seq){
  start <- substr(duplex,regexpr(":",duplex)+1,nchar(duplex)) 
  start <- substr(start,1,regexpr("\\(",start)-1) 
  start <- gsub(" ","",start)
  end <- as.numeric(substr(start,regexpr(",",start)+1,nchar(start))) 
  start <- as.numeric(substr(start,1,regexpr(",",start)-1))
  
  footprint <- as.character(target_seq)
  footprint <- substr(footprint,start,end)
  return(footprint)
}

get_miRNA_structure <- function(duplex,ref_seq){
  start <- get_ref_seq_start(duplex) 
  miR_str <- substr(duplex,1,regexpr("&",duplex)-1)
  delta <- start-1
  delta <- strrep(".", delta)
  miR_str <- paste0(delta,miR_str) 
  delta <- nchar(ref_seq)-nchar(miR_str)
  delta <- strrep(".", delta)
  miR_str <- paste0(miR_str,delta) 
  return(miR_str)
}
