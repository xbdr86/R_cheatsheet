getDuplexEnergy<-function(ref_seq,target_seq){
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
