options(scipen=999)
if(!"seqinr"%in%installed.packages()){
  install.packages("seqinr")
}
library(seqinr)
library(tools)
if(length(commandArgs(TRUE))==0){
  cat("Choose the fasta file.\n")
  pathForFasta<-file.choose()
  cat("Write the start position included.\n")
  startPos<-as.numeric(readLines(con=stdin(),n=1))
  cat("Write the end position included.\n")
  endPos<-as.numeric(readLines(con=stdin(),n=1))
  outputPath<-paste0(dirname(pathForFasta),"/",basename(pathForFasta),"_",startPos,"-",endPos,".fa")
} else{
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript dateOfScript_extractFromFasta.R pathForFasta startPos endPos [outputPath]\n")
    stop()
  }
  pathForFasta<-commandArgs(TRUE)[1]
  startPos<-as.numeric(commandArgs(TRUE)[2])
  endPos<-as.numeric(commandArgs(TRUE)[3])
  if(length(commandArgs(TRUE))>3){
    outputPath<-commandArgs(TRUE)[4]
  } else {
    outputPath<-paste0(dirname(pathForFasta),"/",basename(pathForFasta),"_",startPos,"-",endPos,".fa")
  }
}

if(startPos>endPos){
  stop("Start is greater than end. It is not possible to extract.")
}
cat("Loading fasta...\n")
if (file_ext(pathForFasta)=="gz"){
  chrSeq<-read.fasta(gzfile(pathForFasta), seqtype="DNA", forceDNAtolower=FALSE)[[1]]  
} else{
  chrSeq<-read.fasta(pathForFasta, seqtype="DNA", forceDNAtolower=FALSE)[[1]] 
}
cat("Done.\n")
cat("Writing new fasta...")
write.fasta(chrSeq[startPos:endPos],name=paste0(getName(chrSeq),":",startPos,"-",endPos),outputPath,open='w')
cat("Done.\n")
