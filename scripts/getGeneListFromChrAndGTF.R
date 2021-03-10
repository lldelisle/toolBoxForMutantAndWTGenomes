options(stringsAsFactors=F)
if(!"rtracklayer"%in%installed.packages()){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rtracklayer")
}
library(rtracklayer)
library(tools)
rm(list=ls())

if(length(commandArgs(TRUE))==0){
  cat("Choose the gtf file.\n")
  gtfFile<-file.choose()
  cat("Put the name of the chromosome or the chromosomes separated by comma (for example chrM or chrX,chrY,chrM):\n")
  chrList<-readLines(con=stdin(),n=1)
  cat("Put the path of the folder where you want to put the file:\n")
  outputFolder<-readLines(con=stdin(),n=1)
} else{
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help" || length(commandArgs(TRUE))<3){
    cat("Usage: Rscript dateOfScript_getGeneListFromChrAndGTF.R.R pathForGtf listOfChr outputFolder\n")
    stop()
  }
  gtfFile<-commandArgs(TRUE)[1]
  chrList<-commandArgs(TRUE)[2]
  outputFolder<-commandArgs(TRUE)[3]
}

chrToGet<-strsplit(chrList,",")[[1]]

##################################################
cat("Loading gtf file...")
gtfInput<-readGFF(gtfFile)
cat("Done.\n")
##################################################
subsetGtf<-subset(gtfInput,seqid%in%chrToGet)
if(nrow(subsetGtf)==0){
  cat("There is no gene in these chromosomes.\n")
} else {
  dir.create(outputFolder,recursive=T,showWarnings=F)
  cat(unique(subsetGtf$gene_id),sep="\n",file=paste0(outputFolder,"/genesIn",chrList,"from",basename(file_path_sans_ext(gtfFile)),".txt"))
}
