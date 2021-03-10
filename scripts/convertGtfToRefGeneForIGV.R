options(stringsAsFactors=F)
if(!"rtracklayer"%in%installed.packages()){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rtracklayer")
}
library(rtracklayer)
rm(list=ls())

if(length(commandArgs(TRUE))==0){
  connection<-stdin()
  cat("Choose the gtf file.\n")
  gtfFile<-file.choose()
} else{
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript dateOfScript_checkGtf.R pathForGtf \n")
    stop()
  }
  connection<-"stdin"
  gtfFile<-commandArgs(TRUE)[1]
}

##################################################
cat("Loading gtf file...")
gtfInput<-readGFF(gtfFile)
cat("Done.\n")
###################################################
#The format of RefGene.txt is
#Col1 number (?)
#Col2 ID
#Col3 chr
#Col4 strand
#Col5 start 0-based included
#Col6 end 0-based excluded
#Col7 start of coding (or end) 0-based included
#Col8 end of coding (or end) 0-based excluded
#Col9 nb of exons
#Col10 begin of exons separated by comma, finished by comma (not necessary)
#Col11 end of exons idem
#Col12 0
#Col13 Name to display
#Col14 cmpl, incmpl or unk
#Col15 cmpl, incmpl or unk
#Col16 phase? for each exon 0,1,2 or -1
gtfInput<-subset(gtfInput,type=="exon")
gtfInput<-gtfInput[order(gtfInput$start,gtfInput$end),]
cat("Building the transcripts from the exons.\n")
transcripts<-unique(gtfInput[,c("seqid","strand","gene_name","transcript_id")])
if(sum(duplicated(transcripts$transcript_id))>0){
  print(transcripts[duplicated(transcripts$transcript_id),])
  stop("The gtf file contains different information for the same transcript.")
}
start<-aggregate(gtfInput$start-1,by=list(transcript_id=gtfInput$transcript_id),FUN=min)
transcripts$start<-start[match(transcripts$transcript_id,start$transcript_id),"x"]
end<-aggregate(gtfInput$end,by=list(transcript_id=gtfInput$transcript_id),FUN=max)
transcripts$end<-end[match(transcripts$transcript_id,end$transcript_id),"x"]
nbExons<-aggregate(gtfInput$exon_id,by=list(transcript_id=gtfInput$transcript_id),FUN=length)
transcripts$nbExons<-nbExons[match(transcripts$transcript_id,nbExons$transcript_id),"x"]
startsExons<-aggregate(gtfInput$start-1,by=list(transcript_id=gtfInput$transcript_id),FUN=paste,collapse=",")
transcripts$startsExons<-startsExons[match(transcripts$transcript_id,startsExons$transcript_id),"x"]
endsExons<-aggregate(gtfInput$end,by=list(transcript_id=gtfInput$transcript_id),FUN=paste,collapse=",")
transcripts$endsExons<-endsExons[match(transcripts$transcript_id,endsExons$transcript_id),"x"]
transcripts$endOriented<-transcripts$end
transcripts$endOriented[transcripts$strand=="-"]<-transcripts$start[transcripts$strand=="-"]
transcripts$phaseExons<-sapply(transcripts$nbExons,function(n){paste(rep("-1",n),collapse=",")})
finalDataFrame<-data.frame(nb=1:nrow(transcripts),id=transcripts$transcript_id,
                           chr=transcripts$seqid,strand=transcripts$strand,
                           start=transcripts$start,end=transcripts$end,
                           scod=transcripts$endOriented,ecod=transcripts$endOriented,
                           nbE=transcripts$nbExons,bE=transcripts$startsExons,
                           eE=transcripts$endsExons,v12=rep(0,nrow(transcripts)),
                           name=transcripts$gene_name,v14=rep("unk",nrow(transcripts)),
                           v15=rep("unk",nrow(transcripts)),v16=rep("unk",nrow(transcripts)))
write.table(finalDataFrame,paste0(dirname(gtfFile),"/RefGene.txt"),sep="\t",quote=F,row.names=F,col.names=F)
