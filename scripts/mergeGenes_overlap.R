if(!"rtracklayer"%in%installed.packages()){
  source("https://bioconductor.org/biocLite.R")
  biocLite("rtracklayer")
}
library(rtracklayer)
################################################################################
UIinput <- function(s){
  
  #Ask for user input
  x <- readLines(con = connection,1)
  
  #Check if empty
  if (x==""){
    x<-s
  }
  #Return
  return(x)
}
##################################################################
if(length(commandArgs(TRUE))>0){
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript mergeGenes_overlap.R pathForGtf\n")
    stop()
  }
  connection<-"stdin"
  gtfFile<-commandArgs(TRUE)[1]
} else {
  connection<-stdin()
  cat("Choose the gtf file to filter.\n")
  gtfFile<-file.choose()
}
##################################################
cat("Loading gtf file...")
gtfInput<-readGFF(gtfFile)
cat("Done.\n")
simpleGTF<-unique(gtfInput[!is.na(gtfInput$gene_name),c("gene_id","gene_name")])
duplicatedGeneNames<-unique(with(simpleGTF,gene_name[duplicated(gene_name)]))
newGTF<-gtfInput
change<-data.frame(gene_name=character(),before=character(),after=character())
colsToChange<-grep("gene",colnames(gtfInput))
for(gn in duplicatedGeneNames){
  # cat(gn,"\n")
  while(T){
    smallgrl<-reduce(split(GRanges(newGTF[newGTF$type=="exon"&newGTF$gene_name==gn,]),newGTF$gene_id[newGTF$type=="exon"&newGTF$gene_name==gn]))
    smallTable<-lapply(smallgrl,function(g){table(names(smallgrl)[as.matrix(findOverlaps(g,smallgrl))[,2]])})
    # print(smallTable)
    overlapsExons<-unlist(unname(smallTable),use.names = T)
    ensIDwithOverlap<-names(overlapsExons)[duplicated(names(overlapsExons))]
    if(length(ensIDwithOverlap)==0){
      break
    }
    masterEns<-names(which.max(overlapsExons[names(overlapsExons)%in%ensIDwithOverlap]))
    ensIDToRename<-setdiff(names(which(sapply(smallTable,function(t){masterEns%in%names(t)}))),masterEns)
    if(length(ensIDToRename)>0){
      change<-rbind(change,data.frame(gene_name=rep(gn,length(ensIDToRename)),before=ensIDToRename,after=rep(masterEns,length(ensIDToRename))))
      masterCols<-unique(newGTF[newGTF$gene_id==masterEns,colsToChange])
      if(nrow(masterCols)>1){
        cat("There are inconstistancies in \n")
        print(masterCols)
        cat("Only the first line will be kept\n")
        masterCols<-materCols[1,]
      }
      newGTF[newGTF$gene_id%in%c(ensIDToRename,masterEns),colsToChange]<-masterCols
    } else {
      break
    }
  }
}
name<-gsub(".gtf","",basename(gtfFile))
name<-gsub(".gz","",name)
write.table(change,paste0(dirname(gtfFile),"/renamedEnsIDOf",name,".txt"),sep="\t",quote=F,row.names=F)
cat("Writing the merged gtf...")
rownames(newGTF)<-NULL
export.gff(newGTF,paste0(dirname(gtfFile),"/mergeOverlapGenesOf",name,".gtf"),format="gtf")
cat("Done.\n")
