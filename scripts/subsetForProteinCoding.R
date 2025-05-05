options(stringsAsFactors = FALSE, scipen = 999)

if (!"rtracklayer" %in% installed.packages()) {
  if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
  }
  devtools::install_github("lldelisle/usefulLDfunctions", upgrade = "never")
  library(usefulLDfunctions)
  safelyLoadAPackageInCRANorBioconductor("rtracklayer")
} else {
  library(rtracklayer)
}

rm(list = ls())

library(tools)

if (length(commandArgs(TRUE)) == 0) {
  cat("Choose the gtf file.\n")
  gtfFile <- file.choose()
  cat("Choose the FPKM/count table to subset.\n")
  tableFile <- file.choose()
  cat("Put the name of the column in the FPKM/count table which matches the gene_id:\n")
  geneIdCol <- readLines(con = stdin(), n = 1)
} else {
  if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help" || length(commandArgs(TRUE)) <
    3) {
    cat("Usage: Rscript subsetForProteinCoding.R pathForGtf pathForTable geneIdCol\n")
    stop()
  }
  gtfFile <- commandArgs(TRUE)[1]
  tableFile <- commandArgs(TRUE)[2]
  geneIdCol <- commandArgs(TRUE)[3]
}

##################################################
cat("Loading table file...")
tableToSubset <- read.delim(tableFile, check.names = FALSE)
cat("Done.\n")
if (!geneIdCol %in% colnames(tableToSubset)) {
  stop(paste("The column", geneIdCol, "is not part of column names."))
}
cat("Loading gtf file...")
gtfInput <- readGFF(gtfFile)
cat("Done.\n")
##################################################
if ("gene_biotype" %in% colnames(gtfInput)) {
  subsetGtf <- subset(gtfInput, gene_biotype %in% "protein_coding")
} else {
  subsetGtf <- subset(gtfInput, gene_type %in% "protein_coding")
}
if (nrow(subsetGtf) == 0) {
  stop("There is no gene with gene_biotype protein_coding.")
}
subsettedTable <- tableToSubset[tableToSubset[, geneIdCol] %in% subsetGtf$gene_id, ]

cat("Writting subsetted table...")
write.table(
  subsettedTable,
  paste0(file_path_sans_ext(tableFile), "_onlypc.txt"),
  quote = FALSE, row.names = FALSE, sep = "\t"
)
cat("Done.\n")
