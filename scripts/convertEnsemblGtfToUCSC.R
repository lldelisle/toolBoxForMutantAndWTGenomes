options(stringsAsFactors = F, scipen = 999)

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

if (length(commandArgs(TRUE)) == 0) {
    cat("Choose the alias file.\n")
    alias.file <- file.choose()
    cat("Choose the ensembl gtf.\n")
    original.gtf <- file.choose()
    output.gtf <- gsub(".gtf$", "_UCSC.gtf", gsub(".gz$", "", original.gtf))
} else {
    if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help") {
        cat("Usage: Rscript convertEnsemblGtfToUCSC.R pathForAliasFile pathForInputGtf pathForOutput\n")
        stop()
    }
    alias.file <- commandArgs(TRUE)[1]
    original.gtf <- commandArgs(TRUE)[2]
    output.gtf <- commandArgs(TRUE)[3]
}


##################################################
cat("Loading gtf file...")
input.gr <- readGFFAsGRanges(original.gtf)
cat("Done.\n")
###################################################
alias.df <- read.delim(alias.file)
new.names <- alias.df$X..ucsc
names(new.names) <- alias.df$ensembl
output.gr <- renameSeqlevels(input.gr, new.names[seqlevels(input.gr)])

cat("Writting new gtf file to ", output.gtf, "...")
export.gff(output.gr, output.gtf)
cat("Done.\n")
