# Documentation

## Inputs
Usually you need fasta file which you can download on [UCSC for example](https://hgdownload.soe.ucsc.edu/downloads.html).
And a gtf, you can get on on [GENCODE](https://www.gencodegenes.org/) but I personally prefer the one from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) which are more exhaustive but sometimes too exhaustive so I filter them.

## Extract a small sequence from a big fasta
- Name of the script: `extractFromFasta.R`
- What is needed:
  - A fasta file from which you want to extract a small sequence (it can be gzipped).
  - The start and the end position you want to extract (1-based=the first position is called 1)
- What it will do:
  - It will create a new fasta file next to the original one containing only the sequence between the start and the end position you gave.
  - The name of the file will be the name of the original fasta_startPos-EndPos.fa
  - The name of the sequence (in the fasta file) will be the name of the original sequence : startPos-EndPos.
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the original fasta file.
    - In the console, the script ask you to write the start position and the end position.

## Extract ensembl id list from gtf and chromosome name
- Name of the script: `getGeneListFromChrAndGTF.R`
- What is needed:
  - A gtf file (it can be gzipped).
- What it will do:
  - It will create a file named genesInchrNfromNAMEOFGTF.txt in which you will have all ensembl id present in the chromosome you asked (one per line).
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the gtf file.
    - Then you need to put, in the console of RStudio, the name of the chromosome or chromosomes separated by comma for example chrM for the mitochondria and chrX,chrY,chrM for sexual chromosomes + mitochondria.
    - Then you need to put the folder where you want to have the output. If you have a mac you can copy the folder from finder and paste it in the console of RStudio.

## Create a refgene for IGV
- Why?
  - IGV is able to read gtf but this is really really slow and sometimes IGV crashes.
  - I discovered that there is a format which IGV is using which is much more faster. This is RefGene.txt. And it allows to see the name of the genes!
  - To know it is a RefGene file, IGV needs that the file name is exactly RefGene.txt.
  - Here is a script which can convert a gtf to RefGene.txt with transcripts information (for the moment which part of the transcript is coding is not put in the RefGene.txt).
- Name of the script: `convertGtfToRefGeneForIGV.R`
- What is needed:
  - A gtf file (it can be gzipped).
- What it will do:
  - It will create a file named RefGene.txt which correspond to every transcript which have exon in the original gtf.
  - For the moment, the script is not trying to find where the transcript is coding. It is like every transcript was a non-coding transcript.
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the gtf file.

## Filter a gtf like Anouk
- Why?
  - Anouk has much more experience in RNA-seq than I have.
  - She explained me how readthrough transcripts (transcripts which overlaps at least 2 genes which are not themselves readthrough transcripts) can biais gene quantification.
  - For example, there is a gene which overlap Hoxd3 and Hoxd4 (Gm28230). If you keep this gene in your gtf, STAR and HTSeq-counts will set as "ambiguous" a lot of reads and thus artificially decrease the level of expression of Hoxd3 and Hoxd4.
  - The ensembl annotations are more complete but may lead to biais in quantifications.
  - To solve this issue, she is filtering Ensembl gtf and only use the FilteredTranscript gtf for analyses.
  - Here is a script to filter your gtf like she would do.
- Name of the script: `filterGtfLikeAnouk.R`
- What is needed:
  - A gtf file (it can be gzipped).
  - Optionnal: a list of monoexon biotypes which would not be count for the readthrough filter.
  - The maximum length of transcript to be keep (default is 2.5Mb).
  - If you want to keep only exons (these are the only features necessary for any downstream analysis).
  - If you want to use the UCSC format.
- What it will do:
  - It will create a file named FilteredTranscriptsOf....gtf which will contains only the transcripts which passed the different filters.
  - Filter 1: Readthrough
    - For each transcript, it will count how many genes are overlapped (at the exon levels) which are not monoexonic and which are not overlapping more than 2 genes.
    - All transcripts which overlaps more than 2 genes of this type will be removed.
  - Filter 2: Long transcripts
    - For each transcript, it will measure the length and if it is more than the maximum length (default is 2.5Mb) it will be removed.
  - Filter 3: Non-coding transcripts in a protein-coding gene.
    - All transcripts which are not protein-coding whereas they are part of a protein-coding gene are removed.
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the gtf file.
    - Then in the console, you can put the list of monoexonic biotypes (I recommand you to type enter, this will take the default ones).
    - In the console, you need to put the maximum length of the transcripts or put enter to choose 2.5Mb.
    - In the console, you need to type Y if you want to use UCSC format (use chr1 instead of 1 and keep only the chromosomes with numbers + X, Y, M)
    - In the console you need to type Y if you want to keep only exons  which I recommand (for the moment, I did not find any software which is using the other features which are genes, transcripts, startCodon...).

## Merge genes with same name which overlap
- Why?
  - I noticed that since Ensembl version 92, Hoxd4 has 2 ensembl_id. Then I realized, it was not the only one to have been split into 2.
  - If you keep the two genes in your gtf, STAR and HTSeq-counts will set as "ambiguous" a lot of reads and thus artificially decrease the level of expression Hoxd4, even if you sum the number of reads asigned to each of the ensembl_id.
  - I made some tests and using the merged version of a gtf compared to a unmerged version:
    - Do not change the mapping result with STAR, so the coverage is the same.
    - Do not change the FPKM for isoforms.
    - Just do the sum of FPKM for both ensembl_id.
    - Increase the number of reads falling into this gene because it decreases the number of ambiguous reads.
  - That's why I recommand to use a merge version of a filtered gtf.
- Name of the script: `mergeGenes_overlap.R`
- What is needed:
  - A gtf file (it can be gzipped).
- What it will do:
  - It will look for ensembl genes id which have the same gene_name and overlaps.
  - It will modify the information for the gene which has the smaller number of exons and put the information of the gene with the higher number of exons.
  - The result will be a gtf with only one ensembl_id for genes which have the same name and overlap.
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the gtf file.

## Subset count or FPKM tables to only keep protein coding genes
- Why?
  - Because sometimes you prefer to stick to protein coding genes for downstream analysis (PCA/Clustering/DE).
- Name of the script: `subsetForProteinCoding.R`
- What is needed:
  - A gtf file (it can be gzipped).
  - A count/FPKM table.
  - The name of the column with gene ids.
- What it will do:
  - It will look for genes which are protein_coding.
  - It create a new table where any line that does not correspond to these genes is discarded.
  - The new table will be saved in the same directory as the input table with the suffix `_onlypc.txt`.
- How to launch it?
  - SOURCE ! (open it in RStudio and click on the source button which is top right).
    - A window will automatically open and you need to choose the gtf file.
    - Then, a window will automatically open and you need to choose the count/FPKM table file.
    - Then, in the console you will need to type the column name in your count/FPKM table that matches gene_id in the gtf.
