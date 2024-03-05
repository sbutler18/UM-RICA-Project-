#installing packages
install.packages("BiocManager")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("DECIPHER")
#installing remotesa
install.packages("remotes")
remotes::install_github("omicsCore/SEQprocess")
BiocManager::install("decontam")

#load the libraries 
library(dada2)
library(remotes)
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library("DECIPHER")
library(decontam)


#set your path
path = "/Users/sterling.butler/Downloads/240209_M02476_0584_000000000-LCRB7/Alignment_1/20240211_215111/Fastq"

# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))


#what primeres were used in your study (if you can't find the study, just use 515F and 806R 
#(FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT)

FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...

# to ensure we have the right primers, and the correct orientation of the primers on the reads, we will verify the presence and orientation of these primers in the data.
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences 
#difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other 
#filtering.

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory


#Creating output directory: /Volumes/ROCKET-nano/TMm/PRJNA413177//filtN
#this is going to take awhile
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

# Counts number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))


#direct them to cutadapt
cutadapt <- "/Users/sterling.butler/miniconda3/envs/cutadaptenv/bin/cutadapt"
system2(cutadapt, args = "--version") # Run shell commands from R


#now time to remove the primers; oof this goes on and on
path.cut <- file.path(path, "cutadapt") #adds a directory for the trimmed reads
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#gets the reverse comps
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 


# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags,  "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o",  fnFs.cut[i], # output files
                             fnFs[i])) # input files
}


#let's check if it works 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))


#mine still had some issues, but I honestly don't know what to do about it rather than just note it
#Forward Complement Reverse RevComp
#FWD.ForwardReads       2          0       0       8
#FWD.ReverseReads      53          0       0       0
#REV.ForwardReads       0          0       0       0
#REV.ReverseReads       0          0       0       0


#they are now trimmed and in a folder called cutadapt
# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001.fastq.gz", full.names = TRUE))


# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#get a plot of the quality of the forward reads 
plotQualityProfile(cutFs[7:8])

#Now for the final filtering stage of the process
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, truncLen = 250, 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  
head(out)

#learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE)

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

#construct a squence table
seqtab <- makeSequenceTable(dadaFs)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#inspect distribution sequence lengths
table(nchar(getSequences(seqtab.nochim)))

#########finally get to assign taxonomy with GSR Data Base
unite.ref = "/Users/sterling.butler/Downloads/GSR-DB_V3-V5.fasta"

taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

getwd()
setwd("/Users/sterling.butler/Desktop")

#alternative way to assign taxa with silva database
## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")


## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

## save workspace
save.image(file = "data_tutorial.RData")

#extracting standard goods from dada2
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(taxa.print, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

