#!/usr/bin/env Rscript

library("dada2"); packageVersion("dada2")
library("yaml"); packageVersion("yaml")
#library("phyloseq"); packageVersion("phyloseq")

args <- commandArgs(trailingOnly = TRUE)
# set output paths
output_path  <- args[1]  #"/data/ampliseq/201014_16Sv4_Test2/DADA2_summary" #THIS SHOULD BE THE DADA2_summary folder
#output_path2 <- "/data/ampliseq/201014_16Sv4_Test2/DADA2_outputs" #THIS SHOULD BE THE DADA2_output folder
# create output directories
mainDir <- output_path
#dir.create( mainDir )
subDir <- "DADA2_summary"
dir.create(file.path(mainDir, subDir))
#subDir <- "DADA2_outputs"
#dir.create(file.path(mainDir, subDir))

# list samples
# the directory containing the fastq files
path <- args[2] # "/data/myBaseSpace/201014_16Sv4_Test2/fastq"
list.files(path)

# read a list of parameters from a yaml file
config_path  <- args[3]
config <- yaml.load_file(config_path, as.named.list=TRUE)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME-XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_L001_R1_001.fastq.gz"), `[`, 1)

#Inspect read quality profiles of F and R
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
fns2=cbind(fnFs[1:10], fnRs[1:10])
plotQualityProfile(fns2)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncLen=c(250,180), maxEE=c(0.3,0.3),truncQ=2       c(240,160), trimLeft=31,
# v1v3 -> c(290,290),
# v4 -> c(220,160), maxEE=c(2,2),
# v3v4 -> c(250,230) maxEE=c(3,6)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
	      truncLen=c(config$truncLenR1,config$truncLenR2),
        maxN=0,
	      trimLeft=c(config$trimLeftR1,config$trimLeftR2),
	      maxEE=c(config$maxEER1,config$maxEER2),
	      #minQ=0, don't use
	      truncQ=config$truncQ,
	      rm.phix=TRUE,
        compress=TRUE,
	      verbose=TRUE,
	      multithread=TRUE)
print(out)
# filter out the sample stat matrix reads.out ==0
out <- out[out[, "reads.out"] > 0,]

#write.table(out, "Sample_stat.tsv", sep="\t", quote=F, col.names=NA)

# filter out the empty samples
filtFs <- filtFs[file.exists(filtFs) ]
filtRs <- filtRs[file.exists(filtRs) ]

# get the filtered sample names
sample.names <- sapply(strsplit(basename(filtFs), "_F_filt.fastq.gz"), `[`, 1)
names(filtFs) <- sample.names
names(filtRs) <- sample.names



# generate error model
errF <- learnErrors(filtFs, randomize=TRUE, multithread=TRUE)
errR <- learnErrors(filtRs, randomize=TRUE, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names


# Sample Inference    ="pseudo" POOL=FALSE
#dadaFs <- dada(derep_forward, err=errF, pool=TRUE, multithread=TRUE)
#dadaRs <- dada(derep_reverse, err=errR, pool=TRUE, multithread=TRUE)

# POOL=FALSE
dadaFs <- dada(derep_forward, err=errF, pool=FALSE, multithread=TRUE)
dadaRs <- dada(derep_reverse, err=errR, pool=FALSE, multithread=TRUE)

#selfConsist=TRUE, infer both sample decomposition and its parameters of its parameters from the data
#dadaFs <- dada(derep_forward, NULL, selfConsist=TRUE, pool=TRUE, multithread=TRUE)
#dadaRs <- dada(derep_reverse, NULL, selfConsist=TRUE, pool=TRUE, multithread=TRUE)


dadaFs[[1]]

# merge PE reads
mergers <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, verbose=TRUE, trimOverhang=TRUE,
                       minOverlap=config$minOverlap )  # default =12
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
head (table(nchar(getSequences(seqtab))))


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

if (TRUE) {
#filter ASVs by length for V3V4
#MINLEN <- 350
#MAXLEN <- 410

#filter ASVs by length for V4
MINLEN <- config$ASVMINLEN #200
MAXLEN <- config$ASVMAXLEN #250

seqlens <- nchar(getSequences(seqtab.nochim))
seqtab.nochim.filt <- seqtab.nochim [,seqlens >= MINLEN & seqlens <= MAXLEN]

#filter ASVs by abundance
MINABUND <- 10
abundances <- colSums(seqtab.nochim.filt)
seqtab.nochim <- seqtab.nochim.filt[,abundances >= MINABUND]
}

sum(seqtab.nochim)/sum(seqtab)

# Track the number of reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track, file.path(output_path,"Sample_stats.tsv"), sep="\t", quote=F, col.names=NA)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, args[4],
		       minBoot=config$minBoot, # default 50
		       multithread=TRUE,
		       tryRC=TRUE)

# no species level any more
#taxa <- addSpecies(taxa, "/mnt/efs/scratch/Xmeng/data/MiSeq/silvaDB/silva_species_assignment_v138.fa.gz", tryRC=TRUE, allowMultiple=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

 # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

  # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output_path, "ASVs.fa" )  )

  # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(output_path, "ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)

  # tax table:
  # creating table of taxonomy and setting any that are unclassified as "NA"
rownames(taxa) <- gsub(pattern=">", replacement="", x=asv_headers)
write.table(taxa, file.path(output_path, "ASVs_taxonomy_dada2.tsv"), sep = "\t", quote=F, col.names=NA)


###############################################################################
############## Output summary at each level of taxa    ########################
if (FALSE) {
  taxa <- assignTaxonomy(seqtab.nochim, "/mnt/efs/scratch/Xmeng/data/MiSeq/silvaDB/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, tryRC=TRUE)
  tax_tab <- tax_table(taxa)
  pseq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), tax_tab)


  otut<-otu_table(tax_glom(pseq,taxrank = "Phylum"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Phylum"))
  write.table(taxt,  file.path(output_path,"Phylum_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path,"Phylum_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  tax_ph<- paste(a1,a2,sep=";")
  write.table(otut, file.path(output_path,"Phylum_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)

  #otut2 <- cbind(otut, purity=round(rowMaxs(otut)/rowSums(otut)*100, 1))
  #write.table(otut2, file.path(output_path,"Phylum_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)

  otut<-otu_table(tax_glom(pseq,taxrank = "Class"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Class"))
  write.table(taxt,  file.path(output_path,"Class_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path,"Class_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  tax_ph<- paste(a1,a2,a3,sep=";")
  write.table(otut, file.path(output_path,"Class_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Order"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Order"))
  write.table(taxt,  file.path(output_path,"Order_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path,"Order_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  tax_ph<- paste(a1,a2,a3,a4,sep=";")
  write.table(otut, file.path(output_path,"Order_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Family"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Family"))
  write.table(taxt,  file.path(output_path,"Family_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path,"Family_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  a5<-as.character(a[,5])
  tax_ph<- paste(a1,a2,a3,a4,a5,sep=";")
  write.table(otut, file.path(output_path,"Family_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)


  otut<-otu_table(tax_glom(pseq,taxrank = "Genus"))
  taxt<-tax_table(tax_glom(pseq,taxrank = "Genus"))
  write.table(taxt,  file.path(output_path,"Genus_taxa.txt"), sep='\t', row.names=TRUE, col.names=TRUE)
  a <- read.table(file.path(output_path,"Genus_taxa.txt"))

  a1<-as.character(a[,1])
  a2<-as.character(a[,2])
  a3<-as.character(a[,3])
  a4<-as.character(a[,4])
  a5<-as.character(a[,5])
  a6<-as.character(a[,6])
  tax_ph<- paste(a1,a2,a3,a4,a5,a6,sep=";")
  write.table(otut, file.path(output_path,"Genus_summary.txt"), sep='\t', quote=F, row.names=TRUE, col.names=tax_ph)

}
  if (FALSE) {
      library(DECIPHER); packageVersion("DECIPHER")

      dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
      load("/data/ampliseq/silvaDB/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
      ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=FALSE) # use all processors

      ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
      asv_tax <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
      }))
      colnames(asv_tax) <- ranks
      rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

      #colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

      write.table(asv_tax, "ASVs_taxonomy_decipher.tsv", sep = "\t", quote=F, col.names=NA)


      #Removing likely contaminants
#      library("phyloseq"); packageVersion("phyloseq")
#      library(decontam); packageVersion("decontam")


#      ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
#                     tax_table(taxa))
#      ps
  }
#}
