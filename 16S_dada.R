library(dada2)
library(Biostrings)
library(ShortRead)

setwd("/home/cmacfarland/ProjTest1/16S_Parasite_Microbiome/dada_algorithm/FASTQ_Generation")
fasqraw = getwd()

fnFs = sort(list.files(fasqraw, pattern = "R1_001.fastq.gz", full.names = T))
fnRs = sort(list.files(fasqraw, pattern = "R2_001.fastq.gz", full.names = T))

#cut primers
FWD = "GTGYCAGCMGCCGCGGTA"
REV = "GGACTACHVGGGTWTCTAAT"

fnFs.filtN = file.path(fasqraw, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN = file.path(fasqraw, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 2, multithread = 8) 

cutadapt = "/home/cmacfarland/miniconda3/bin/cutadapt" #cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

path.cut = file.path(fasqraw, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut = file.path(path.cut, basename(fnFs))
fnRs.cut = file.path(path.cut, basename(fnRs))

FWD.RC = dada2:::rc(FWD)
REV.RC = dada2:::rc(REV)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-m", 150, # attempt to filter reads that are shorter than 150bp (nothing should be shorter than ~450)
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files (the files with Ns removed)
}

cutFs = sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs = sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))

get.sample.name = function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names = unname(sapply(cutFs, get.sample.name))

#filter and trim
filtFs = file.path(path.cut, "filtered", basename(cutFs))
filtRs = file.path(path.cut, "filtered", basename(cutRs))

out = filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(220,190),
                    maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=8) 

errF = learnErrors(filtFs, multithread=8)
errR = learnErrors(filtRs, multithread=8)

#Dereplication
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) = sample.names
names(derepRs) = sample.names

#apply core sample inference algorithm
dadaFs = dada(derepFs, err=errF, multithread=8, pool = "pseudo")
dadaRs = dada(derepRs, err=errR, multithread=8, pool = "pseudo")

#merge
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#seqtab
seqtab = makeSequenceTable(mergers)
seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% seq(252,295)]

#remove chimera
seqtab.nochim = removeBimeraDenovo(seqtab2, method="consensus", multithread=8, verbose=TRUE)

#assign taxonomy
dadataxa16s = assignTaxonomy(seqtab.nochim, "/home/cmacfarland/ProjTest1/16S_Parasite_Microbiome/silva_nr_v132_train_set.fa.gz", tryRC = T, multithread = 8)
#assign species 
dadataxasp16s = addSpecies(dadataxa16s, "/home/cmacfarland/ProjTest1/16S_Parasite_Microbiome/silva_species_assignment_v132.fa.gz")

taxa.print = dadataxasp16s # Removing sequence rownames for display only
rownames(taxa.print) = NULL

# giving ourselves text files to pull off of synergy, can load back in to R on local machine
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "dadaspASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "dadaspASVs_counts.txt", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- dadataxasp16s
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "dadaspASVs_taxonomy.txt", sep="\t", quote=F, col.names=NA)