
# P. Murat, MRC Laboratory of Molecular Biology, February 2022
# Code for generating Figures 5 and S7 
# in Guilbaud et al. Determination of human DNA replication origin position and efficiency reveals principles of initiation zone organisation

library(dplyr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gplots)
library(ggpubr)
library(caret)
library(factoextra)

#########################################################################################################
#########################################################################################################
# Load origin coordinates, format datasets and define useful functions
#########################################################################################################
#########################################################################################################

###############
# Origin coordinate and efficiencies

ori <- read.table("./Dataset/Final_HL_enrchiment_of_2_for_3h_All_042019.gDNA_rm_dup_single_nt_position_by_100_bp_windows.bed", header = F)
colnames(ori) <- c("chr", "start", "end", "EFF")
nrow(ori) # 23,905 origins
# Prepare EFF 3-quantiles
quantile(ori$EFF, c(0.333,0.666))
#  33.3%  66.6% 
#  0.8418 0.9065 
# Add origin id
ori <- ori %>% mutate(chr = paste("chr", chr, sep = ""), ori.id = paste(chr, ":", start, "-", end, sep = ""))
# Define origin classes
ori.low <- ori %>% filter(EFF < 0.8418)
ori.medium <- ori %>% filter(EFF >= 0.8418 & EFF < 0.9065)
ori.high <- ori %>% filter(EFF >= 0.9065)
write.table(ori.low, "./Dataset/ori.low.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ori.medium, "./Dataset/ori.medium.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ori.high, "./Dataset/ori.high.bed", sep = "\t", quote = F, col.names = F, row.names = F)

###############
# Repeat and DNA structure features

# non B-DNA structure are form Non-B DB v2.0
# https://ncifrederick.cancer.gov/bids/ftp/?nonb
# reference: Cer et al. NAR 2013, 41, D95-D100

# 7 classes of sequences/structures:
# APR : A-phased repeats
# DR : Direct repeats
# GQ : G-quadruplexes
# IR : Inverted repeats
# MR : Mirror repeats
# STR : short tandem repeats
# Z : Z-DNA

# .tsv files of individual chromosomes were downloaded and combined to generate .csv files reporting the genome-wide distribution of each structural feature
# After formatting, the relevant bed files were used to generate bedgraph summaries of feature coverage using the bedtools genomecov function
# Finally bedgraph files were converted into bigwig files

# The following commands show as an example how the data for G-quadruplexes (GQ) were processed
GQ.files <- list.files("./Dataset/nonB-DB_hg38/GQ/") # Folder containing tsv for individual chromosomes
GQ.hg38 <- NULL
for (i in GQ.files) {
  file <- paste("./Dataset/nonB-DB_hg38/GQ/", i, sep = "")
  p <- read.csv(file, sep = "\t")
  GQ.hg38 <- rbind(GQ.hg38, p)
}
write.csv(GQ.hg38, "./Dataset/nonB-DB_hg38/GQ.hg38.csv", row.names = F)
# Prepare GQ data for coverage analysis
GQ.hg38.bed <- GQ.hg38 %>% dplyr::select(Sequence_name, Start, Stop)
GQ.hg38.bed <- data.frame(lapply(GQ.hg38.bed.2, function(x) {gsub("chrMT", "chrM", x)}))
write.table(GQ.hg38.bed, "./Dataset/nonB-DB_hg38/GQ.hg38.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# Use bedtools to prepare bedgraph files (command run from the terminal) - hg38.chrom.sizes available from https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/
bedtools genomecov -bga -i GQ.hg38.bed -g hg38.chrom.sizes > GQ.hg38.bedgraph
# Convert to bigwig file (commands run from the terminal)
sort -k1,1 -k2,2n GQ.hg38.bedgraph > GQ.hg38.sorted.bedgraph
bedGraphToBigWig GQ.hg38.sorted.bedgraph hg38.chrom.sizes GQ.hg38.bw

# A similar pipeline was used to generate MR.hg38.bw, IR.hg38.bw, DR.hg38.bw, STR.hg38.bw and Z.hg38.bw

###############
# Base composition and skews in nucleotides

# Define ori +- 2100 bp domains
ori.2kb <- ori %>% mutate(chr = gsub("chr", "", chr)) %>% mutate(start = as.integer(start - 2100), end = as.integer(end + 2100)) %>% filter(start > 0)
# Save bed
write.table(ori.2kb, "./Dataset/BaseComp/ori.2kb.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# tile origins with 100 nt bins using bedtools makewindows function (run from the terminal)
sort -k1,1 -k2,2n ori.2kb.bed > ori.2kb.sorted.bed
bedtools makewindows -b ori.2kb.sorted.bed -w 100 > ori.2kb.100nt.windows.bed
# Compute base composition using bedtools nuc function and hg38 (run from the terminal)
bedtools nuc -fi ./Genomes/Homo_sapiens.GRCh38.dna.primary_assembly -bed ori.2kb.20nt.windows.bed > ori.2kb.20nt.nuc.bed
# Open result
ori.nuc <- read.csv("Dataset/BaseComp/ori.2kb.100nt.nuc.bed", sep = "\t")
colnames(ori.nuc) <- c("chr", "start", "end", "AT", "GC", "A", "C", "G", "T", "N", "other", "length")
# Compute skews
ori.nuc <- ori.nuc %>% mutate(GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T))
# Prepare bed files
ori.GC.bed <- ori.nuc %>% mutate(chr = paste("chr", chr, sep = ""), id = paste(chr, start, sep = "_")) %>% dplyr::select(chr, start, end, GC, id) %>% distinct(id, .keep_all= TRUE) %>% dplyr::select(-id)
ori.GC.skew.bed <- ori.nuc %>% mutate(chr = paste("chr", chr, sep = ""), id = paste(chr, start, sep = "_")) %>% dplyr::select(chr, start, end, GC.skew, id) %>% distinct(id, .keep_all= TRUE) %>% dplyr::select(-id) %>% drop_na()
ori.AT.skew.bed <- ori.nuc %>% mutate(chr = paste("chr", chr, sep = ""), id = paste(chr, start, sep = "_")) %>% dplyr::select(chr, start, end, AT.skew, id) %>% distinct(id, .keep_all= TRUE) %>% dplyr::select(-id) %>% drop_na()
write.table(ori.GC.bed, "./Dataset/BaseComp/ori.GC.conversion.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ori.GC.skew.bed, "./Dataset/BaseComp/ori.GC.skew.conversion.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ori.AT.skew.bed, "./Dataset/BaseComp/ori.AT.skew.conversion.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# Use bedtools to prepare bigwig files (run from the terminal)
sort -k1,1 -k2,2n ori.GC.conversion.bed > ori.GC.conversion.sorted.bed
sort -k1,1 -k2,2n ori.GC.skew.conversion.bed > ori.GC.skew.conversion.sorted.bed
sort -k1,1 -k2,2n ori.AT.skew.conversion.bed > ori.AT.skew.conversion.sorted.bed
bedGraphToBigWig ori.GC.conversion.sorted.bed hg38.chrom.sizes GC.hg38.bw
bedGraphToBigWig ori.GC.skew.conversion.sorted.bed hg38.chrom.sizes GC.skew.hg38.bw
bedGraphToBigWig ori.AT.skew.conversion.sorted.bed hg38.chrom.sizes AT.skew.hg38.bw

###############
# CpG islands

# Coordinate of CpG islands were recovered from UCSC
# from table browser (http://genome.ucsc.edu/cgi-bin/hgTables) to generate cpgIslandExt.txt
# Then formated as followed
CpG.island <- read.csv("./Dataset/BaseComp/cpgIslandExt.txt", sep = "\t", header = F) %>% dplyr::select(V2, V3, V4)
# Select chromosome
chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
         "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
         "chr21", "chr22", "chrX", "chrY", "chrM")
CpG.island <- CpG.island %>% filter(V2 %in% chr)
write.table(CpG.island, "./Dataset/BaseComp/CpG.island.hg38.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# Convert to bigwig file (run from the terminal)
sort -k1,1 -k2,2n CpG.island.hg38.bed > CpG.island.hg38.sorted.bed
bedtools genomecov -bga -i CpG.island.hg38.sorted.bed -g hg38.chrom.sizes > CpG.island.hg38.bedgraph
sort -k1,1 -k2,2n CpG.island.hg38.bedgraph > CpG.island.hg38.sorted.bedgraph
bedGraphToBigWig CpG.island.hg38.sorted.bedgraph hg38.chrom.sizes CpG.island.hg38.bw

###############
# Defining useful function

# coefficient of determination
r2 <- function(x,y) {
  r2 <- 1 - (sum((x-y )^2)/sum((x-mean(x))^2))
  return(r2)
}

# coverage function
cov <- function(x) { # input is a dataframe with two columns (start and end)
  if(nrow(x) == 0) {
    return()
  } else {
    return(unlist(lapply(1:nrow(x), function(i) c(x$start[i]:x$end[i]))))
  }
}

#########################################################################################################
#########################################################################################################
# Assess local changes in base composition at origins
#########################################################################################################
#########################################################################################################

# Prepare bed files reporting 20 kb domains centered on origins binned in 100 bp windows (200 windows)

# Prepare genomic coordinates of origin centers +- 10kb

ori.low.10kb <- ori.low %>% mutate(center = round((start+end)/2)) %>% mutate(start = as.integer(center-10000), end = as.integer(center+10000)) %>% select(chr, start, end, ori.id, EFF)
ori.medium.10kb <- ori.medium %>% mutate(center = round((start+end)/2)) %>% mutate(start = as.integer(center-10000), end = as.integer(center+10000)) %>% select(chr, start, end, ori.id, EFF)
ori.high.10kb <- ori.high %>% mutate(center = round((start+end)/2)) %>% mutate(start = as.integer(center-10000), end = as.integer(center+10000)) %>% select(chr, start, end, ori.id, EFF)

# Save bed files

write.table(ori.low.10kb, "./BaseComp/ori.low.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(ori.medium.10kb, "./BaseComp/ori.medium.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)
write.table(ori.high.10kb, "./BaseComp/ori.high.10kb.bed", sep="\t", col.names = F, row.names = F, quote = F)

# Use bedools to prepare coordinate of non-overlapping 100 nt windows around origins
# from the terminal
sort -k1,1 -k2,2n ori.low.10kb.bed > ori.low.10kb.sorted.bed
sort -k1,1 -k2,2n ori.medium.10kb.bed > ori.medium.10kb.sorted.bed
sort -k1,1 -k2,2n ori.high.10kb.bed > ori.high.10kb.sorted.bed
bedtools makewindows -b ori.low.10kb.sorted.bed -w 100 -i winnum > ori.low.10kb.100nt.split.bed
bedtools makewindows -b ori.medium.10kb.sorted.bed -w 100 -i winnum > ori.medium.10kb.100nt.split.bed
bedtools makewindows -b ori.high.10kb.sorted.bed -w 100 -i winnum > ori.high.10kb.100nt.split.bed

# Import results as Granges object

IniSeq.low.10kb.gr <- import("./BaseComp/ori.low.10kb.100nt.split.bed")
IniSeq.medium.10kb.gr <- import("./BaseComp/ori.medium.10kb.100nt.split.bed")
IniSeq.high.10kb.gr <- import("./BaseComp/ori.high.10kb.100nt.split.bed")

# Compute base composition

IniSeq.low.views <- Views(Hsapiens, IniSeq.low.10kb.gr)
IniSeq.low.bc <- letterFrequency(IniSeq.low.views, c("A", "C", "G", "T"), as.prob = FALSE)
IniSeq.low.bc.df <- cbind.data.frame(bin = IniSeq.low.10kb.gr$name, IniSeq.low.bc)
IniSeq.low.summary <- IniSeq.low.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

IniSeq.medium.views <- Views(Hsapiens, IniSeq.medium.10kb.gr)
IniSeq.medium.bc <- letterFrequency(IniSeq.medium.views, c("A", "C", "G", "T"), as.prob = FALSE)
IniSeq.medium.bc.df <- cbind.data.frame(bin = IniSeq.medium.10kb.gr$name, IniSeq.medium.bc)
IniSeq.medium.summary <- IniSeq.medium.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

IniSeq.high.views <- Views(Hsapiens, IniSeq.high.10kb.gr)
IniSeq.high.bc <- letterFrequency(IniSeq.high.views, c("A", "C", "G", "T"), as.prob = FALSE)
IniSeq.high.bc.df <- cbind.data.frame(bin = IniSeq.high.10kb.gr$name, IniSeq.high.bc)
IniSeq.high.summary <- IniSeq.high.bc.df %>% group_by(bin) %>% summarise_at(c("A", "C", "G", "T"), sum) %>% mutate(dist = (as.numeric(bin)-100)*100, GC.skew = (G-C)/(G+C), AT.skew = (A-T)/(A+T)) %>% 
  mutate(A.freq = A/(A+C+T+G), C.freq = C/(A+C+T+G), G.freq = G/(A+C+T+G), T.freq = T/(A+C+T+G))

# Plot results for IniSeq origins binned by efficiency

IniSeq.df1 <- IniSeq.low.summary %>% mutate(freq = (G.freq+C.freq)/(A.freq+C.freq+G.freq+T.freq)) %>% select(dist, freq) %>% mutate(class = "low")
IniSeq.df2 <- IniSeq.medium.summary %>% mutate(freq = (G.freq+C.freq)/(A.freq+C.freq+G.freq+T.freq)) %>% select(dist, freq) %>% mutate(class = "medium")
IniSeq.df3 <- IniSeq.high.summary %>% mutate(freq = (G.freq+C.freq)/(A.freq+C.freq+G.freq+T.freq)) %>% select(dist, freq) %>% mutate(class = "high")
IniSeq.plot <- rbind(IniSeq.df1, IniSeq.df2, IniSeq.df3)

GC.plot <- ggplot(IniSeq.plot, aes(x=dist, y=freq, color=class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("GC content") + ggtitle("%GC") +
  theme_bw() + theme(aspect.ratio=1)

IniSeq.df4 <- IniSeq.low.summary %>% select(dist, skew = GC.skew) %>% mutate(class = "low")
IniSeq.df5 <- IniSeq.medium.summary %>% select(dist, skew = GC.skew) %>% mutate(class = "medium")
IniSeq.df6 <- IniSeq.high.summary %>% select(dist, skew = GC.skew) %>% mutate(class = "high")
IniSeq.plot.2 <- rbind(IniSeq.df4, IniSeq.df5, IniSeq.df6)

GC.skew.plot <- ggplot(IniSeq.plot.2, aes(x=dist, y=skew, color=class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("GC nucleotide skew") + ggtitle("GC skew") +
  theme_bw() + theme(aspect.ratio=1)

IniSeq.df7 <- IniSeq.low.summary %>% select(dist, skew = AT.skew) %>% mutate(class = "low")
IniSeq.df8 <- IniSeq.medium.summary %>% select(dist, skew = AT.skew) %>% mutate(class = "medium")
IniSeq.df9 <- IniSeq.high.summary %>% select(dist, skew = AT.skew) %>% mutate(class = "high")
IniSeq.plot.3 <- rbind(IniSeq.df7, IniSeq.df8, IniSeq.df9)

AT.skew.plot <- ggplot(IniSeq.plot.3, aes(x=dist, y=skew, color=class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("AT nucleotide skew") + ggtitle("AT skew") +
  theme_bw() + theme(aspect.ratio=1)

pdf("./Rplot/Figure_4A_4B_S4B.pdf", width=15, height=4, useDingbats=FALSE)
ggarrange(GC.plot, GC.skew.plot, AT.skew.plot, ncol = 3, nrow = 1)
dev.off()

#########################################################################################################
#########################################################################################################
# Assess the distribution of non-B DNA structures and repeats around origins
#########################################################################################################
#########################################################################################################

###############
# Select origin centres

ori.center <- ori %>% mutate(ori.width = end - start, start = round((start+end)/2), end = start+1, id = paste(chr, start, sep = "_")) %>% dplyr::select(-ori.id) %>%  
  mutate(class = case_when(EFF < 0.8418 ~ "Low", 
                           EFF >= 0.8418 & EFF < 0.9065 ~ "Medium",
                           TRUE ~ "High"))
write.table(ori.center, "./Dataset/ori.center.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# From the terminal
sort -k1,1 -k2,2n ori.center.bed > ori.center.sorted.bed

###############
# Compute distances from structure/repeats to origins centres

# Using the bedtools closest function, run in the terminal

sort -k1,1 -k2,2n APR.hg38.bed > APR.hg38.sorted.bed
sort -k1,1 -k2,2n DR.hg38.bed > DR.hg38.sorted.bed
sort -k1,1 -k2,2n GQ.hg38.bed > GQ.hg38.sorted.bed
sort -k1,1 -k2,2n IR.hg38.bed > IR.hg38.sorted.bed
sort -k1,1 -k2,2n MR.hg38.bed > MR.hg38.sorted.bed
sort -k1,1 -k2,2n STR.hg38.bed > STR.hg38.sorted.bed
sort -k1,1 -k2,2n Z.hg38.bed > Z.hg38.sorted.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/APR.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/APR.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/DR.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/DR.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/GQ.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/GQ.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/IR.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/IR.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/MR.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/MR.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/STR.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/STR.ori.closest.bed
bedtools closest -a ./Dataset/nonB-DB_hg38/Z.hg38.sorted.bed -b ./Dataset/ori.center.sorted.bed > ./Struct_analysis/Z.ori.closest.bed

###############
# Analyse the distribution of structures/repeats

# APR
APR.ori.closest <- read.csv("./Struct_analysis/APR.ori.closest.bed", sep = "\t", header = F)
colnames(APR.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "Tracts", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select APR close to origins
APR.ori.closest.2.5kb <- APR.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
APR.ori.low.closest.2.5kb <- APR.ori.closest.2.5kb %>% filter(class == "Low")
APR.ori.low.cov <- cov(APR.ori.low.closest.2.5kb)
APR.ori.low.dens <- as.data.frame(table(APR.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(APR.ori.low.cov)), APR.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, APR.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
APR.ori.medium.closest.2.5kb <- APR.ori.closest.2.5kb %>% filter(class == "Medium")
APR.ori.medium.cov <- cov(APR.ori.medium.closest.2.5kb)
APR.ori.medium.dens <- as.data.frame(table(APR.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(APR.ori.medium.cov)), APR.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, APR.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
APR.ori.high.closest.2.5kb <- APR.ori.closest.2.5kb %>% filter(class == "High")
APR.ori.high.cov <- cov(APR.ori.high.closest.2.5kb)
APR.ori.high.dens <- as.data.frame(table(APR.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(APR.ori.high.cov)), APR.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, APR.dens) %>% 
  mutate(class = "High")
# Bind
APR.ori.df <- rbind(APR.ori.low.dens, APR.ori.medium.dens, APR.ori.high.dens)
# Plot
APR.plot <- APR.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = APR.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("APR Percentage coverage") + ggtitle("A-phased repeats") + ylim(0,1) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# DR
DR.ori.closest <- read.csv("./Struct_analysis/DR.ori.closest.bed", sep = "\t", header = F)
colnames(DR.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "Spacer", "Repeated", "Subset", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select DR close to origins
DR.ori.closest.2.5kb <- DR.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
DR.ori.low.closest.2.5kb <- DR.ori.closest.2.5kb %>% filter(class == "Low")
DR.ori.low.cov <- cov(DR.ori.low.closest.2.5kb)
DR.ori.low.dens <- as.data.frame(table(DR.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(DR.ori.low.cov)), DR.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, DR.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
DR.ori.medium.closest.2.5kb <- DR.ori.closest.2.5kb %>% filter(class == "Medium")
DR.ori.medium.cov <- cov(DR.ori.medium.closest.2.5kb)
DR.ori.medium.dens <- as.data.frame(table(DR.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(DR.ori.medium.cov)), DR.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, DR.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
DR.ori.high.closest.2.5kb <- DR.ori.closest.2.5kb %>% filter(class == "High")
DR.ori.high.cov <- cov(DR.ori.high.closest.2.5kb)
DR.ori.high.dens <- as.data.frame(table(DR.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(DR.ori.high.cov)), DR.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, DR.dens) %>% 
  mutate(class = "High")
# Bind
DR.ori.df <- rbind(DR.ori.low.dens, DR.ori.medium.dens, DR.ori.high.dens)
# Plot
DR.plot <- DR.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = DR.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("DR Percentage coverage") + ggtitle("Direct repeats") + ylim(0,10) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# GQ
GQ.ori.closest <- read.csv("./Struct_analysis/GQ.ori.closest.bed", sep = "\t", header = F)
colnames(GQ.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "nIslands.nRuns.maxGQ", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select GQ close to origins
GQ.ori.closest.2.5kb <- GQ.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
GQ.ori.low.closest.2.5kb <- GQ.ori.closest.2.5kb %>% filter(class == "Low")
GQ.ori.low.cov <- cov(GQ.ori.low.closest.2.5kb)
GQ.ori.low.dens <- as.data.frame(table(GQ.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(GQ.ori.low.cov)), GQ.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, GQ.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
GQ.ori.medium.closest.2.5kb <- GQ.ori.closest.2.5kb %>% filter(class == "Medium")
GQ.ori.medium.cov <- cov(GQ.ori.medium.closest.2.5kb)
GQ.ori.medium.dens <- as.data.frame(table(GQ.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(GQ.ori.medium.cov)), GQ.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, GQ.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
GQ.ori.high.closest.2.5kb <- GQ.ori.closest.2.5kb %>% filter(class == "High")
GQ.ori.high.cov <- cov(GQ.ori.high.closest.2.5kb)
GQ.ori.high.dens <- as.data.frame(table(GQ.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(GQ.ori.high.cov)), GQ.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, GQ.dens) %>% 
  mutate(class = "High")
# Bind
GQ.ori.df <- rbind(GQ.ori.low.dens, GQ.ori.medium.dens, GQ.ori.high.dens)
# Plot
GQ.plot <- GQ.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = GQ.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("GQ Percentage coverage") + ggtitle("G-quadruplexes") + ylim(0,10) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# IR
IR.ori.closest <- read.csv("./Struct_analysis/IR.ori.closest.bed", sep = "\t", header = F)
colnames(IR.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "Spacer", "Permutations", "Subset", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select IR close to origins
IR.ori.closest.2.5kb <- IR.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
IR.ori.low.closest.2.5kb <- IR.ori.closest.2.5kb %>% filter(class == "Low")
IR.ori.low.cov <- cov(IR.ori.low.closest.2.5kb)
IR.ori.low.dens <- as.data.frame(table(IR.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(IR.ori.low.cov)), IR.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, IR.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
IR.ori.medium.closest.2.5kb <- IR.ori.closest.2.5kb %>% filter(class == "Medium")
IR.ori.medium.cov <- cov(IR.ori.medium.closest.2.5kb)
IR.ori.medium.dens <- as.data.frame(table(IR.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(IR.ori.medium.cov)), IR.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, IR.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
IR.ori.high.closest.2.5kb <- IR.ori.closest.2.5kb %>% filter(class == "High")
IR.ori.high.cov <- cov(IR.ori.high.closest.2.5kb)
IR.ori.high.dens <- as.data.frame(table(IR.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(IR.ori.high.cov)), IR.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, IR.dens) %>% 
  mutate(class = "High")
# Bind all results
IR.ori.df <- rbind(IR.ori.low.dens, IR.ori.medium.dens, IR.ori.high.dens)
# Plot
IR.plot <- IR.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = IR.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("IR Percentage coverage") + ggtitle("Inverted repeats") + ylim(2.5,11) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# MR
MR.ori.closest <- read.csv("./Struct_analysis/MR.ori.closest.bed", sep = "\t", header = F)
colnames(MR.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "Spacer", "Permutations", "Subset", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select MR close to origins
MR.ori.closest.2.5kb <- MR.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
MR.ori.low.closest.2.5kb <- MR.ori.closest.2.5kb %>% filter(class == "Low")
MR.ori.low.cov <- cov(MR.ori.low.closest.2.5kb)
MR.ori.low.dens <- as.data.frame(table(MR.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(MR.ori.low.cov)), MR.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, MR.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins)
MR.ori.medium.closest.2.5kb <- MR.ori.closest.2.5kb %>% filter(class == "Medium")
MR.ori.medium.cov <- cov(MR.ori.medium.closest.2.5kb)
MR.ori.medium.dens <- as.data.frame(table(MR.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(MR.ori.medium.cov)), MR.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, MR.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
MR.ori.high.closest.2.5kb <- MR.ori.closest.2.5kb %>% filter(class == "High")
MR.ori.high.cov <- cov(MR.ori.high.closest.2.5kb)
MR.ori.high.dens <- as.data.frame(table(MR.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(MR.ori.high.cov)), MR.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, MR.dens) %>% 
  mutate(class = "High")
# Bind all results
MR.ori.df <- rbind(MR.ori.low.dens, MR.ori.medium.dens, MR.ori.high.dens)
# Plot
MR.plot <- MR.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = MR.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("MR Percentage coverage") + ggtitle("Mirror repeats") + ylim(0,25) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# STR
STR.ori.closest <- read.csv("./Struct_analysis/STR.ori.closest.bed", sep = "\t", header = F)
colnames(STR.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "Spacer", "Repeated", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select STR close to origins
STR.ori.closest.2.5kb <- STR.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
STR.ori.low.closest.2.5kb <- STR.ori.closest.2.5kb %>% filter(class == "Low")
STR.ori.low.cov <- cov(STR.ori.low.closest.2.5kb)
STR.ori.low.dens <- as.data.frame(table(STR.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(STR.ori.low.cov)), STR.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, STR.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
STR.ori.medium.closest.2.5kb <- STR.ori.closest.2.5kb %>% filter(class == "Medium")
STR.ori.medium.cov <- cov(STR.ori.medium.closest.2.5kb)
STR.ori.medium.dens <- as.data.frame(table(STR.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(STR.ori.medium.cov)), STR.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, STR.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
STR.ori.high.closest.2.5kb <- STR.ori.closest.2.5kb %>% filter(class == "High")
STR.ori.high.cov <- cov(STR.ori.high.closest.2.5kb)
STR.ori.high.dens <- as.data.frame(table(STR.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(STR.ori.high.cov)), STR.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, STR.dens) %>% 
  mutate(class = "High")
# Bind all results
STR.ori.df <- rbind(STR.ori.low.dens, STR.ori.medium.dens, STR.ori.high.dens)
# Plot
STR.plot <- STR.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = STR.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("STR Percentage coverage") + ggtitle("Short tandem Repeats") + ylim(0,5) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# Z
Z.ori.closest <- read.csv("./Struct_analysis/Z.ori.closest.bed", sep = "\t", header = F)
colnames(Z.ori.closest) <- c("chr", "start", "end", "Strand", "Length", "Repeat", "KVScore", "Subset", "Composition", "Sequence", "chr.1", "ori.pos", "ori.pos.1", "EFF", "ori.width", "ori.id", "class")
# Select Z close to origins
Z.ori.closest.2.5kb <- Z.ori.closest %>% mutate(start = start - ori.pos, end = end - ori.pos) %>% filter(start >= -2500 & start <= 2500)
# Compute coverage plot for
# Low efficient origins
Z.ori.low.closest.2.5kb <- Z.ori.closest.2.5kb %>% filter(class == "Low")
Z.ori.low.cov <- cov(Z.ori.low.closest.2.5kb)
Z.ori.low.dens <- as.data.frame(table(Z.ori.low.cov)) %>% mutate(dist = as.numeric(as.character(Z.ori.low.cov)), Z.dens = (Freq/7948)*100) %>% dplyr::select(dist, Freq, Z.dens) %>% 
  mutate(class = "Low")
# Medium efficient origins
Z.ori.medium.closest.2.5kb <- Z.ori.closest.2.5kb %>% filter(class == "Medium")
Z.ori.medium.cov <- cov(Z.ori.medium.closest.2.5kb)
Z.ori.medium.dens <- as.data.frame(table(Z.ori.medium.cov)) %>% mutate(dist = as.numeric(as.character(Z.ori.medium.cov)), Z.dens = (Freq/7972)*100) %>% dplyr::select(dist, Freq, Z.dens) %>% 
  mutate(class = "Medium")
# High efficient origins
Z.ori.high.closest.2.5kb <- Z.ori.closest.2.5kb %>% filter(class == "High")
Z.ori.high.cov <- cov(Z.ori.high.closest.2.5kb)
Z.ori.high.dens <- as.data.frame(table(Z.ori.high.cov)) %>% mutate(dist = as.numeric(as.character(Z.ori.high.cov)), Z.dens = (Freq/7985)*100) %>% dplyr::select(dist, Freq, Z.dens) %>% 
  mutate(class = "High")
# Bind all results
Z.ori.df <- rbind(Z.ori.low.dens, Z.ori.medium.dens, Z.ori.high.dens)
# Plot
Z.plot <- Z.ori.df %>% filter(class == "Low" | class == "Medium" | class == "High") %>% ggplot(aes(x = dist, y = Z.dens, color = class)) +
  geom_line(size = 0.5) +
  scale_color_manual(values=c("#F21A00", "#3B9AB2", "#EBCC2A")) +
  xlab("Distance from origin (bp)") + ylab("Z Percentage coverage") + ggtitle("Z-DNA") + ylim(0,1) + xlim(-1500,1500) +
  theme_bw() + theme(aspect.ratio=1)

# Save plot

pdf("./Rplot/Figure_S4B.pdf", width=20, height=8, useDingbats=FALSE)
ggarrange(APR.plot, DR.plot, GQ.plot, IR.plot, MR.plot, STR.plot, Z.plot, ncol = 4, nrow = 2)
dev.off()

#########################################################################################################
#########################################################################################################
# Compute feature coverage or distribution at origins
#########################################################################################################
#########################################################################################################

# Coverage/distribution were computed with deepTools computeMatrix function

# Here we report an example for computing the coverage associated with DNase-seq around origins

# using the dataset from ENCODE
# DNase-seq in H9 cells (G1 phase) - hg38
# ENCFF172SIB.bigWig

# The following commands are run from the terminal

computeMatrix scale-regions -S ./Dataset/DNase-seq/ENCFF172SIB.bigWig \
  -R ./Dataset/ori.high.bed ./Dataset/ori.medium.bed ./Dataset/ori.low.bed \
  --beforeRegionStartLength 10000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 10000 \
  --skipZeros -o ./deepTools/DNase-seq_H9.mat.gz
# The result can be visualised:
plotHeatmap -m ./deepTools/DNase-seq_H9.mat.gz \
  --colorMap YlOrBr \
  --heatmapHeight 10 \
  --startLabel "start" \
  --endLabel "end" \
  --zMin 0 \
  --zMax 1 \
  -out ./deepTools/DNase-seq_H9.pdf

# The results from plotHeatmap were used to build Figure 4C and S4C

# Similar analyses are run (using the same parameters) using the following datasets:

# ENCODE H3K4me1 in H9 cells - hg38 - ENCFF749WSP.bigWig
# ENCODE H3K79me2 in H9 cells - hg38 - ENCFF974LDP.bigWig
# ENCODE H3K9ac in H9 cells - hg38 - ENCFF936ETS.bigWig
# ENCODE H3K36me3 in H9 cells - hg38 - ENCFF739CKD.bigWig
# ENCODE H2AFZ in H9 cells - hg38 - ENCFF905SSM.bigWig
# ENCODE H3K27ac in H9 cells - hg38 - ENCFF429TIE.bigWig
# ENCODE H4K20me1 in H9 cells - hg38 - ENCFF548ENQ.bigWig
# ENCODE H3K27me3 in H9 cells - hg38 - ENCFF342DVJ.bigWig
# ENCODE H3K9me3 in H9 cells - hg38 - ENCFF456MYY.bigWig
# non B-DNA structure - hg38 - GQ.hg38.bw
# non B-DNA structure - hg38 - MR.hg38.bw
# non B-DNA structure - hg38 - IR.hg38.bw
# non B-DNA structure - hg38 - DR.hg38.bw
# non B-DNA structure - hg38 - STR.hg38.bw
# non B-DNA structure - hg38 - Z.hg38.bw

# The parameters for base composition and skews in nucleotides features were slightly different to consider shorter flanking sequences
# Example for GC content

computeMatrix scale-regions -S ./Dataset/BaseComp/GC.hg38.bw \
  -R ./Dataset/ori.high.bed ./Dataset/ori.medium.bed ./Dataset/ori.low.bed \
  --beforeRegionStartLength 2000 \
  --regionBodyLength 3000 \
  --afterRegionStartLength 2000 \
  --skipZeros -o ./deepTools/GC.mat.gz
 
# Similar analyses are run (using the same parameters) using the following datasets:

# Nucleotide composition - GC.skew.hg38.bw
# Nucleotide composition - AT.skew.hg38.bw
# Nucleotide composition - CpG.island.hg38.bw

#########################################################################################################
#########################################################################################################
# Compute correlation coefficients for all pairwise combinations of feature
#########################################################################################################
#########################################################################################################

# Open deepTools results and compute mean values
# The regions selected for computing mean values have been defined by visual inspection of the deepTools plotHeatmap results

# DNase
DNase <- read.table(gzfile("./deepTools/DNase-seq_H9.mat.gz"), skip = 1)   
DNase.val <- DNase %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
DNase.val.ori <- DNase.val[,1000:1500]
DNase.val.ori.2 <- rowMeans(DNase.val.ori)
# prepare df
DNase.df <- cbind.data.frame(ori.id = DNase$V4, DNase = DNase.val.ori.2)

# H2AZ
H2AZ <- read.table(gzfile("./deepTools/H2AZ_H9.mat.gz"), skip = 1)   
H2AZ.val <- H2AZ %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at origin boundaries (start +- 2000, end +- 2000)
H2AZ.val.ori <- H2AZ.val[,c(800:1200, 1300:1700)]
H2AZ.val.ori.2 <- rowMeans(H2AZ.val.ori)
# prepare df
H2AZ.df <- cbind.data.frame(ori.id = H2AZ$V4, H2AZ = H2AZ.val.ori.2)

# H3K9ac
H3K9ac <- read.table(gzfile("./deepTools/H3K9ac_H9.mat.gz"), skip = 1)   
H3K9ac.val <- H3K9ac %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
H3K9ac.val.ori <- H3K9ac.val[,1000:1500]
H3K9ac.val.ori.2 <- rowMeans(H3K9ac.val.ori)
# prepare df
H3K9ac.df <- cbind.data.frame(ori.id = H3K9ac$V4, H3K9ac = H3K9ac.val.ori.2)

# H3K27ac
H3K27ac <- read.table(gzfile("./deepTools/H3K27ac_H9.mat.gz"), skip = 1)   
H3K27ac.val <- H3K27ac %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at origin boundaries (start +- 2000, end +- 2000)
H3K27ac.val.ori <- H3K27ac.val[,c(800:1200, 1300:1700)]
H3K27ac.val.ori.2 <- rowMeans(H3K27ac.val.ori)
# prepare df
H3K27ac.df <- cbind.data.frame(ori.id = H3K27ac$V4, H3K27ac = H3K27ac.val.ori.2)

# H4K20me1
H4K20me1 <- read.table(gzfile("./deepTools/H4K20me1_H9.mat.gz"), skip = 1)   
H4K20me1.val <- H4K20me1 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
H4K20me1.val.ori <- H4K20me1.val[,1000:1500]
H4K20me1.val.ori.2 <- rowMeans(H4K20me1.val.ori)
# prepare df
H4K20me1.df <- cbind.data.frame(ori.id = H4K20me1$V4, H4K20me1 = H4K20me1.val.ori.2)

# H3K4me1
H3K4me1 <- read.table(gzfile("./deepTools/H3K4me1_H9.mat.gz"), skip = 1)   
H3K4me1.val <- H3K4me1 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at origin boundaries (start +- 2000, end +- 2000)
H3K4me1.val.ori <- H3K4me1.val[,c(800:1200, 1300:1700)]
H3K4me1.val.ori.2 <- rowMeans(H3K4me1.val.ori)
# prepare df
H3K4me1.df <- cbind.data.frame(ori.id = H3K4me1$V4, H3K4me1 = H3K4me1.val.ori.2)

# H3K27me3
H3K27me3 <- read.table(gzfile("./deepTools/H3K27me3_H9.mat.gz"), skip = 1)   
H3K27me3.val <- H3K27me3 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at origin boundaries (start +- 2000, end +- 2000)
H3K27me3.val.ori <- H3K27me3.val[,c(800:1200, 1300:1700)]
H3K27me3.val.ori.2 <- rowMeans(H3K27me3.val.ori)
# prepare df
H3K27me3.df <- cbind.data.frame(ori.id = H3K27me3$V4, H3K27me3 = H3K27me3.val.ori.2)

# H3K79me2
H3K79me2 <- read.table(gzfile("./deepTools/H3K79me2_H9.mat.gz"), skip = 1)   
H3K79me2.val <- H3K79me2 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at origin boundaries (start +- 2000, end +- 2000)
H3K79me2.val.ori <- H3K79me2.val[,c(800:1200, 1300:1700)]
H3K79me2.val.ori.2 <- rowMeans(H3K79me2.val.ori)
# prepare df
H3K79me2.df <- cbind.data.frame(ori.id = H3K79me2$V4, H3K79me2 = H3K79me2.val.ori.2)

# H3K36me3
H3K36me3 <- read.table(gzfile("./deepTools/H3K36me3_H9.mat.gz"), skip = 1)   
H3K36me3.val <- H3K36me3 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
H3K36me3.val.ori <- H3K36me3.val[,1000:1500]
H3K36me3.val.ori.2 <- rowMeans(H3K36me3.val.ori)
# prepare df
H3K36me3.df <- cbind.data.frame(ori.id = H3K36me3$V4, H3K36me3 = H3K36me3.val.ori.2)

# H3K9me3
H3K9me3 <- read.table(gzfile("./deepTools/H3K9me3_H9.mat.gz"), skip = 1)   
H3K9me3.val <- H3K9me3 %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
H3K9me3.val.ori <- H3K9me3.val[,1000:1500]
H3K9me3.val.ori.2 <- rowMeans(H3K9me3.val.ori)
# prepare df
H3K9me3.df <- cbind.data.frame(ori.id = H3K9me3$V4, H3K9me3 = H3K9me3.val.ori.2)

# GQ
GQ <- read.table(gzfile("./deepTools/GQ.mat.gz"), skip = 1)   
GQ.val <- GQ %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
GQ.val.ori <- GQ.val[,200:500]
GQ.val.ori.2 <- rowMeans(GQ.val.ori)
# prepare df
GQ.df <- cbind.data.frame(ori.id = GQ$V4, GQ = GQ.val.ori.2)

# MR
MR <- read.table(gzfile("./deepTools/MR.mat.gz"), skip = 1)   
MR.val <- MR %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
MR.val.ori <- MR.val[,200:500]
MR.val.ori.2 <- rowMeans(MR.val.ori)
# prepare df
MR.df <- cbind.data.frame(ori.id = MR$V4, MR = MR.val.ori.2)

# IR
IR <- read.table(gzfile("./deepTools/IR.mat.gz"), skip = 1)   
IR.val <- IR %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
IR.val.ori <- IR.val[,200:500]
IR.val.ori.2 <- rowMeans(IR.val.ori)
# prepare df
IR.df <- cbind.data.frame(ori.id = IR$V4, IR = IR.val.ori.2)

# DR
DR <- read.table(gzfile("./deepTools/DR.mat.gz"), skip = 1)   
DR.val <- DR %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
DR.val.ori <- DR.val[,200:500]
DR.val.ori.2 <- rowMeans(DR.val.ori)
# prepare df
DR.df <- cbind.data.frame(ori.id = DR$V4, DR = DR.val.ori.2)

# STR
STR <- read.table(gzfile("./deepTools/STR.mat.gz"), skip = 1)   
STR.val <- STR %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
STR.val.ori <- STR.val[,200:500]
STR.val.ori.2 <- rowMeans(STR.val.ori)
# prepare df
STR.df <- cbind.data.frame(ori.id = STR$V4, STR = STR.val.ori.2)

# Z
Z <- read.table(gzfile("./deepTools/Z.mat.gz"), skip = 1)   
Z.val <- Z %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values between start and end of origins
Z.val.ori <- Z.val[,200:500]
Z.val.ori.2 <- rowMeans(Z.val.ori)
# prepare df
Z.df <- cbind.data.frame(ori.id = Z$V4, Z = Z.val.ori.2)

# GC content
GC <- read.table(gzfile("./deepTools/GC.mat.gz"), skip = 1)   
GC.val <- GC %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average values at the center of origins
GC.val.ori <- GC.val[,300:400]
GC.val.ori.2 <- rowMeans(GC.val.ori)
# prepare df
GC.df <- cbind.data.frame(ori.id = GC$V4, GC = GC.val.ori.2)

# GC skew
GC.skew <- read.table(gzfile("./deepTools/GC.skew.mat.gz"), skip = 1)   
GC.skew.val <- GC.skew %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute the difference between average values computed at the beginning and end of origins
GC.skew.val.ori.1 <- GC.skew.val[-1,200:500]
# Fit line and extract slope
GC.skew.slope <- vector()
for (i in 1:nrow(GC.skew.val.ori.1)) {
  print(i/nrow(GC.skew.val.ori.1))
  df <- cbind.data.frame(dist = c(1:301), val = as.numeric(as.vector(GC.skew.val.ori.1[i,])))
  linearMod <- lm(dist ~ val, data=df)
  p <- linearMod$coefficients[[2]]
  GC.skew.slope <- append(GC.skew.slope, p)
}
# prepare df
GC.skew.df <- cbind.data.frame(ori.id = GC.skew[-1,]$V4, `GC skew` = GC.skew.slope)

# AT skew
AT.skew <- read.table(gzfile("./deepTools/AT.skew.mat.gz"), skip = 1)   
AT.skew.val <- AT.skew %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute the difference between average values computed at the beginning and end of origins
AT.skew.val.ori.1 <- AT.skew.val[-1,200:500]
# Fit line and extract slope
AT.skew.slope <- vector()
for (i in 1:nrow(AT.skew.val.ori.1)) {
  print(i/nrow(AT.skew.val.ori.1))
  df <- cbind.data.frame(dist = c(1:301), val = as.numeric(as.vector(AT.skew.val.ori.1[i,])))
  linearMod <- lm(dist ~ val, data=df)
  p <- linearMod$coefficients[[2]]
  AT.skew.slope <- append(AT.skew.slope, p)
}
# prepare df
AT.skew.df <- cbind.data.frame(ori.id = AT.skew[-1,]$V4, `AT skew` = -AT.skew.slope)

# CpG island
CpG <- read.table(gzfile("./deepTools/CpG.island.mat.gz"), skip = 1)   
CpG.val <- CpG %>% dplyr::select(-V1,-V2,-V3,-V4,-V5,-V6)
# compute average value between start and end of origins
CpG.val.ori.1 <- CpG.val[,200:500]
CpG.val.ori.2 <- rowMeans(CpG.val.ori.1)
# prepare df
CpG.df <- cbind.data.frame(ori.id = CpG$V4, `CpG islands` = CpG.val.ori.2)

# Combine all information

ori.chrom <- ori %>% left_join(DNase.df, by = "ori.id") %>% left_join(H2AZ.df, by = "ori.id") %>% left_join(H3K9ac.df, by = "ori.id") %>%
  left_join(H3K27ac.df, by = "ori.id") %>% left_join(H4K20me1.df, by = "ori.id") %>% left_join(H3K4me1.df, by = "ori.id") %>%
  left_join(H3K27me3.df, by = "ori.id") %>% left_join(H3K79me2.df, by = "ori.id") %>% left_join(H3K36me3.df, by = "ori.id") %>%
  left_join(H3K9me3.df, by = "ori.id") %>% left_join(GQ.df, by = "ori.id") %>% left_join(MR.df, by = "ori.id") %>%
  left_join(IR.df, by = "ori.id") %>% left_join(DR.df, by = "ori.id") %>% left_join(STR.df, by = "ori.id") %>% 
  left_join(Z.df, by = "ori.id") %>% left_join(GC.df, by = "ori.id") %>% left_join(GC.skew.df, by = "ori.id") %>%
  left_join(AT.skew.df, by = "ori.id") %>% left_join(CpG.df, by = "ori.id")

# Transform coverage values in log scale to optimise distributions toward gaussians

ori.chrom.log <- ori.chrom %>% mutate(DNase = log10(DNase), H2AZ = log10(H2AZ), H3K9ac = log10(H3K9ac), H3K27ac = log10(H3K27ac),
                                      H4K20me1 = log10(H4K20me1), H3K4me1 = log10(H3K4me1), H3K27me3 = log10(H3K27me3), H3K79me2 = log10(H3K79me2),
                                      H3K36me3 = log10(H3K36me3), H3K9me3 = log10(H3K9me3))

# Save table

write.csv(ori.chrom.log, "./Ori.coverage.log.csv", row.names = F)

# Load

ori.chrom <- read.csv("./Ori.coverage.log.csv")

# Replace infinite values by NA

ori.chrom[sapply(ori.chrom, is.infinite)] <- NA

# Test correlations

ori.chrom.cor <- ori.chrom %>% dplyr::select(-chr, -start, -end, -ori.id)
res <- cor(ori.chrom.cor, use = "complete.obs") # Pearson correlation

# Plot individual correlations

a <- ori.chrom %>% ggplot(aes(x = EFF, y = GC)) +
  geom_point(alpha = 0.05) +
  geom_smooth(colour = "#F21A00") + ylim(0.5,0.9) + xlab("Efficiency") + ylab ("%GC") + ggtitle("Rho = 0.677148569") +
  theme_bw() + theme(aspect.ratio=1)
b <- ori.chrom %>% ggplot(aes(x = EFF, y = GQ)) +
  geom_point(alpha = 0.05) +
  geom_smooth(colour = "#F21A00") + ylim(0,0.3) + xlab("Efficiency") + ylab ("GQ density") + ggtitle("Rho = 0.2000336691") +
  theme_bw() + theme(aspect.ratio=1)
c <- ori.chrom %>% ggplot(aes(x = EFF, y = H3K36me3)) +
  geom_point(alpha = 0.05) +
  geom_smooth(colour = "#3B9AB2") + ylim(-1.5,1) + xlab("Efficiency") + ylab ("H3K36me3 coverage (log10)") + ggtitle("Rho = -0.375923598") +
  theme_bw() + theme(aspect.ratio=1)
# save plots
pdf("./Rplot/Figure_4D.pdf", width=10, height=4, useDingbats=FALSE)
ggarrange(a, b, c, nrow = 1, ncol = 3)
dev.off()

# Plot heatmap reporting pairwise correlations

# Define color palette
col <- colorRampPalette(c("#3B9AB2", "white", "#F21A00"))(n = 299)
# Define the color breaks manually
col_breaks = c(seq(-0.5,-0.1,length=100),  # for blue
               seq(-0.099,0.099,length=100), # for white
               seq(0.1,1,length=100)) # for red
# Define distance and clustering methods
distance <- dist(res, method = "manhattan")
cluster <- hclust(distance, method = "single")
# save heatmap
pdf("./Rplot/Figure_4E.pdf", width=10, height=10, useDingbats=FALSE)
heatmap.2(x = res, col = col, breaks= col_breaks, symm = TRUE, trace="none",
          Rowv = as.dendrogram(cluster), # apply default clustering method
          Colv = as.dendrogram(cluster)) # apply default clustering method
dev.off()

#########################################################################################################
#########################################################################################################
# Principal component analysis
#########################################################################################################
#########################################################################################################

set.seed(825)

# Centering, scaling and splitting the data 

ori.chrom.2 <- ori.chrom %>% dplyr::select(-chr, -start, -end, -ori.id) %>% na.omit()
preProcValues <- preProcess(ori.chrom.2, method = c("center", "scale"))
mydata <- predict(preProcValues, ori.chrom.2)

# Principal component analysis

mydata.2 <- mydata %>% dplyr::select(-EFF)
data.pca <- prcomp(mydata.2, scale = TRUE)

# Contributions of the variables to the principal components

e <- fviz_pca_var(data.pca, col.var="contrib") +
  scale_color_gradient2(low="white", mid="blue", 
                        high="red", midpoint=5) + theme_minimal()

# Graph of individuals

pca.individual.df <- cbind.data.frame(PC1 = data.pca$x[,1], PC2 = data.pca$x[,2], EFF = mydata$EFF)
f <- pca.individual.df %>% 
  ggplot(aes(x=PC1, y=PC2, color = EFF) ) +
  geom_point(alpha = 0.8, size = 1.5) + xlim (-6,5) + ylim(-8.5,7.5) + xlab("PC1 - 19.1%") + ylab("PC2 - 11.8%") +
  scale_colour_gradient(low = "#3B9AB2", high = "#F21A00", na.value = NA) +
  theme_bw() + theme(aspect.ratio=1)

# Save plots

pdf("./Rplot/Figure_4F.pdf", width=10, height=5, useDingbats=FALSE)
ggarrange(e, f, ncol = 2, nrow = 1)
dev.off()

#########################################################################################################
#########################################################################################################
# Statistical Model
#########################################################################################################
#########################################################################################################

# Identify and remove higlhy correlated predictors (cutoff = 0.85)

descrCor <- cor(mydata)
highlyCorDescr <- findCorrelation(descrCor, cutoff = .85)
#  --> Nothing filtered

# Identify Linear Dependencies

comboInfo <- findLinearCombos(mydata)  ## --> NULL

# Centering, Scaling and Splitting the data 

split <- createDataPartition(y=mydata$EFF, p = 0.7, list=FALSE)
train <- mydata[split,]
rest <- mydata[-split,]
split.2 <- createDataPartition(y=rest$EFF, p = 0.333, list=FALSE)
test.1 <- rest[split.2,]
rest.1 <- rest[-split.2,]
split.3 <- createDataPartition(y=rest.1$EFF, p = 0.5, list=FALSE)
test.2 <- rest.1[split.3,]
test.3 <- rest.1[-split.3,]
all.test <- rbind(test.1, test.2, test.3)
# Different testing sets were generated for model selection
# Selected a SVM algorithm based on performance and accuracy

# Optimise a Support Vector Machine with Radial Basis Function Kernel algorithm on all 20 predictors

# Found optimal sigma (Kernel function) and C (penalty parameter) values

svmGrid <-  expand.grid(sigma = 0.05*c(0.1, 0.2, 0.3, 0.4, 0.5), 
                          C = c(1.30, 1.35, 1.4, 1.45, 1.5))
svmFit <- train(EFF ~., data = train,
                  method = "svmRadial",
                  trControl = trainControl(method="cv", number=10),
                  tuneGrid = svmGrid,
                  verbose = FALSE)
# The final values used for the model were sigma = 0.015 and C = 1.3

# Save model

saveRDS(svmFit, file = "./svmFit.rds") # This is the final model

# Load model

svmFit <- readRDS(file = "./svmFit.rds")

# Plot model result

train.model <- cbind.data.frame(Experimental = train$EFF, Predicted = predict(svmFit, newdata = train), class = "train")
test.model <- cbind.data.frame(Experimental = all.test$EFF, Predicted = predict(svmFit, newdata = all.test), class = "test")
model.df <- rbind(train.model, test.model)

cor.test(train.model$Experimental, train.model$Predicted) # Rho = 0.8507047
cor.test(test.model$Experimental, test.model$Predicted) # Rho = 0.8418099
r2(train.model$Experimental, train.model$Predicted) # r2 = 0.7234197
r2(test.model$Experimental, test.model$Predicted) # r2 = 0.7072211

model.plot <- model.df %>%
  ggplot(aes(x=Experimental, y=Predicted, color = class) ) +
  scale_color_manual(values=c("#F21A00", "#999999")) +
  geom_point(alpha = 0.4, size = 1.5) + xlim (-2.05,2.05) + ylim(-2.05,2.05) + xlab("Experimental efficiency") + ylab("Predicted efficiency") +
  ggtitle("Train r2 = 0.7234197 | Rho = 0.8507047\nTest r2 = 0.7072211 | Rho = 0.8418099") +
  theme_bw() + theme(aspect.ratio=1) + geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# Plot contribution

svmImp <- varImp(svmFit, scale = TRUE)
cont <- as.data.frame(svmImp$importance)
cont$variable <- row.names(cont)
cont$class <- c("accessibility", "accessibility", "active", "active", "active", "active", "inactive", "inactive", "active", "inactive", "structure", "structure", "structure", "structure", "structure", "structure", "base composition", "base composition", "base composition", "base composition")

contribution.plot <- cont %>%
  ggplot(aes(x=reorder(variable, Overall), y=Overall, fill = class) ) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values=c("#3B9AB2", "#78B7C5", "#EBCC2A", "#F21A00", "#E1AF00")) + ylab("Contribution") +
  theme_bw() + theme(aspect.ratio=1) +
  coord_flip() +
  theme(axis.title.y = element_blank())

# Save Figures

pdf("./Rplot/Figure_4G_4H.pdf", width=12, height=5, useDingbats=FALSE)
ggarrange(model.plot, contribution.plot, nrow = 1, ncol = 2)
dev.off()



























