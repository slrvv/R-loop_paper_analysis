##plotting the % of rDNA in fastq file

library(readr)
library(ggplot2)
##1. gather the % of rDNA in the  fastq file from stats files

read_stats <- function(sample){
  filename <- paste0("/project/ChromGroup_Seq_data/Rloop_paper/rDNA_quantification/data/",
                     sample,
                     "_matched_stats.txt")
  line <- read_lines(filename, 
                     skip=2, n_max=1)
  percentatge <- strsplit(line, "\t")[[1]][3]
  return(percentatge)
}

sampleL <- c("Flag.mo", "Cntsi.Flag.m", "Topo1si.fl", "RNAseH1si.fl")
repL <- c("_1", "_2")
percentL <- c()
samplerepL <- c()
for (exp in sampleL){
  for (rep in repL){
    sample <- paste0(exp, rep)
    percentatge <- read_stats(sample)
    percentL <- c(percentL, percentatge)
    samplerepL <- c(samplerepL, sample)
  } 
}

percentDat <- data.frame(Sample=samplerepL, 
                         Percentatge_rDNA=percentL)

percent_cont <- unlist(strsplit(percentDat$Percentatge_rDNA, "%"))

percentDat$Percentatge_rDNA <- as.numeric(percent_cont)
percentDat$SamplenoRep <- sapply(strsplit(percentDat$Sample, "_"), "[[", 1)
ggplot(percentDat, aes(x=Sample, y=Percentatge_rDNA, fill= SamplenoRep)) + 
  geom_bar(stat = "identity")+
  theme_light() +
  scale_fill_brewer(palette = "Blues")
  theme(axis.title.x=element_blank(),
        axis.text=element_text(size=14),
        axis.title.y = element_text(size=16))+
  guides(fill="none") +
  ylab("% of reads matching rDNA")
