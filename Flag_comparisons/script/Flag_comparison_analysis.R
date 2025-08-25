################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Analysis                                                                     #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: comparison of the different Flag experimen                          #
################################################################################


suppressMessages(library("GenomicRanges"))
library(ChIPseeker)
library("ChIPpeakAnno")
library("eulerr")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

c.flag.path <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_complete/peakCalling/SEACR/HBD.flag_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"

c.flag <- read.table(c.flag.path)

s.flag.mo.path <- "/project/ChromGroup_Seq_data/Serkan/2024_09/peakCalling/SEACR/Flag.mo_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"
s.flag.rab.path <- "/project/ChromGroup_Seq_data/Serkan/2024_09/peakCalling/SEACR/Flag.rab_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"
s.flag.mo <- read.table(s.flag.mo.path)
s.flag.rab <- read.table(s.flag.rab.path)

c.flag.peaks <- nrow(c.flag)
s.flag.mo.peaks <- nrow(s.flag.mo)
s.flag.rab.peaks <- nrow(s.flag.rab)

print(paste0("Flag celeste has ", c.flag.peaks, " peaks"))
print(paste0("Flag Serkan Mouse has ", s.flag.mo.peaks, " peaks"))
print(paste0("Flag Serkan Rabbit has ", s.flag.rab.peaks, " peaks"))

##overlaps of Serkan's Flag mouse and Flag rabbit

s.flag.mo.gr = GRanges(s.flag.mo$V1, IRanges(start = s.flag.mo$V2, end = s.flag.mo$V3), strand = "*")
s.flag.rab.gr = GRanges(s.flag.rab$V1, IRanges(start = s.flag.rab$V2, end = s.flag.rab$V3), strand = "*")

overlaps <- findOverlapsOfPeaks(GRangesList(FLAG.MOUSE= s.flag.mo.gr , FLAG.RABBIT = s.flag.rab.gr),
                                minoverlap = 0.25,
                                connectedPeaks = "keepAll")

overlaps$venn_cnt

overlapsdat <- data.frame(OverlapName = c("FLAG.MOUSE", 
                                        "FLAG.RABBIT", 
                                        "FLAG.MOUSE&FLAG.RABBIT"),
                        OverlapCounts = c(4560, 19, 3309))

write.table(overlapsdat,
            "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/overlaps_Flagmo_Flagrab.txt",
            row.names = F)

## Overlap of Celeste's Flag and Serkan's Flag mouse
c.flag.gr = GRanges(c.flag$V1, IRanges(start = c.flag$V2, end = c.flag$V3), strand = "*")

overlaps <- findOverlapsOfPeaks(GRangesList(FLAG.S= s.flag.mo.gr , FLAG.C = c.flag.gr),
                                minoverlap = 0.25, connectedPeaks = "keepAll")

overlaps$venn_cnt

overlapsdat <- data.frame(OverlapName = c("FLAG.S", 
                                          "FLAG.C",
                                          "FLAG.S&FLAG.C"),
                          OverlapCounts = c(3711, 553, 4158))
write.table(overlapsdat,
            "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/overlaps_FlagS_FlagC.txt",
            row.names = F)


## Genomic Locations of both flag S and flag C
peakAnnoSerkan <- annotatePeak(s.flag.mo.gr, TxDb=txdb, annoDb="org.Hs.eg.db",
                            level = "transcript", 
                            genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic") )

peakAnnoCeleste <- annotatePeak(c.flag.gr, TxDb=txdb, annoDb="org.Hs.eg.db",
                               level = "transcript", 
                               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic") )
peakAnnoCeleste <- peakAnnoCeleste@annoStat
peakAnnoCeleste$Sample <- rep("Flag.C", nrow(peakAnnoCeleste))

peakAnnoAll <- rbind(peakAnnoSerkan, peakAnnoCeleste)

write.table(peakAnnoAll,
            "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/genomiclocations_FlagS_FlagC.txt",
            row.names = F)

plotAnnoBar(list(BRIGHTER1.Flag = peakAnnoCeleste,
                 BRIGHTER2.Flag = peakAnnoSerkan), 
            title = "Genomic location distribution")



## Annotation of peaks Topo1si & Flag 

cntsi.flag.path <- "/project/ChromGroup_Seq_data/Serkan/2024_09/peakCalling/SEACR/Cntsi.Flag.m_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"
topo1si.path <- "/project/ChromGroup_Seq_data/Serkan/2024_09/peakCalling/SEACR/Topo1si.fl_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"
rnaseh1si.path <- "/project/ChromGroup_Seq_data/Serkan/2024_09/peakCalling/SEACR/RNAseH1si.fl_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed"

cntsi.flag <- read.table(cntsi.flag.path)
topo1si <- read.table(topo1si.path)
rnaseh1si <- read.table(rnaseh1si.path)

cntsi.flag.gr = GRanges(cntsi.flag$V1, IRanges(start = cntsi.flag$V2, end = cntsi.flag$V3), strand = "*")
topo1si.gr = GRanges(topo1si$V1, IRanges(start = topo1si$V2, end = topo1si$V3), strand = "*")
rnaseh1si.gr = GRanges(rnaseh1si$V1, IRanges(start = rnaseh1si$V2, end = rnaseh1si$V3), strand = "*")

cntsi.peakanno <- annotatePeak(cntsi.flag.gr, TxDb=txdb, annoDb="org.Hs.eg.db",
                               level = "transcript", 
                               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic") )

topo1si.peakanno <- annotatePeak(topo1si.gr, TxDb=txdb, annoDb="org.Hs.eg.db",
                                 level = "transcript", 
                                 genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic") )

rnaseh1si.peakanno <- annotatePeak(rnaseh1si.gr, TxDb=txdb, annoDb="org.Hs.eg.db",
                                 level = "transcript", 
                                 genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic") )

plotAnnoBar(list(BRIGHTER2.CNTsi = cntsi.peakanno,
     BRIGHTER2.Topo1si = topo1si.peakanno,
     BRIGHTER2.RNASseH1si = rnaseh1si.peakanno), 
     title = "Genomic location distribution")
