################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Analysis                                                                     #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: comparison of the different Flag experimen  plots                   #
################################################################################


library(ggplot2)
library(viridis)
##Venn Flag mo v Flag rab
overlaps_mo_rab.path <- "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/overlaps_Flagmo_Flagrab.txt"

overlaps_mo_rab <- read.table(overlaps_mo_rab.path, header=T)

countvec <- overlaps_mo_rab$OverlapCounts
names(countvec) <- overlaps_mo_rab$OverlapName
venn <- euler(countvec)
plot(venn, quantities = TRUE)

## Venn Flag Serkan v Flag Celeste
overlaps_S_C.path <- "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/overlaps_FlagS_FlagC.txt"

overlaps_S_C<- read.table(overlaps_S_C.path, header=T)

countvec <- overlaps_S_C$OverlapCounts
names(countvec) <- overlaps_S_C$OverlapName
venn <- euler(countvec)
plot(venn, quantities = TRUE)


## Stacked barplot genomic locations

locationpeaks.path <-  "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/genomiclocations_FlagS_FlagC.txt"

locationpeaks <- read.table(locationpeaks.path, header=T)
locationpeaks
ggplot(locationpeaks, aes(fill=Feature, y=Frequency, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  theme_light() +
  theme(axis.title.x=element_blank())+
  ylab("% peaks overlapping regions")+
  scale_fill_viridis(discrete=T) 


locationpeakssi.path <- "/project/ChromGroup_Seq_data/Rloop_paper/Flag_comparisons/data/silencing_peak_annot.txt"
locationpeakssi <- read.table(locationpeakssi.path, header = T)
locationpeakssi
ggplot(locationpeakssi, aes(fill=Feature, y=Frequency, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  theme_light() +
  theme(axis.title.x=element_blank())+
  ylab("Percentatge (%)")+
  scale_fill_viridis(discrete=T) 
