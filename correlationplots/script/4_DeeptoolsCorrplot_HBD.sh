################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Reproducibility assessment                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create deeptools heatmap for reproducibility assessment for         #
# Celeste"s & Serkan's data                                                    #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
namesummary=brigther1_2
nameplot=brighter1_2
outPath=/project/ChromGroup_Seq_data/Celeste/analysis
#-------------------------Script-----------------------------------------------#

#deeptools correlation map

bamfiles=()
samples=()
for SampleName in HBD.igg_1 HBD.igg_2 HBD.his_1 HBD.his_2 HBD.rpa70_2 HBD.h3k4me3_1 \
HBD.h3k4me3_2
do
   bamfilefiltered=$projPathC/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bamfilefiltered"
done 


for SampleName in Flag.mo_1 Flag.mo_2 
do
   bamfilefiltered=$projPathS/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bamfilefiltered"
done 


multiBamSummary bins --bamfiles $bamfiles -l IgG_1 IgG_2 BRIGHTER1_His_1 BRIGHTER1_His_2 RPA70_2 \
H3K4me3_1 H3K4me3_2 BRIGHTER2_Flag_1 BRIGHTER2_Flag_2 \
 --binSize 2000 \
  --blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
  -o $outPath/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/${namesummary}_2000bp.npz -c pearson  \
--colorMap bwr --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/${nameplot}_2000bp.cormatrix.txt \
