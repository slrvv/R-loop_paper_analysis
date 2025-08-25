################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Reproducibility assessment                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create deeptools heatmap to show differences between Flag, GFP and  #
# His between WKK, WT, HBD                                                     #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09

#-------------------------Script-----------------------------------------------#

namesummary=hbdhis_hm_experiments_bam
nameplot=hbdhis_hm_experiments_corr_heatmap
outPath=/project/ChromGroup_Seq_data/Rloop_paper/correlationplots

#populate bam files array
bamfiles=()
samples=()
for SampleName in HBD.his_1 HBD.his_2 HBD.igg_1 
do
   bam=$projPathC/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bam"
   samples+=" ${SampleName}"
done

for SampleName in H3K4me1_1 H3K4me1_2 H3K4me3_1 H3K4me3_2 H3K27ac_1 H3K27ac_2
do
   bam=$projPathS/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bam"
   samples+=" ${SampleName}"
done


echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
  -o $outPath/data/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/data/${namesummary}_2000bp.npz -c pearson  \
--colorMap bwr --plotNumbers --removeOutliers --plotHeight 18 --plotWidth 20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/plot/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/data/${nameplot}_2000bp.cormatrix.txt 


