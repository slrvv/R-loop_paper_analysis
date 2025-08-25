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
experimentSummaryC=$projPathC/experiment_summary_align_formatted.csv
namesummary=flag_his_gfp_experiments_bam_Celeste
nameplot=flag_his_gfp_experiments_corr_heatmap_Celeste
outPath=/project/ChromGroup_Seq_data/Rloop_paper/correlationplots
#-------------------------Script-----------------------------------------------#

#populate bam files array
bamfiles=()
samples=()
for SampleName in HBD.flag_1 HBD.flag_2 HBD.gfp_1 HBD.gfp_2 HBD.his_1 HBD.his_2 \
WKK.flag_1 WKK.flag_2 WKK.gfp_1 WKK.gfp_2 WKK.his_1 WKK.his_2 \
WT.flag_1 WT.flag_2 WT.gfp_1 WT.gfp_2 WT.his_1 WT.his_2 \
HBD.igg_1 WKK.igg_1 WT.igg_1
do
   bam=$projPathC/alignment/bam/$SampleName.mapped.sorted.filtered.bam
   bamfiles+=" $bam"
   samples+=" ${SampleName}"
done


echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
 --blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
  -o $outPath/data/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/data/${namesummary}_2000bp.npz -c pearson  \
--colorMap bwr --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/plot/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/data/${nameplot}_2000bp.cormatrix.txt\


