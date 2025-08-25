################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Reproducibility assessment                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create deeptools heatmap to show differences between Flag, GFP and  #
# His                                                                          #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
experimentSummaryC=$projPathC/experiment_summary_align_formatted.csv
namesummary=brighter_flag_his_gfp
nameplot=brighter_flag_his_gfp
outPath=/project/ChromGroup_Seq_data/Rloop_paper/correlationplots
#-------------------------Script-----------------------------------------------#

#populate bam files array
bamfiles=()

for SampleName in HBD.flag_1 HBD.flag_2 HBD.gfp_1 HBD.gfp_2 HBD.his_1 HBD.his_2 \
HBD.igg_1 HBD.igg_2; do
bam=$projPathC/alignment/bam/$SampleName.mapped.sorted.filtered.bam
bamfiles+=" $bam"
done

samples=("BRIGHTER1_Flag_1" "BRIGHTER1_Flag_2" "BRIGHTER1_GFP_1" "BRIGHTER1_GFP_2" \
"BRIGHTER1_His_1" "BRIGHTER1_His_2" "BRIGHTER1_IgG_1" "BRIGHTER1_IgG_2")

multiBamSummary bins --bamfiles $bamfiles -l ${samples[@]} --binSize 2000 \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-o $outPath/data/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/data/${namesummary}_2000bp.npz -c pearson  \
--colorMap bwr --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/plot/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/data/${nameplot}_2000bp.cormatrix.txt\