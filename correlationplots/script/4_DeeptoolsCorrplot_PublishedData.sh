################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Reproducibility assessment                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Correlation plot with reprocessed S9.6 data                         #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
projPathPl=/project/ChromGroup_Seq_data/External_data/S9.6_cutntag/Reprocessing
namesummary=published_comp_bam
nameplot=published_comp_bam
outPath=/project/ChromGroup_Seq_data/Rloop_paper/correlationplots
#-------------------------Script-----------------------------------------------#
bamsuf="mapped.sorted.filtered.bam"

#populate bam files array
bamfiles=("$projPathC/alignment/bam/HBD.his_1.$bamsuf" \
"$projPathC/alignment/bam/HBD.his_2.$bamsuf" \
"$projPathC/alignment/bam/HBD.igg_1.$bamsuf" \
"$projPathC/alignment/bam/HBD.igg_2.$bamsuf" \
"$projPathS/alignment/bam/Flag.mo_1.$bamsuf" \
"$projPathS/alignment/bam/Flag.mo_2.$bamsuf" \
"$projPathPl/alignment/bam/GSM7009544.Rloop.CTR_1.mapped.sorted.bam" \
"$projPathPl/alignment/bam/GSM7009545.Rloop.CTR_2.mapped.sorted.bam")

samples=("BRIGHTER1.His_1" "BRIGHTER1.His_2" "BRIGHTER1.IgG_1" "BRIGHTER1.IgG_2" \
"BRIGHTER2.Flag_1" "BRIGHTER2.Flag_2" " GSM7009544.S9.6_1" "GSM7009545.S9.6_2")

echo ${bamfiles[@]}

multiBamSummary bins --bamfiles ${bamfiles[@]} -l ${samples[@]} --binSize 2000 \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
  -o $outPath/data/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/data/${namesummary}_2000bp.npz -c pearson  \
--colorMap bwr --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/plot/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/data/${nameplot}_2000bp.cormatrix.txt\


