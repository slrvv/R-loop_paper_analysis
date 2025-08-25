################################################################################
#                                                                              #
# cut & tag analysis R-loop project                                            #
#                                                                              #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Comparison to the reprocessed s9.6 data                             #
################################################################################


#----------------------------Paths---------------------------------------------#
projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathCbw=$projPathC/alignment/bigwig
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
projPathSbw=$projPathS/alignment/bigwig
projPathPlbw=/project/ChromGroup_Seq_data/External_data/S9.6_cutntag/Reprocessing/alignment/bigwig
projPathPl=/project/ChromGroup_Seq_data/External_data/S9.6_cutntag/Reprocessing
namesummary=published_comp
nameplot=published_com
outpath=/project/ChromGroup_Seq_data/Rloop_paper/heatmaps
gencode=/project/ChromGroup_Seq_data/Celeste/analysis/gencode.v44.annotation.gtf

#-------------------------Script-----------------------------------------------#
bwsuf=bowtie2.fragments.normalized.filtered.bw
peaksuf=merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed
balcklist=$projPathC/GRCh38_unified_blacklist.bed

samples=("$projPathSbw/Flag.mo_$bwsuf" "$projPathCbw/HBD.his_$bwsuf"  \
"$projPathPlbw/GSM7009544.Rloop.CTR_bowtie2.fragments.normalized.bw" "$projPathCbw/HBD.igg_$bwsuf")

echo ${samples[@]}

bedops -u $projPathS/peakCalling/SEACR/Flag.mo_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed \
$projPathPl/peakCalling/SEACR/GSM7009544.Rloop.CTR_merged_seacr_norm_control.peaks.relaxed.bed > $outpath/data/Flag_S9.6.bed

computeMatrix reference-point --referencePoint center -S ${samples[@]} \
-R $projPathS/peakCalling/SEACR/Flag.mo_merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed \
$projPathPl/peakCalling/SEACR/GSM7009544.Rloop.CTR_merged_seacr_norm_control.peaks.relaxed.bed \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 3000 \
-b 3000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_peaks.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_peaks.mat.gz \
--colorList 'white,#0000cc' 'white,#006666' 'white,#4c0099' 'white,#000000' \
--samplesLabel BRIGHTER2.Flag BRIGHTER1.His GSM7009545_S9.6 BRIGHTER1.IgG \
--regionsLabel "BRIGHTER2.Flag peaks" "GSM7009545_S9.6 peaks" \
--yMin 0.0 0.0 0.0 0.0 \
--yMax 8 0.4 2100 0.4 \
--zMin 0.0 0.0 0.0 0.0 \
--zMax 7 0.3 auto 0.4 \
--boxAroundHeatmaps no \
--heatmapWidth 5 \
--heatmapHeight 20 \
-out $outpath/plot/${nameplot}_center_heatmap_peaks.png \