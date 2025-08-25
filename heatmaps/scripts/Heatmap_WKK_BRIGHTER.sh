################################################################################
#                                                                              #
# cut & tag analysis R-loop project                                            #
#                                                                              #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Heatmap HBD.his and WKK.his comparison                              #
################################################################################


#----------------------------Paths---------------------------------------------#
projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathCbw=$projPathC/alignment/bigwig
namesummary=brighter1_wkk_comparison
nameplot=brighter1_wkk_comparison
outpath=/project/ChromGroup_Seq_data/Rloop_paper/heatmaps
gencode=/project/ChromGroup_Seq_data/Celeste/analysis/gencode.v44.annotation.gtf

#-------------------------Script-----------------------------------------------#
bwsuf=bowtie2.fragments.normalized.filtered.bw
peaksuf=merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed
balcklist=$projPathC/GRCh38_unified_blacklist.bed

samples=("$projPathCbw/HBD.his_$bwsuf"  "$projPathCbw/WKK.his_$bwsuf" \
"$projPathCbw/WT.his_$bwsuf" "$projPathCbw/HBD.igg_$bwsuf")

echo ${samples[@]}

bedtools intersect -a $projPathC/peakCalling/SEACR/HBD.his_$peaksuf \
-b $projPathC/peakCalling/SEACR/WKK.his_$peaksuf > $outpath/data/HBD_and_WKK_peaks.bed


bedtools intersect -v -a $projPathC/peakCalling/SEACR/HBD.his_$peaksuf \
-b $projPathC/peakCalling/SEACR/WKK.his_$peaksuf > $outpath/data/HBD_only_peaks.bed

bedtools intersect -v -a $projPathC/peakCalling/SEACR/WKK.his_$peaksuf   \
-b $projPathC/peakCalling/SEACR/HBD.his_$peaksuf> $outpath/data/WKK_only_peaks.bed

##Matrix over peaks of one and the other
computeMatrix reference-point --referencePoint center -S ${samples[@]} \
-R $outpath/data/HBD_only_peaks.bed $outpath/data/HBD_and_WKK_peaks.bed $outpath/data/WKK_only_peaks.bed \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 3000 \
-b 3000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_peaks.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_peaks.mat.gz \
--colorList 'white,#006666' 'white,#00cccc' 'white,#009999' 'white,#000000' \
--samplesLabel BRIGHTER1_His WKK_His WT_His BRIGHTER1_IgG \
--boxAroundHeatmaps no \
--regionsLabel "BRIGHTER1_His peaks" "BRIGHTER1_His + WKK_His peaks" "WKK_His peaks" \
-out $outpath/plot/${nameplot}_heatmap_refpoint_peaks.png \

