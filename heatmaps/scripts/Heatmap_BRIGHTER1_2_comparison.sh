################################################################################
#                                                                              #
# cut & tag analysis R-loop project                                            #
#                                                                              #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Brighter 1 and 2 comparisons                                        #
################################################################################


#----------------------------Paths---------------------------------------------#
projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathCbw=$projPathC/alignment/bigwig
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
projPathSbw=$projPathS/alignment/bigwig
namesummary=brighter1_2_comparison
nameplot=brighter1_2_comparison
outpath=/project/ChromGroup_Seq_data/Rloop_paper/heatmaps
gencode=/project/ChromGroup_Seq_data/Celeste/analysis/gencode.v44.annotation.gtf

#-------------------------Script-----------------------------------------------#
bwsuf=bowtie2.fragments.normalized.filtered.bw
peaksuf=merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed
balcklist=$projPathC/GRCh38_unified_blacklist.bed

samples=("$projPathCbw/HBD.his_$bwsuf" "$projPathSbw/Flag.mo_$bwsuf" "$projPathCbw/HBD.h3k4me3_$bwsuf" "$projPathSbw/H3K4me3_$bwsuf" "$projPathCbw/HBD.rpa70_2_$bwsuf" "$projPathCbw/HBD.igg_$bwsuf")

echo ${samples[@]}
echo ${labels}

computeMatrix reference-point --referencePoint TSS -S ${samples[@]} \
-R $gencode \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 5000 \
-b 5000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_genes.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_genes.mat.gz \
--colorList 'white,#006666' 'white,#0000cc' 'white,#fcba03' 'white,#b266ff' 'white,#960096' 'white,#000000' \
--samplesLabel BRIGHTER1.His BRIGHTER2.Flag BRIGHTER1.H3K4me3 BRIGHTER2.H3K4me3 BRIGHTER1.RPA70_2 BRIGHTER1.IgG \
--boxAroundHeatmaps no \
--yMin 0.0 0.0 0.0 0.0 0.0 0.0 \
--yMax 0.4 2.5 1.5 200 0.4 0.4 \
--zMin 0.0 0.0 0.0 0.0 0.0 0.0 \
--zMax 0.2 2 1.5 200 0.2 0.2 \
--heatmapWidth 5 \
--heatmapHeight 20 \
-out $outpath/plot/${nameplot}_tss_heatmap_genes.png \



##Matrix over genes
computeMatrix scale-regions -S ${samples[@]} \
-R $gencode \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros -o $outpath/data/${namesummary}_genes.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_genes.mat.gz \
--samplesLabel HBD.BRIGHTER1 HBD.BRIGHTER2 HBD.BRIGHTER1.H3K4me3 HBD.BRIGHTER1.RPA70_2 HBD.BRIGHTER1.IgG \
--heatmapWidth 5 \
--heatmapHeight 20 \
-out $outpath/plot/${nameplot}_heatmap_genes.png \



##Matrix over genes
computeMatrix scale-regions -S ${samples[@]} \
-R $gencode \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros -o $outpath/data/${namesummary}_genes.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_genes.mat.gz \
--samplesLabel HBD.BRIGHTER1 HBD.BRIGHTER2 HBD.BRIGHTER1.H3K4me3 HBD.BRIGHTER1.RPA70_2 HBD.BRIGHTER1.IgG \
--heatmapWidth 5 \
--heatmapHeight 20 \
-out $outpath/plot/${nameplot}_heatmap_genes.png \

##Matrix over peaks of one and the other
computeMatrix reference-point --referencePoint center -S ${samples[@]} \
-R $projPathC/peakCalling/SEACR/HBD.his_$peaksuf $projPathS/peakCalling/SEACR/Flag.mo_$peaksuf \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 3000 \
-b 3000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_peaks.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_peaks.mat.gz \
--colorList 'white,#006666' 'white,#0000cc' 'white,#fcba03' 'white,#b266ff' 'white,#960096' 'white,#000000' \
--samplesLabel BRIGHTER1.His BRIGHTER2.Flag BRIGHTER1.H3K4me3 BRIGHTER2.H3K4me3 BRIGHTER1.RPA70_2 BRIGHTER1.IgG \
--boxAroundHeatmaps no \
--regionsLabel "BRIGHTER1.his peaks" "BRIGHTER2.flag peaks" \
--heatmapWidth 5 \
--heatmapHeight 20 \
--yMin 0.0 0.0 0.0 0.0 0.0 0.0 \
--yMax 0.5 10 4 500 0.5 0.4 \
--zMin 0.0 0.0 0.0 0.0 0.0 0.0 \
--zMax 0.5 8 4 500 0.3 0.2 \
-out $outpath/plot/${nameplot}_heatmap_refpoint_peaks.png \

