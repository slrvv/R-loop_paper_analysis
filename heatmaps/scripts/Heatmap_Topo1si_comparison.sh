################################################################################
#                                                                              #
# cut & tag analysis R-loop project                                            #
#                                                                              #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Topo1si and RNaseH2si experiments                                   #
################################################################################


#----------------------------Paths---------------------------------------------#
projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
projPathCbw=$projPathC/alignment/bigwig
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
projPathSbw=$projPathS/alignment/bigwig

outpath=/project/ChromGroup_Seq_data/Rloop_paper/heatmaps
gencode=/project/ChromGroup_Seq_data/Celeste/analysis/gencode.v44.annotation.gtf

#-------------------------Script-----------------------------------------------#
bwsuf=bowtie2.fragments.normalized.filtered.bw
peaksuf=merged_seacr_norm_control.peaks.relaxed.filtered.rmblck.bed
balcklist=$projPathC/GRCh38_unified_blacklist.bed

samples=("$projPathSbw/Cntsi.Flag.m_$bwsuf" "$projPathSbw/Topo1si.fl_$bwsuf" \
"$projPathSbw/RNAseH1si.fl_$bwsuf"\  
"$projPathCbw/HBD.igg_$bwsuf")

computeMatrix reference-point --referencePoint TSS -S ${samples[@]} \
-R $gencode \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 5000 \
-b 5000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_genes.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_genes.mat.gz \
--colorList 'white,#000066' 'white,#ff7f50' 'white,#990000' 'white,#000000' \
--samplesLabel BRIGHTER2.CNTsi BRIGHTER2.Topo1si BRIGHTER2.RNAseH1si BRIGHTER1.IgG \
--boxAroundHeatmaps no \
-out $outpath/plot/${nameplot}_tss_heatmap_genes.png \
#

#Matrix over genes
computeMatrix scale-regions -S ${samples[@]} \
-R $gencode \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros -o $outpath/data/${namesummary}_genes.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_genes.mat.gz \
--samplesLabel BRIGHTER2.Topo1si BRIGHTER2.RNAseH1si BRIGHTER2.CNTsi BRIGHTER1.IgG \
-out $outpath/plot/${nameplot}_heatmap_genes.png \

namesummary=topo1si_fl_comparison
nameplot=topo1_fl_comparison

bedtools intersect -a $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf \
-b $projPathS/peakCalling/SEACR/Topo1si.fl_$peaksuf > $outpath/data/cntsi_and_Topo_peaks.bed


bedtools intersect -v -a $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf \
-b $projPathS/peakCalling/SEACR/Topo1si.fl_$peaksuf > $outpath/data/cntsi_only_peaks.bed

bedtools intersect -v -a $projPathS/peakCalling/SEACR/Topo1si.fl_$peaksuf \
-b $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf > $outpath/data/topo_only_peaks.bed

samples=("$projPathSbw/Cntsi.Flag.m_$bwsuf" "$projPathSbw/Topo1si.fl_$bwsuf" \
"$projPathCbw/HBD.igg_$bwsuf")


computeMatrix reference-point --referencePoint center -S ${samples[@]} \
-R $outpath/data/cntsi_only_peaks.bed \
$outpath/data/cntsi_and_Topo_peaks.bed \
$outpath/data/topo_only_peaks.bed \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 3000 \
-b 3000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_peaks.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_peaks.mat.gz \
--colorList 'white,#000066' 'white,#ff7f50' 'white,#000000' \
--samplesLabel BRIGHTER2_CNTsi BRIGHTER2_Topo1si BRIGHTER1_IgG \
--regionsLabel "B2_CNTsi peaks" "B2_CNTsi+B2_Topo1si peaks" "B2_Topo1si peaks" \
--heatmapWidth 5 \
--heatmapHeight 30 \
--boxAroundHeatmaps no \
-out $outpath/plot/${nameplot}_center_heatmap_peaks.png \




namesummary=rnaseh1si_fl_comparison
nameplot=rnaseh1si_fl_comparison

bedtools intersect -a $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf \
-b $projPathS/peakCalling/SEACR/RNAseH1si.fl_$peaksuf > $outpath/data/cntsi_and_rnaseh1_peaks.bed


bedtools intersect -v -a $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf \
-b $projPathS/peakCalling/SEACR/RNAseH1si.fl_$peaksuf > $outpath/data/cntsi_notrnaseh1_peaks.bed

bedtools intersect -v -a $projPathS/peakCalling/SEACR/RNAseH1si.fl_$peaksuf \
-b $projPathS/peakCalling/SEACR/Cntsi.Flag.m_$peaksuf > $outpath/data/rnase_only_peaks.bed

samples=("$projPathSbw/Cntsi.Flag.m_$bwsuf" \
"$projPathSbw/RNAseH1si.fl_$bwsuf" \  
"$projPathCbw/HBD.igg_$bwsuf")


computeMatrix reference-point --referencePoint center -S ${samples[@]} \
-R $outpath/data/cntsi_notrnaseh1_peaks.bed \
$outpath/data/cntsi_and_rnaseh1_peaks.bed \
$outpath/data/rnase_only_peaks.bed \
--blackListFileName $projPathC/GRCh38_unified_blacklist.bed \
-a 3000 \
-b 3000 \
--missingDataAsZero \
--skipZeros -o $outpath/data/${namesummary}_refpoint_peaks.mat.gz

plotHeatmap -m $outpath/data/${namesummary}_refpoint_peaks.mat.gz \
--colorList 'white,#000066' 'white,#990000' 'white,#000000' \
--samplesLabel BRIGHTER2_CNTsi BRIGHTER2_RNAseH1si BRIGHTER1_IgG \
--regionsLabel "B2_CNTsi peaks" "B2_CNTsi+B2_RNAseH1si peaks" "B2_RNAseH1 peaks" \
--heatmapWidth 5 \
--heatmapHeight 30 \
--boxAroundHeatmaps no \
-out $outpath/plot/${nameplot}_center_heatmap_peaks.png \
