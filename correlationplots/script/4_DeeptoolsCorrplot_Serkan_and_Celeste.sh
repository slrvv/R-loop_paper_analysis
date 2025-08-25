################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# Reproducibility assessment                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create deeptools heatmap for reproducibility assessment for         #
# Celeste"s & Serkan's data. Compare all of their data                         #
################################################################################

#-----------------------------Paths--------------------------------------------#

projPathC=/project/ChromGroup_Seq_data/Celeste/2024_summer_complete
experimentSummaryC=$projPathC/experiment_summary_align_formatted.csv
projPathS=/project/ChromGroup_Seq_data/Serkan/2024_09
experimentSummaryS=$projPathS/experiment_summary_alignment.csv
namesummary=all_experiments_bam_Celeste_Serkan
nameplot=all_experiments_corr_heatmap_Celeste_Serkan
outPath=/project/ChromGroup_Seq_data/Rloop_paper
#-------------------------Script-----------------------------------------------#

#populate bam files array
bamfiles=()
samples=()
while IFS=, read -r SampleName R1 R2
do
   bam=$projPathC/alignment/bam/$SampleName.mapped.bam.sorted
   bamfiles+=" $bam"
   samples+=" ${SampleName}_C"
done < <(tail -n +2 $experimentSummaryC)

while IFS=, read -r SampleName R1 R2
do
   bam=$projPathS/alignment/bam/$SampleName.mapped.bam.sorted
   bamfiles+=" $bam"
   samples+=" ${SampleName}_S"
done < <(tail -n +2 $experimentSummaryS)

echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples \
  -o $outPath/$namesummary.npz

plotCorrelation -in $outPath/$namesummary.npz -c pearson  \
--colorMap RdYlBu --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/$nameplot.pdf \
--outFileCorMatrix $outPath/$nameplot.cormatrix.txt\

#populate bam files array
bamfiles=()
samples=()
while IFS=, read -r SampleName R1 R2
do
   bam=$projPathC/alignment/bam/$SampleName.mapped.bam.sorted
   bamfiles+=" $bam"
   samples+=" ${SampleName}_C"
done < <(tail -n +2 $experimentSummaryC)

while IFS=, read -r SampleName R1 R2
do
   bam=$projPathS/alignment/bam/$SampleName.mapped.bam.sorted
   bamfiles+=" $bam"
   samples+=" ${SampleName}_S"
done < <(tail -n +2 $experimentSummaryS)

echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 5000 \
  -o $outPath/${namesummary}_5000bp.npz

plotCorrelation -in $outPath/${namesummary}_5000bp.npz -c pearson  \
--colorMap RdYlBu --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/${nameplot}_5000bp.pdf \
--outFileCorMatrix $outPath/${nameplot}_5000bp.cormatrix.txt\


echo $samples
multiBamSummary bins --bamfiles $bamfiles -l $samples --binSize 2000 \
  -o $outPath/${namesummary}_2000bp.npz

plotCorrelation -in $outPath/${namesummary}_2000bp.npz -c pearson  \
--colorMap RdYlBu --plotNumbers --removeOutliers --plotHeight 18 --plotWidth  20 \
--plotFileFormat "pdf" -p heatmap -o $outPath/${nameplot}_2000bp.pdf \
--outFileCorMatrix $outPath/${nameplot}_2000bp.cormatrix.txt\