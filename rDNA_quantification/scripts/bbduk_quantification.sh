#RDNA quantification of silencing experiments in 

bbduk=/home/lopez_s/bbmap/bbduk.sh
ref_rDNA=/project/ChromGroup_Seq_data/Rloop_paper/rDNA_quantification/data/Homo_sapiens_rRNA.fa
RAWROOT=/project/avitidata/pipelined/20240917_AV234501_A-PE75-Ad-FS-HO/Samples/V
outpath=/project/ChromGroup_Seq_data/Rloop_paper/rDNA_quantification/ 
  

bash $bbduk in1=$RAWROOT/mpimg_L31496_POOL-SM-1-lane1-007-HBD-Flag-mo-r1/mpimg_L31496_POOL-SM-1-lane1-007-HBD-Flag-mo-r1_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31496_POOL-SM-1-lane1-007-HBD-Flag-mo-r1/mpimg_L31496_POOL-SM-1-lane1-007-HBD-Flag-mo-r1_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Flag.mo_1_matched.fastq stats=$outpath/data/Flag.mo_1_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31496_POOL-SM-1-lane1-008-HBD-Flag-mo-r2/mpimg_L31496_POOL-SM-1-lane1-008-HBD-Flag-mo-r2_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31496_POOL-SM-1-lane1-008-HBD-Flag-mo-r2/mpimg_L31496_POOL-SM-1-lane1-008-HBD-Flag-mo-r2_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Flag.mo_2_matched.fastq stats=$outpath/data/Flag.mo_2_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-001-Cntsi-flag-m-r1/mpimg_L31497_POOL-SM-2-lane1-001-Cntsi-flag-m-r1_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-001-Cntsi-flag-m-r1/mpimg_L31497_POOL-SM-2-lane1-001-Cntsi-flag-m-r1_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Cntsi.Flag.m_1_matched.fastq stats=$outpath/data/Cntsi.Flag.m_1_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-002-Cntsi-flag-m-r2/mpimg_L31497_POOL-SM-2-lane1-002-Cntsi-flag-m-r2_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-002-Cntsi-flag-m-r2/mpimg_L31497_POOL-SM-2-lane1-002-Cntsi-flag-m-r2_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Cntsi.Flag.m_2_matched.fastq stats=$outpath/data/Cntsi.Flag.m_2_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-005-Topo1si-fl-r1/mpimg_L31497_POOL-SM-2-lane1-005-Topo1si-fl-r1_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-005-Topo1si-fl-r1/mpimg_L31497_POOL-SM-2-lane1-005-Topo1si-fl-r1_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Topo1si.fl_1_matched.fastq stats=$outpath/data/Topo1si.fl_1_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-006-Topo1si-fl-r2/mpimg_L31497_POOL-SM-2-lane1-006-Topo1si-fl-r2_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-006-Topo1si-fl-r2/mpimg_L31497_POOL-SM-2-lane1-006-Topo1si-fl-r2_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/Topo1si.fl_2_matched.fastq stats=$outpath/data/Topo1si.fl_2_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-003-RNAseH1si-fl-r1/mpimg_L31497_POOL-SM-2-lane1-003-RNAseH1si-fl-r1_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-003-RNAseH1si-fl-r1/mpimg_L31497_POOL-SM-2-lane1-003-RNAseH1si-fl-r1_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/RNAseH1si.fl_1_matched.fastq  stats=$outpath/data/RNAseH1si.fl_1_matched_stats.txt

bash $bbduk in1=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-004-RNAseH1si-fl-r2/mpimg_L31497_POOL-SM-2-lane1-004-RNAseH1si-fl-r2_R1.fastq.gz \
in2=$RAWROOT/mpimg_L31497_POOL-SM-2-lane1-004-RNAseH1si-fl-r2/mpimg_L31497_POOL-SM-2-lane1-004-RNAseH1si-fl-r2_R2.fastq.gz \
ref=$ref_rDNA k=31 hdist=2 outm=$outpath/data/RNAseH1si.fl_2_matched.fastq stats=$outpath/data/RNAseH1si.fl_2_matched_stats.txt
