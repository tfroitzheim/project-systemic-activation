#!/bin/bash
#
#SBATCH -p shared 
#SBATCH -c 36 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH -t 10:00:00 # time (D-HH:MM)
#SBATCH -o CTL-ATAC.%N.%j.out # STDOUT
#SBATCH -e CTL-ATAC.%N.%j.err # STDERR

export FQ=/n/whited_lab/Lab/data/230413_A00794_0824_BHWT3NDSX5_SUB14089_singleindex/fastq
export INDEX=/n/whited_lab/Everyone/resources/AmexG_v6.0-DD.chunked.bowtie2/ambMex60DD.chunked.bowtie2index
export CHROMESIZE=/n/whited_lab/Everyone/resources/AmexG_v6.0-DD.chunked.bowtie2/chrom.sizes.tab
export THREADS=32
module load gcc/12.1.0-fasrc01
module load bowtie2/2.3.2-fasrc02
module load samtools/1.10-fasrc01
module load jdk/1.8.0_172-fasrc01
module load ucsc/20150820-fasrc01

for SAMPLE in CTL-1 CTL-2 CTL-3 CTL-4 CTL-5; do
	export BAM=/n/holyscratch01/whited_lab/sjblair/SAC4NATAC/noUnplaced_BAM/$SAMPLE
	mkdir $BAM
	bowtie2 -q -p $THREADS --no-mixed -X 2000 --dovetail --no-discordant -x $INDEX -1 $FQ/$SAMPLE'_R1.trimmed.fq.gz' -2 $FQ/$SAMPLE'_R2.trimmed.fq.gz' | samtools view -bSu - > $BAM/$SAMPLE'.unsorted.bam'
	samtools sort -o $BAM/$SAMPLE'.bam' $BAM/$SAMPLE'.unsorted.bam'
	samtools index -b $BAM/$SAMPLE'.bam'
	rm $BAM/$SAMPLE'.unsorted.bam'
	# Filter aligned BAM
	samtools sort -@ $THREADS -o $BAM/$SAMPLE'.total.bam' $BAM/$SAMPLE'.bam'
	samtools index -b -@ $THREADS $BAM/$SAMPLE'.total.bam'
	samtools view -@ $THREADS -h $BAM/$SAMPLE'.total.bam' | grep -v chrM | samtools sort -O bam -o $BAM/$SAMPLE'.rmChrM.bam'
	samtools view -@ $THREADS -h -q 30 $BAM/$SAMPLE'.rmChrM.bam' > $BAM/$SAMPLE'.rmChrm.unsorted.bam'
	samtools sort -@ $THREADS -o $BAM/$SAMPLE'_uniq.bam' $BAM/$SAMPLE'.rmChrm.unsorted.bam'
	samtools index -b -@ $THREADS $BAM/$SAMPLE'_uniq.bam'
	rm $BAM/$SAMPLE'.rmChrm.unsorted.bam'
	$PICARD picard MarkDuplicates INPUT=$BAM/$SAMPLE'_uniq.bam' OUTPUT=$BAM/$SAMPLE'_uniq.rmdup.unsorted.bam' METRICS_FILE=$BAM/$SAMPLE'_picard.rmDup.txt' ASSUME_SORTED=true REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.2
	samtools sort -n -@ $THREADS -o $BAM/$SAMPLE'_uniq.rmdup.bam' $BAM/$SAMPLE'_uniq.rmdup.unsorted.bam'
	samtools sort -@ $THREADS -o $BAM/$SAMPLE'_uniq.rmdup_forigv.bam' $BAM/$SAMPLE'_uniq.rmdup.unsorted.bam'
	samtools index -b -@ $THREADS $BAM/$SAMPLE'_uniq.rmdup_forigv.bam'
	rm $BAM/$SAMPLE'_uniq.rmdup.unsorted.bam'
	$PICARD CollectInsertSizeMetrics HISTOGRAM_FILE=$BAM/$SAMPLE'_fragsize.pdf' OUTPUT=$BAM/$SAMPLE'_frag_size.txt' METRIC_ACCUMULATION_LEVEL=ALL_READS INCLUDE_DUPLICATES=false INPUT=$BAM/$SAMPLE'_uniq.rmdup.bam'
	samtools view -@ $THREADS -c $BAM/$SAMPLE'.total.bam' > $BAM/$SAMPLE'_cnt_trimmed.txt'
	samtools view -@ $THREADS -c -F 4 $BAM/$SAMPLE'.total.bam'> $BAM/$SAMPLE'_cnt_mapped.txt'
	samtools view -@ $THREADS -c $BAM/$SAMPLE'.rmChrM.bam' > $BAM/$SAMPLE'_cnt_rmChrM.txt'
	samtools view -@ $THREADS -c $BAM/$SAMPLE'_uniq.bam' > $BAM/$SAMPLE'_cnt_uniq.txt'
	samtools view -@ $THREADS -c $BAM/$SAMPLE'_uniq.rmdup.bam' > $BAM/$SAMPLE'_cnt_uniq.rmdup.txt'
	# samtools sort -n -@ $THREADS -o $BAM/$SAMPLE'.sorted.bam' $BAM/$SAMPLE'_uniq.rmdup.bam'
	$DEEPTOOLS bamCoverage -b $BAM/$SAMPLE'.total.bam' -o $BAM/$SAMPLE'.bw' --outFileFormat=bigwig --normalizeUsing CPM -p 16 --binSize 20
done