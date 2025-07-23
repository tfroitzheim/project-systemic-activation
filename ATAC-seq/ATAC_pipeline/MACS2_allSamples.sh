singularity exec /n/holylabs/LABS/whited_lab/Lab/singularity/macs2:2.2.9.1--py39hf95cd2a_0 macs2 callpeak -t bams/BLA-1_uniq.rmdup.bam bams/BLA-2_uniq.rmdup.bam bams/BLA-3_uniq.rmdup.bam bams/BLA-4_uniq.rmdup.bam bams/BLA-5_uniq.rmdup.bam bams/CTL-1_uniq.rmdup.bam bams/CTL-2_uniq.rmdup.bam bams/CTL-3_uniq.rmdup.bam bams/CTL-4_uniq.rmdup.bam bams/CTL-5_uniq.rmdup.bam bams/INT-1_uniq.rmdup.bam bams/INT-2_uniq.rmdup.bam bams/INT-3_uniq.rmdup.bam bams/INT-4_uniq.rmdup.bam bams/INT-5_uniq.rmdup.bam \
  -g 3.2e+10 \
  --nomodel -n allsamppeaks \
  &> macs2.log 