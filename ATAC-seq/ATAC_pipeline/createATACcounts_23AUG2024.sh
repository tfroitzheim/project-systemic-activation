#!/bin/bash

# Define common variables
SAF_FILE="/n/home09/sjblair/allsamppeaks_peaks.narrowPeak.saf"
SINGULARITY_IMAGE="/n/holylabs/LABS/whited_lab/Lab/singularity/subread_v2.0.3.sing"
BASE_DIR="/n/holyscratch01/whited_lab/sjblair/ATAC_MAR5_2024/BAM"
OUTPUT_DIR="/n/holyscratch01/whited_lab/sjblair/ATAC_MAR5_2024/MACS2_COUNTS"

# Define prefixes and corresponding output file prefixes
prefixes=("BLA" "CTL" "INT")
output_prefixes=("bl" "cl" "in")

# Iterate over prefixes
for i in "${!prefixes[@]}"; do
  prefix="${prefixes[i]}"
  output_prefix="${output_prefixes[i]}"

  # Iterate over numbers 1 to 5
  for num in {1..5}; do
    bam_file="$BASE_DIR/$prefix-$num/$prefix-$num"_uniq.rmdup.bam
    output_file="$OUTPUT_DIR/$output_prefix$num.atac.counts"

    ln -s "$bam_file"
    singularity exec "$SINGULARITY_IMAGE" featureCounts -p -F SAF -a "$SAF_FILE" --fracOverlap 0.2 -o "$output_file" "$bam_file"
  done
done