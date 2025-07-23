#!/bin/bash
SING=/n/holylabs/LABS/whited_lab/Lab/singularity # location of singularity library 
BAM=/n/holyscratch01/whited_lab/sjblair/ATAC_MAR5_2024/MACS2_input # Folder containing BAM files
BASE_OUT=/n/holyscratch01/whited_lab/sjblair/ATAC_MAR5_2024/MACS2_output_PE # The Base path to output files
# As of Mar 29 we have MACS2 and MACS3.  To see what versions are available check 'ls $SING/macs*'
MACS=macs2:2.2.9.1--py39hf95cd2a_0 # version of MACS2
GENOME_SIZE=3.2e+10 # Genome size
QVALUE=0.01 # set minimum FDR, 0.01 was used in Wei, 2021
THREADS=5 # 5 jobs will do an entire condition at once

for BAM_FILE in "$BAM"/*.bam; do
# Loop through all files ending with ".bam" in the bam directory
  # Extract sample name from filename
  SAMPLE=${BAM_FILE##*/}
  SAMPLE=${SAMPLE%%_*}

  # Output directory for the current sample
  OUT="$BASE_OUT/${SAMPLE}"

  # Create output directory if it doesn't exist
  if [ ! -d "$OUT" ]; then
    mkdir -p "$OUT"
    echo "Folder '$OUT' created successfully."
  fi

  # Run MACS2 peak calling for the current sample
  singularity exec $SING/$MACS macs2 callpeak -t "$BAM_FILE" \
  -f BAMPE -g $GENOME_SIZE -q $QVALUE \
  --nomodel -n $SAMPLE --outdir $OUT \
  &> "$OUT"/$SAMPLE.macs2.log 

  echo "Finished processing sample: $SAMPLE"
done