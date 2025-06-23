#!/bin/bash
#SBATCH -c 1
#SBATCH -t 1-12:00
#SBATCH -p sapphire
#SBATCH --mem=800G
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err



# This script produces a reference index, a transcript-to-gene table, and a spliced transcript FASTA (All available on the Harvard dataverse)
# Inputs are a) a GTF file, and b) the reference genome sequence (also, both available on the Harvard dataverse)
# !! note that the GTF file used here was modified by stripping the gene feature line. 

# Activate environment
source activate scanpy_env

# Run the kb ref command 
kb ref -i ../output/kb_ref/index.idx -g ../output/kb_ref/t2g.tsv -f1 ../output/kb_ref/cDNA.fa \
       ../ref/GCF_040938575.1_UKY_AmexF1_1_genomic.fa \
       ../ref/GCF_040938575.1_UKY_AmexF1_1_genomic.noGene.gtf