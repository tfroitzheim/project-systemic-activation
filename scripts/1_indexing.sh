#!/bin/bash
#SBATCH --job-name=kb_indexing
#SBATCH -c 1
#SBATCH -t 1-12:00
#SBATCH -p sapphire
#SBATCH --mem=800G
#SBATCH -o ../log/1_index_%j.out
#SBATCH -e ../log/1_index_%j.err




# This script produces a reference index, a transcript-to-gene table, and a spliced transcript FASTA (All available on the Harvard dataverse)
# Inputs are a) a GTF file, and b) the reference genome sequence (also, both available on the Harvard dataverse)
# !! note that the GTF file used here was modified by stripping the gene feature line. 

# Activate environment
# Note: create "scanpy_env" using ./bin/scanpy_env.yml
module load miniconda3/py39_4.11.0-linux_x64-ncf
source activate scanpy_env

# some users may need to run 'kb compile all', to allow Kallisto to run without errors:
       #module load zlib
       #module load autoconf
       #module load cmake
       #kb compile all

# Run the kb ref command 
kb ref -i ../output/kb_ref/index.idx -g ../output/kb_ref/t2g.tsv -f1 ../output/kb_ref/cDNA.fa \
       ../ref/GCF_040938575.1_UKY_AmexF1_1_genomic.fa \
       ../ref/GCF_040938575.1_UKY_AmexF1_1_genomic.noGene.gtf