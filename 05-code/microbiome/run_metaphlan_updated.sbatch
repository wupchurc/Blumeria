#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu
#SBATCH --job-name=blu-mphlan


cd /scratch.global/wupchurc/Blumeria-Metagenomics || { echo "Directory not found"; exit 1; }

module load metaphlan/4.0.6 parallel || { echo "Failed to load modules"; exit 1; }

mkdir -p 04_metaphlan logs
: > metaphlan.cmd

for r1 in 02_trimmed/*_R1.trimmed.fastq.gz; do
    r2="${r1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}"
    sample=$(basename "${r1%_R1.trimmed.fastq.gz}")


    echo "metaphlan \"$r1,$r2\" --input_type fastq --bowtie2db /scratch.global/wupchurc/metaphlan_db \
	    --nproc 8 --bowtie2out 04_metaphlan/${sample}.bowtie2out.bz2 -o 04_metaphlan/${sample}_profile.txt \
	    > logs/${sample}.metaphlan.log 2>&1" >> metaphlan.cmd
done

parallel -j 4 < metaphlan.cmd

# ---------------------------------------------------------------------
# Step 2: Merge Abundance Tables
# ---------------------------------------------------------------------
# Define the final results directory
RESULTS_DIR="/projects/standard/szprisco/shared/wupchurc_projects/2026/Rat_Stool_Microbiome_Blumeria/04-results"

# Create the directory just in case it doesn't exist yet
mkdir -p "$RESULTS_DIR"

# Merge all profiles and save directly to the final directory
merge_metaphlan_tables.py 04_metaphlan/*_profile.txt > "$RESULTS_DIR"/metagenomics_merged_abundance_table.txt

echo "MetaPhlAn pipeline completed successfully"
