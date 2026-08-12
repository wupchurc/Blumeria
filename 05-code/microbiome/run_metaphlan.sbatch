#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu
#SBATCH --job-name=metaphlan
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

cd /scratch.global/wupchurc/Metagenomics || exit 1

module load metaphlan/4.0.6 parallel

mkdir -p 03_metaphlan logs
: > metaphlan.cmd

for r1 in 02_trimmed/*_R1.trimmed.fastq.gz; do
    r2="${r1/_R1.trimmed.fastq.gz/_R2.trimmed.fastq.gz}"
    sample=$(basename "${r1%_R1.trimmed.fastq.gz}")
    echo "metaphlan \"$r1\",\"$r2\" --input_type fastq --bowtie2db /scratch.global/wupchurc/metaphlan_db --nproc 8 --bowtie2out 03_metaphlan/${sample}.bowtie2out.bz2 -o 03_metaphlan/${sample}_profile.txt > logs/${sample}.metaphlan.log 2>&1" >> metaphlan.cmd
done

parallel -j 4 < metaphlan.cmd
