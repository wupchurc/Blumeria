#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32  # 4 jobs * 8 cores each
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu
#SBATCH --job-name=meta_test

# raw data directory
cd /scratch.global/wupchurc/Metagenomics || { echo "Directory not found"; exit 1; }

module load cutadapt bowtie2/2.4.4 parallel || { echo "Failed to load modules"; exit 1; }

: > run_bt2.cmd
: > run_cutadapt.cmd

#make log directory
mkdir -p logs

# remove rodent host sequences
mkdir -p 01_host_removed
: > run_bt2.cmd

for i in *_R1.fastq.gz; do
	echo "bowtie2 -p 8 -x /common/bioref/ensembl/main/Rattus_norvegicus-115/GRCr8/bowtie2/genome -1 \"$i\" -2 \"${i/_R1_/_R2_}\" --un-conc-gz 01_host_removed/${i/_R1_/_host_removed} > logs/${i%_R1.fastq.gz}.bowtie2.log 2>&1" >> run_bt2.cmd
done

parallel -j 4 < run_bt2.cmd

mkdir -p 02_trimmed
for i in 01_host_removed/*_R1.fastq.1.gz; do
	r2="${i/1.gz/2.gz}"
	sample=$(basename "${i%_R1.fastq.1.gz}")
	echo  "cutadapt --cores 8 --minimum-length 50 -a TACACTCTTTCCCTACACGACGCTCTTCCGATCT -A GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o 02_trimmed/${sample}_R1.trimmed.fastq.gz -p 02_trimmed/${sample}_R2.trimmed.fastq.gz \"$i\" \"$r2\" > logs/${sample}.cutadapt.log 2>&1" >> run_cutadapt.cmd
done

parallel -j 4 < run_cutadapt.cmd

echo "Pipeline completed successfully"
