#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32  # 4 jobs * 8 cores each
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu
#SBATCH --job-name=blum_meta_pre

# Change to raw data directory
cd /scratch.global/wupchurc/Blumeria-Metagenomics || { echo "Directory not found"; exit 1; }

# Load modules
module load cutadapt bowtie2/2.4.4 parallel fastqc || { echo "Failed to load modules"; exit 1; }

# Clean/create output and log architectures
mkdir -p logs 01_host_removed 02_trimmed 03_fastqc
: > run_bt2.cmd
: > run_cutadapt.cmd
: > run_fastqc.cmd

# ---------------------------------------------------------------------
# Step 1: Remove rodent host sequences
# ---------------------------------------------------------------------
# Assumes your files are strictly named SampleName_R1.fastq.gz
for i in *_R1.fastq.gz; do
	# Safely isolate the long base name by stripping the trailing suffix
	sample="${i%_R1.fastq.gz}"
	r2="${sample}_R2.fastq.gz"

	# Build the command. Note: Bowtie2 automatically appends .1.gz and .2.gz
	# to the --un-conc-gz path for paired-end unmapped reads.
	echo "bowtie2 -p 8 -x /common/bioref/ensembl/main/Rattus_norvegicus-115/GRCr8/bowtie2/genome -1 \"$i\" -2 \"$r2\" \
	       	--un-conc-gz 01_host_removed/${sample}_host_removed.fastq.gz -S /dev/null > logs/${sample}.bowtie2.log 2>&1" >> run_bt2.cmd
done 

parallel -j 4 < run_bt2.cmd

# ---------------------------------------------------------------------
# Step 2: Trim Adapters
# ---------------------------------------------------------------------
# Bowtie2 automatically inserts ".1" and ".2" before ".gz" for paired-end unmapped files
for i in 01_host_removed/*_host_removed.fastq.1.gz; do
	r2="${i/.1.gz/.2.gz}"
	sample=$(basename "${i%_host_removed.fastq.1.gz}")

	echo "cutadapt --cores 8 --minimum-length 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	       -o 02_trimmed/${sample}_R1.trimmed.fastq.gz -p 02_trimmed/${sample}_R2.trimmed.fastq.gz \"$i\" \"$r2\" > logs/${sample}.cutadapt.log 2>&1" \
	       >> run_cutadapt.cmd
done

parallel -j 4 < run_cutadapt.cmd

# ---------------------------------------------------------------------
# Step 3: FastQC Quality Check
# ---------------------------------------------------------------------
# Iterate over the newly trimmed R1 files
for i in 02_trimmed/*_R1.trimmed.fastq.gz; do
        sample=$(basename "${i%_R1.trimmed.fastq.gz}")
        r2="02_trimmed/${sample}_R2.trimmed.fastq.gz"

	# -t 8 matches your 8 threads per job. -o specifies the output directory.
        # FastQC handles both R1 and R2 in the same command.
        echo "fastqc -t 8 -o 03_fastqc \"$i\" \"$r2\" > logs/${sample}.fastqc.log 2>&1" >> run_fastqc.cmd
done

parallel -j 4 < run_fastqc.cmd

echo "Pipeline completed successfully"
