#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu

#raw data directory
cd /scratch.global/wupchurc/ITS2

module load cutadapt parallel || { echo "Load cutadapt/parallel modules"; exit 1; }

# Step 1: Remove adapters -> 01_adapter/
mkdir -p 01_adapter
for i in *_R1_001.fastq.gz; do
	echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 \
	       	-a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA \
		-g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT \
		-A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT \
		-o 01_adapter/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} \
		-p 01_adapter/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} \
		${i} ${i//_R1_/_R2_} \
       		> 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" \
		>> run_cutadapt_adapters.sh
done
chmod +x run_cutadapt_adapters.sh
parallel < run_cutadapt_adapters.sh

cd 01_adapter
mkdir -p 01_logs
mv cutadapt*.log.txt 01_logs/ 2>/dev/null || true
grep "passing" 01_logs/*.log.txt > summary_adapter_trimming.txt || true
echo "Adapter trimming summaries in summary_adapter_trimming.txt"

# Step 2: Remove ITS2 primers -> 02_filtered/
mkdir -p ../02_filtered
for i in *_R1_001.fastq.gz; do
       	echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed \
	       	-g TCGATGAAGAACGCAGCG -G TCCTCCGCTTATTGATATGC \
	       	-o ../02_filtered/${i} -p ../02_filtered/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} \
	       	${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.primer.log.txt}" \
	       	>> run_cutadapt_primers.sh
done
chmod +x run_cutadapt_primers.sh
parallel < run_cutadapt_primers.sh

cd ../02_filtered/
mkdir -p 02_logs
mv *.primer.log.txt 02_logs/ 2>/dev/null || true  # Renamed logs for clarity
grep "passing" 02_logs/*.log.txt > summary_primer_trimming.txt || true
echo "Primer trimming summaries in summary_primer_trimming.txt"
echo "Cleaned paired-end files ready in 02_filtered/ for DADA2 or QIIME2."
