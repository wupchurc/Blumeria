#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu

cd /scratch.global/wupchurc/ITS2
module load cutadapt fastqc parallel || { echo "Failed to load modules"; exit 1; }

# 1. Setup Environment
mkdir -p 01_adapter 02_filtered 03_post_qc
: > run_trimming.sh	# clear file

# 2. Define Sequences for UMGC Methods
ADAPTER="CTGTCTCTTATACACATCT"		# Nextera Adapter
FWD_PRIMER="TCGATGAAGAACGCAGCG"		# 5.8SR (Forward Primer)
REV_PRIMER="TCCTCCGCTTATTGATATGC"	# ITS4 (Reverse)

# 3. Build the Command List
for i in *_R1_001.fastq.gz; do
	R1=$i
	R2=${i/_R1_/_R2_}
	SAMPLE=${i%_R1_001.fastq.gz}

	# Remove Nextera Adapters first, then trim primers and discard non-fungal reads
	echo "cutadapt --cores 1 -a $ADAPTER -A $ADAPTER -o 01_adapter/${R1} -p 01_adapter/${R2} ${R1} ${R2} > 01_adapter/${SAMPLE}.adapter.log && \
	cutadapt --cores 1 -g $FWD_PRIMER -G $REV_PRIMER --discard-untrimmed --pair-filter=any --minimum-length 50 -o 02_filtered/${R1} -p 02_filtered/${R2} 01_adapter/${R1} 01_adapter/${R2} > 02_filtered/${SAMPLE}.primer.log" >> run_trimming.sh
done

# 4. Execute
parallel -j 8 < run_trimming.sh
fastqc -t 8 02_filtered/*.fastq.gz -o 03_post_qc/

