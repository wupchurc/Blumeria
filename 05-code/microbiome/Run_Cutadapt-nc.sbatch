#!/bin/bash -l
#SBATCH --time=4:00:00
#SBATCH --ntasks=8
#SBATCH --mem=120g
#SBATCH --tmp=40g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wupchurc@umn.edu

cd /scratch.global/wupchurc/Blumeria-ITS2 || exit 1
module load cutadapt fastqc parallel || { echo "Failed to load modules"; exit 1; }

# 1. Setup Environment
mkdir -p 01_adapter-nc 02_filtered-nc 03_post_qc-nc
: > run_trimming.sh

# 2. Define Sequences for UMGC Methods
ADAPTER="CTGTCTCTTATACACATCT"		# Nextera Adapter
FWD_PRIMER="TCGATGAAGAACGCAGCG"		# 5.8SR (Forward Primer)
REV_PRIMER="TCCTCCGCTTATTGATATGC"	# ITS4 (Reverse)
RC_FWD="CGCTGCGTTCTTCATCGA"       
RC_REV="GCATATCAATAAGCGGAGGA"

# 3. Build the Command List
for i in *_R1_001.fastq.gz; do
	R1=$i
	R2=${i/_R1_/_R2_}
	SAMPLE=${i%_R1_001.fastq.gz}

	cat <<EOF >> run_trimming.sh
	cutadapt --cores 1 \
	-a $ADAPTER -A $ADAPTER \
	-o 01_adapter-nc/${R1} -p 01_adapter-nc/${R2} \
	${R1} ${R2} > 01_adapter-nc/${SAMPLE}.adapter.log && \
	cutadapt --cores 1 \
	-g $FWD_PRIMER -a $RC_REV \
	-G $REV_PRIMER -A $RC_FWD \
	-n 2 --discard-untrimmed --pair-filter=any --minimum-length 50 \
	-o 02_filtered-nc/${R1} -p 02_filtered-nc/${R2} \
	01_adapter-nc/${R1} 01_adapter-nc/${R2} > 02_filtered-nc/${SAMPLE}.primer.log
EOF
done

# 4. Execute
parallel -j 8 < run_trimming.sh
fastqc -t 8 02_filtered-nc/*.fastq.gz -o 03_post_qc-nc/

