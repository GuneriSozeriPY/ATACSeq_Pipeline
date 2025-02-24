#!/bin/bash

#SBATCH -J atac_lengthQC

#SBATCH --account=investor

#SBATCH --qos=mid_investor

#SBATCH --partition=mid_investor
#SBATCH --cpus-per-task=10

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/QC/Slurm/lengthQC-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/QC/Slurm/lengthQC-%j.err

#SBATCH --array=0-6

#conda activate atac to activate picard

module load picard-2.20.8-gcc-9.2.0-hehlgdl

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/QC"
bam_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample $sample , QC starts"

picard CollectInsertSizeMetrics \
 I="$bam_path/${sample}.blacklist-filtered.bam" \
 O="$output_path/${sample}.blacklist-filtered.fraglen.stats" \
 H="$output_path/${sample}.blacklist-filtered.fraglen.pdf" M=0.5

echo "Done"