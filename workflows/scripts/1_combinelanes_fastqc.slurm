#!/bin/bash

#SBATCH -J fastqc

#SBATCH --account=investor

#SBATCH --qos=long_investor

#SBATCH --partition=long_investor

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc/Slurm/combinelanes_fastqc-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc/Slurm/combinelanes_fastqc-%j.err
#SBATCH --array=0-6

#conda activate fastp

sample_path="/cta/groups/adebali/ATAC_part2/proj/sancarlb/HTSF/250109_UNC52-VH01297_223_22277GHNX"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
fastqc_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc"

sample_list=($(find "$sample_path" -maxdepth 1 -name "*_L001_R1_001.fastq.gz" | grep "NHF1" | cut -d "/" -f10 | cut -d "_" -f1-6))
sample=${sample_list[$SLURM_ARRAY_TASK_ID]}

echo $SLURM_ARRAY_TASK_ID
echo $sample

# Combine lanes
zcat $sample_path/${sample}_L001_R1_001.fastq.gz $sample_path/${sample}_L002_R1_001.fastq.gz > $output_path/${sample}_R1.fastq
zcat $sample_path/${sample}_L001_R2_001.fastq.gz $sample_path/${sample}_L002_R2_001.fastq.gz > $output_path/${sample}_R2.fastq

sample_list2=($(find "$sample_path" -maxdepth 1 -name "*_L001_R1_001.fastq.gz" | grep "NHF1" | cut -d "/" -f10 | cut -d "_" -f1-4))
newsample=${sample_list2[$SLURM_ARRAY_TASK_ID]}
echo $newsample

# Rename files
mv $output_path/${sample}_R1.fastq $output_path/${newsample}_R1.fastq
mv $output_path/${sample}_R2.fastq $output_path/${newsample}_R2.fastq

# Run FastQC on the renamed files
fastqc -t 32 $output_path/${newsample}_R1.fastq --outdir $fastqc_path
fastqc -t 32 $output_path/${newsample}_R2.fastq --outdir $fastqc_path




