#!/bin/bash

#SBATCH -J FastQC

#SBATCH --account=investor

#SBATCH --qos=long_investor

#SBATCH --partition=long_investor

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc/Slurm/fastqc_afterFastpadap-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc/Slurm/fastqc_afterFastpadap-%j.err

#SBATCH --array=0-6

#conda activate fastp

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/adapter"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/adapter"
fastqc_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/fastqc/afterFastp"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample $sample"    

fastqc --t 32 $output_path/${sample}_R1.fastq.gz --outdir $fastqc_path
fastqc --t 32 $output_path/${sample}_R2.fastq.gz --outdir $fastqc_path

echo "Fastqc after trimming $sample done"


