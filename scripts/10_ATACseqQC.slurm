#!/bin/bash

#SBATCH -J ATACseqQC

#SBATCH --account=investor

#SBATCH --qos=long_investor

#SBATCH --partition=long_investor
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/Slurm/ATACseqQC_0-6_rep2-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/Slurm/ATACseqQC-0-6_rep2%j.err

#SBATCH --array=0-6

module load r-4.1.0-gcc-9.2.0-wq6f22q

export TMPDIR=/cta/users/ygsuzeri/TMP1  

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}


echo "Processing sample $sample , ATACseqQC starts"
#run bowtie2
Rscript /cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/atacQC_NHF1_rep2.R $sample

echo "Done"