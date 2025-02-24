#!/bin/bash

#SBATCH -J ATAC_trimming

#SBATCH --account=investor

#SBATCH --qos=short_investor

#SBATCH --cpus-per-task=10

#SBATCH --partition=short_investor

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/Slurm/Fastp_dedupadapter-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/Slurm/Fastp_deduppadapter-%j.err

#SBATCH --array=0-6

#conda activate fastp

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/adapter"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample $sample"    

fastp --thread 16 --in1 $sample_path/${sample}_R1.fastq.gz \
    --in2 $sample_path/${sample}_R2.fastq.gz \
    --out1 $output_path/${sample}_fastp_dedup_adap_R1.fastq.gz \
    --out2 $output_path/${sample}_fastp_dedup_adap_R2.fastq.gz --adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=CTGTCTCTTATACACATCT -g --dedup \
    -h $output_path/${sample}_fastp_dedup_adap_trimming.html -R ${sample}_fastp_dedup_adap_trimming.out

echo "Trimming $sample with removing duplicates is done"


