#!/bin/bash

#SBATCH -J atac_bowtie2

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor
#SBATCH --cpus-per-task=32

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/bowtie2-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/bowtie2-%j.err

#SBATCH --array=0-6

#conda activate atac

module load samtools-1.9-gcc-9.2.0-w7pulwi

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"
trimmed_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/adaptor_trimmimg/adapter"


filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx="/cta/users/ygsuzeri/ATAC/alignment/reference/GRCh38.p14.genome.build"

echo "Processing sample $sample , alignment starts"
#run bowtie2
bowtie2 -p 32 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $bt2idx \
    -1 $trimmed_path/${sample}_fastp_dedup_adap_R1.fastq.gz -2 $trimmed_path/${sample}_fastp_dedup_adap_R2.fastq.gz | samtools view -bS - > $output_path/${sample}.bam

echo "Alignment done, sorting and indexing starts"
#Sort the output bam file by coordinate
samtools sort $output_path/${sample}.bam -o $output_path/${sample}_sorted.bam 

echo "Sorting done, indexing starts"
#Generate an index file
samtools index $output_path/${sample}_sorted.bam

echo "Done"