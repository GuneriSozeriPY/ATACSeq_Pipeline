#!/bin/bash

#SBATCH -J post_aln

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor
#SBATCH --cpus-per-task=10

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/removeblacklistRegions-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/removeblacklistRegions-%j.err

#SBATCH --array=0-6

#conda activate atac

module load samtools-1.9-gcc-9.2.0-w7pulwi
module load bedtools2-2.27.1-gcc-9.2.0-joqljqy

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

# echo "Processing sample $sample, removing blacklist regions"

# #Remove reads within the blacklist regions
# bedtools intersect -nonamecheck -v -abam $output_path/${sample}.filtered.bam -b /cta/users/ygsuzeri/hg38-blacklist.v2.bed > $output_path/${sample}.tmp.bam

# #Sort and index the bam file
# samtools sort -O bam -o $output_path/${sample}.blacklist-filtered.bam $output_path/${sample}.tmp.bam
# samtools index $output_path/${sample}.blacklist-filtered.bam

# rm $output_path/${sample}.tmp.bam
# samtools index $output_path/${sample}.blacklist-filtered.bam

echo "removing unnecessary bamfiles"

rm -r $output_path/${sample}.filtered.bam.bai
rm -r $output_path/${sample}.filtered.bam

echo "done"