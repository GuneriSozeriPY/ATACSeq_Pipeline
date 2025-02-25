#!/bin/bash

#SBATCH -J post_aln

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor
#SBATCH --cpus-per-task=10

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/post_alignment-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/post_alignment-%j.err

#SBATCH --array=0-6

#conda activate atac

module load samtools-1.9-gcc-9.2.0-w7pulwi

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

# echo "Processing sample $sample, post alignment processing starts"

# #Generate the idxstats report
# samtools idxstats $output_path/${sample}_sorted.bam > $output_path/${sample}_sorted.idxstats

# #Generate the flagstat report
# samtools flagstat $output_path/${sample}_sorted.bam > $output_path/${sample}_sorted.flagstat

# #Remove reads aligned to the mitochondria
# samtools view -h $output_path/${sample}_sorted.bam | grep -v chrM | samtools sort -O bam -o $output_path/${sample}_sorted.rmChrM.bam -T .

# echo "reheading the bam file"

# samtools addreplacerg -@ 16 -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o $output_path/${sample}_sorted.rmChrM.reheaded.bam $output_path/${sample}_sorted.rmChrM.bam

# echo "marking duplicates"
# picard MarkDuplicates QUIET=true INPUT=$output_path/${sample}_sorted.rmChrM.reheaded.bam OUTPUT=$output_path/${sample}.marked.bam \
#     METRICS_FILE=$output_path/${sample}.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

# #View the % of duplicates
# head -n 8 $output_path/${sample}.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2 > $output_path/${sample}.dup.txt

# echo "removing duplicates and filtering"
# samtools view -h -b -f 2 -F 1548 -q 30 $output_path/${sample}.marked.bam | samtools sort -o $output_path/${sample}.filtered.bam

# samtools index $output_path/${sample}.filtered.bam

# echo "removing unnecessary bamfiles"

rm -r $output_path/${sample}_sorted.rmChrM.bam
rm -r $output_path/${sample}.marked.bam
rm -r $output_path/${sample}.marked.bai
rm -r $output_path/${sample}.bam
rm -r $output_path/${sample}_sorted.rmChrM.reheaded.bam
rm -r $output_path/${sample}_sorted.bam
rm -r $output_path/${sample}_sorted.bam.bai

echo "done"