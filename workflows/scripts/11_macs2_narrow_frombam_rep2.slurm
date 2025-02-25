#!/bin/bash

#SBATCH -J macs2

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor
#SBATCH --cpus-per-task=10

#SBATCH --output=/cta/users/ygsuzeri/ATAC/macs2/Slurm/macs2_narrow_frombam_rep2-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/macs2/Slurm/macs2_narrow_frombam_rep2-%j.err

#SBATCH --array=0-6

#conda activate atac

module load py-macs2-2.2.7.1-gcc-9.2.0-ed2mm3j

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"
macs2_path="/cta/users/ygsuzeri/ATAC/macs2"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

# echo "Processing sample $sample, bamtobed"

echo "Processing sample $sample, finding atac-seq peaks"
#Call peaks
macs2 callpeak --format BAM --bdg --SPMR --nomodel --shift -37 --extsize 73 -g 2701495711 --keep-dup all -q 0.05 \
    -n ${sample}_narrow_frombam -t $output_path/${sample}.blacklist-filtered.bam \
    --outdir $macs2_path/ 2> $macs2_path/${sample}_macs2_narrow_frombam.log

echo "done"