#!/bin/bash

#SBATCH -J bamtobw

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor
#SBATCH --cpus-per-task=10

#SBATCH --output=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/bamtobw-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/Slurm/bamtobw-%j.err

#SBATCH --array=0-6

#conda activate deeptools

sample_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/data"
output_path="/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment"

filelist=($(find $sample_path -maxdepth 1 -type f -name "*_R1.fastq.gz" | awk -F '/' '{print $NF}' | sed 's/_R1\.fastq\.gz$//'))
sample=${filelist[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample $sample, bam to bw"

bamCoverage --numberOfProcessors 10 --binSize 10 --normalizeUsing RPKM \
  --effectiveGenomeSize 2701495711 --bam $output_path/${sample}.blacklist-filtered.bam -o $output_path/${sample}_coverage_RPKM.bw

echo "done"