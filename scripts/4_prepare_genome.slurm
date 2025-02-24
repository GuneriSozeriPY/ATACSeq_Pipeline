#!/bin/bash

#SBATCH -J prepare_genome

#SBATCH --account=investor

#SBATCH --qos=longer_investor

#SBATCH --partition=longer_investor

#SBATCH --output=/cta/users/ygsuzeri/ATAC/alignment/Slurm/prepare_genome-%j.out
#SBATCH --error=/cta/users/ygsuzeri/ATAC/alignment/Slurm/prepare_genome-%j.err

#SBATCH --cpus-per-task=32
#conda activate atac

module load samtools-1.9-gcc-9.2.0-w7pulwi

cd /cta/users/ygsuzeri/ATAC/alignment/reference

#to download the more recent GRCh38
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.p14.genome.fa.gz -P /cta/users/ygsuzeri/ATAC/alignment/reference/.
gunzip GRCh38.p14.genome.fa.gz

#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx="/cta/users/ygsuzeri/ATAC/alignment/reference/"

bowtie2-build --threads 32 $bt2idx/GRCh38.p14.genome.fa $bt2idx/GRCh38.p14.genome.build

#Generate an index file
samtools faidx $bt2idx/GRCh38.p14.genome.fa

