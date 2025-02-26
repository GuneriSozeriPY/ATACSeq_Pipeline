---
Generating an ATAC-Seq Pipeline using snakemake
---
## General Info 

This repo will be used to write an ATAC-Seq Pipeline, using the scripts under the "scripts" directory.


#### Create a conda environment with the defined packages.

```
conda create -c conda-forge -c bioconda -c r -n atac_snakemake snakemake=8.11.3
conda activate atac_snakemake
pip install snakemake-executor-plugin-cluster-generic
```

```
screen -S atac_snakemake
screen -xS atac_snakemake 
```

```
cd ATACSeq_Pipeline
snakemake --profile config/slurm/ --use-conda --conda-frontend conda --rerun-incomplete
```