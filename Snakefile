configfile: "config/config.yaml"


include: "workflows/rules/common.smk"


onstart:
    shell("mkdir -p logs/slurm")  # sbatch writes job logs here


rule all:
    input:
        allInput(),


include: "workflows/rules/prepare_genome.smk"
include: "workflows/rules/fastqc.smk"
include: "workflows/rules/fastp.smk"
include: "workflows/rules/alignment.smk"
include: "workflows/rules/peaks.smk"
include: "workflows/rules/qc.smk"
include: "workflows/rules/differential.smk"
