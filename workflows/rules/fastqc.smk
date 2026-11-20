################### Raw reads: clean-named links and FastQC ###################

# Links resources/fastq/{sample}_R{1,2}.fastq.gz -> the FASTQ named in the sample sheet,
# so every downstream file (and FastQC/MultiQC) uses the clean sample name.
localrules: link_fastq


rule link_fastq:
    input:
        lambda w: SAMPLES[w.sample][f"fastq_{w.read}"],
    output:
        "resources/fastq/{sample}_R{read}.fastq.gz",
    shell:
        "ln -sfn $(realpath {input}) {output}"


rule fastqc:
    input:
        rules.link_fastq.output,
    output:
        html="results/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/fastqc/{sample}_R{read}_fastqc.zip",
    threads: 2
    resources:
        mem_mb=4000,
        runtime=240,
    log:
        "logs/fastqc/{sample}_R{read}.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC on raw reads..." &&
        fastqc -t {threads} -o results/fastqc {input} &&
        echo "`date -R`: FastQC done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
