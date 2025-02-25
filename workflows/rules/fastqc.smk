rule fastqc:
    input:
        "results/{sample}.fastq",
    output:
        "results/{sample}_fastqc.html",
        "results/{sample}_fastqc.zip",
    threads: 16
    log:
        "logs/fastqc/{sample}.log",
    benchmark:
        "logs/fastqc/{sample}.benchmark.log",
    conda:
        "../envs/fastqc.yaml",
    shell:
        """
        fastqc -t {threads} -o results/ {input} > {log} 2>&1
        """ 