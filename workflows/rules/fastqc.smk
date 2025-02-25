rule fastqc:
    input:
        read1="resources/samples/{sample}_1.fastq.gz", 
        read2="resources/samples/{sample}_2.fastq.gz",
    output:
        html="results/processed_files/{sample}_1_fastqc.html", 
        zip="results/processed_files/{sample}_1_fastqc.zip", 
        html2="results/processed_files/{sample}_2_fastqc.html", 
        zip2="results/processed_files/{sample}_2_fastqc.zip", 
    wildcard_constraints:
        sample='|'.join(config["meta"].keys()),
    log:
        "logs/rule/fastqc/{sample}.log",
    benchmark:
        "logs/rule/fastqc/{sample}.benchmark.txt",
    log:
        "logs/fastqc/{sample}.log",
    benchmark:
        "logs/fastqc/{sample}.benchmark.log",
    conda:
        "../envs/fastqc.yaml",
    shell:
        """
       (echo "`date -R`: fastqc starts..." &&
        fastqc -t {threads} -o results/processed_files {input.read1} &&
        fastqc -t {threads} -o results/processed_files {input.read2} &&
        echo "`date -R`: fastqc is successful!" || 
        (echo "`date -R`: Process failed..."; exit 1)) \
        >> {log} 2>&1
        """

