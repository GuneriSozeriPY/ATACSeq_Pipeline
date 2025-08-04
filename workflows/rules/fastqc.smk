rule fastqc:
    input:
        read1="resources/samples/{sample}_{rep}.fastq.gz", 
    output:
        html="results/processed_files/{sample}_{rep}_fastqc.html", 
        zip="results/processed_files/{sample}_{rep}_fastqc.zip",  
    wildcard_constraints:
        sample='|'.join(config["meta"].keys()),
    threads: 16
    log:
        "logs/rule/fastqc/{sample}_{rep}.log",
    benchmark:
        "logs/rule/fastqc/{sample}_{rep}.benchmark.log",
    conda:
        "../envs/fastqc.yaml",
    shell:
        """
       (echo "`date -R`: fastqc starts..." &&
        fastqc -t {threads} -o results/processed_files {input.read1} &&
        echo "`date -R`: fastqc is successful!" || 
        (echo "`date -R`: Process failed..."; exit 1)) \
        > {log} 2>&1
        """
