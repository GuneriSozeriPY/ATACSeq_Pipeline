################### Adapter trimming (fastp) and FastQC on trimmed reads ###################

# Duplicates are NOT removed here: Picard MarkDuplicates handles them after alignment,
# which keeps library-complexity estimates meaningful.
rule fastp:
    input:
        read1="resources/fastq/{sample}_R1.fastq.gz",
        read2="resources/fastq/{sample}_R2.fastq.gz",
    output:
        trimmed1=temp("results/trimmed/{sample}_R1.trimmed.fastq.gz"),
        trimmed2=temp("results/trimmed/{sample}_R2.trimmed.fastq.gz"),
        html="results/trimmed/{sample}.fastp.html",
        json="results/trimmed/{sample}.fastp.json",
    params:
        adapter_r1=config["trimming"]["adapter_r1"],
        adapter_r2=config["trimming"]["adapter_r2"],
        extra=config["trimming"]["extra"],
    threads: 16
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/fastp/{sample}.log",
    benchmark:
        "logs/fastp/{sample}.benchmark.txt"
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        (echo "`date -R`: fastp starts..." &&
        fastp --thread {threads} \
            --in1 {input.read1} --in2 {input.read2} \
            --out1 {output.trimmed1} --out2 {output.trimmed2} \
            --adapter_sequence={params.adapter_r1} --adapter_sequence_r2={params.adapter_r2} \
            {params.extra} \
            -h {output.html} -j {output.json} &&
        echo "`date -R`: fastp is successful!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_R{read}.trimmed.fastq.gz",
    output:
        html="results/fastqc/{sample}_R{read}.trimmed_fastqc.html",
        zip="results/fastqc/{sample}_R{read}.trimmed_fastqc.zip",
    threads: 2
    resources:
        mem_mb=4000,
        runtime=240,
    log:
        "logs/fastqc/{sample}_R{read}.trimmed.log",
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        (echo "`date -R`: FastQC on trimmed reads..." &&
        fastqc -t {threads} -o results/fastqc {input} &&
        echo "`date -R`: FastQC done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
