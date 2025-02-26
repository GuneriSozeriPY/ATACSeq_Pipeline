rule fastp:
    input:
        read1="resources/samples/{sample}_1.fastq.gz", 
        read2="resources/samples/{sample}_2.fastq.gz",
    output:
        trimmed1="results/processed_files/{sample}_fastp_dedup_adap_1.fastq.gz", 
        trimmed2="results/processed_files/{sample}_fastp_dedup_adap_2.fastq.gz", 
        htmlfastp="results/processed_files/{sample}_fastp_dedup_adap_trimming.html", 
    threads: 16
    log:
        "logs/fastp/{sample}.log",  
    benchmark:
        "logs/fastp/{sample}.benchmark.log",  
    conda:
        "../envs/fastp.yaml",
    shell:
        """
        (echo "`date -R`: fastp starts..." &&
        fastp --thread {threads} --in1 {input.read1} \
        --in2 {input.read2} \
        --out1 {output.trimmed1} \
        --out2 {output.trimmed2} --adapter_sequence=CTGTCTCTTATACACATCT --adapter_sequence_r2=CTGTCTCTTATACACATCT -g --dedup \
        -h {output.htmlfastp} &&
        echo "`date -R`: fastp is successful!" || 
        (echo "`date -R`: Process failed..."; exit 1)) \
        > {log} 2>&1
        """

