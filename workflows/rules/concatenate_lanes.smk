rule concatenate_lanes:
    input:
        i1="resources/{sample}_L001_R1_001.fastq.gz",        
        i2="resources/{sample}_L002_R1_001.fastq.gz",        
    output:
        "results/{sample}_R1_001.fastq",
    threads: 16
    log:
        "logs/concatenate_lanes/{sample}.log",
    shell:
        """
        ( echo "`date -R`: concatenate_lanes file..." && 
        zcat {input.i1} {input.i2} > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """