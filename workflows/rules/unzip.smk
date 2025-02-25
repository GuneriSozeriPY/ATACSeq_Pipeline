rule unzip:
    input:
        "resources/{sample}.fastq.gz",
    output:
        "results/{sample}.fastq",
    threads: 16
    log:
        "logs/rule1/{sample}.log",
    shell:
        """
        ( echo "`date -R`: Unzip file..." && 
        zcat {input} > {output} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """