rule macs2:
    input:
        blacklistfiltered_bam=rules.filter_blacklist.output.blacklistfiltered_bam,  
    output:
        summits="results/processed_files/{sample}_{build}_narrow_summits.bed",
    params:
        extra="--format BAM --bdg --SPMR --nomodel --shift -100 --extsize 200 -g 2701495711 --keep-dup all -q 0.05",
        name="{sample}_{build}_narrow",
    threads: 5
    log:
        "logs/rule/macs2/{sample}_{build}.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak {params.extra} \
            -t {input.blacklistfiltered_bam} \
            -n {params.name} \
            --outdir results/processed_files/ \
            > {log} 2>&1 
        """
