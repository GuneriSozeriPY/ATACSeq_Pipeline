rule macs2:
    input:
        blacklistfiltered_bam=rules.filter_blacklist.output.blacklistfiltered_bam,
    output:
        bed="/results/processed_files/{sample}_{build}_narrow_peaks.narrowPeak",
        outdir="/results/processed_files/",
        name="{sample}_{build}_narrow",
    params:
        extra_macs2="--format BAM --bdg --SPMR --nomodel --shift -37 --extsize 73 -g 2701495711 --keep-dup all -q 0.05",
    threads: 32  
    log:
        "logs/rule/macs2/{sample}_{build}_macs2.log",
    benchmark:
        "logs/rule/macs2/{sample}_{build}_macs2.benchmark.txt",
    conda:
        "../envs/macs2.yaml"
    shell:  
        """
        (echo "`date -R`: Finding Peaks..." &&
        macs2 callpeak {params.extra_macs2} -n {output.name} -t {input.blacklistfiltered_bam} --outdir {output.outdir} &&
        echo "`date -R`: Success! Finding Peaks is done." || 
        (echo "`date -R`: Process failed..."; exit 1)) > {log} 2>&1
        """