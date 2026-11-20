################### Peak calling (MACS3), IDR, pooled replicates ###################

import math


def group_peaks_for_idr(wildcards):
    return expand(
        "results/peaks/macs3/{sample}_{build}_peaks.idr_input.narrowPeak",
        sample=GROUPS[wildcards.group],
        build=wildcards.build,
    )


def group_final_bams(wildcards):
    return expand(
        "results/alignment/{sample}/{sample}_{build}.final.bam",
        sample=GROUPS[wildcards.group],
        build=wildcards.build,
    )


rule macs3_callpeak:
    """Relaxed peaks (p < pvalue) on Tn5 insertion sites, ENCODE ATAC style."""
    input:
        bam=rules.tn5_shift.output.bam,
        bai=rules.tn5_shift.output.bai,
    output:
        peaks="results/peaks/macs3/{sample}_{build}_peaks.narrowPeak",
        summits="results/peaks/macs3/{sample}_{build}_summits.bed",
        xls="results/peaks/macs3/{sample}_{build}_peaks.xls",
    params:
        extra=config["peaks"]["macs3"],
        pvalue=config["peaks"]["pvalue"],
        gsize=config["genome"]["effective_genome_size"],
        name="{sample}_{build}",
        outdir="results/peaks/macs3",
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/macs3/{sample}_{build}.log",
    conda:
        "../envs/macs3.yaml"
    shell:
        # -f BAM (not BAMPE): each mate is one insertion site, as in the ENCODE ATAC pipeline.
        """
        macs3 callpeak -t {input.bam} -f BAM -g {params.gsize} -p {params.pvalue} \
            {params.extra} -n {params.name} --outdir {params.outdir} > {log} 2>&1
        """


rule peaks_for_idr:
    """Drop blacklisted peaks, sort by p-value (IDR requirement), keep the top N."""
    input:
        peaks=rules.macs3_callpeak.output.peaks,
        blacklist=config["genome"]["blacklist"],
    output:
        "results/peaks/macs3/{sample}_{build}_peaks.idr_input.narrowPeak",
    params:
        max_peaks=config["peaks"]["max_peaks"],
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "logs/macs3/{sample}_{build}_idr_input.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        # awk instead of head: head would close the pipe early and trip pipefail.
        """
        (bedtools intersect -nonamecheck -v -a {input.peaks} -b {input.blacklist} |
        LC_ALL=C sort -k8,8gr |
        awk 'NR <= {params.max_peaks}' > {output} &&
        echo "$(wc -l < {output}) peaks kept for IDR") > {log} 2>&1
        """


rule idr:
    input:
        peaks=group_peaks_for_idr,
    output:
        values="results/peaks/idr/{group}_{build}.idr_values.txt",
        optimal="results/peaks/idr/{group}_{build}.idr.optimal.narrowPeak",
    params:
        threshold=config["peaks"]["idr_threshold"],
        # IDR column 5 = min(int(-125 * log2(IDR)), 1000)
        min_score=int(-125 * math.log2(config["peaks"]["idr_threshold"])),
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/idr/{group}_{build}.log",
    conda:
        "../envs/idr.yaml"
    shell:
        # Uses the first two replicates of the group.
        """
        (idr --samples {input.peaks[0]} {input.peaks[1]} \
            --input-file-type narrowPeak --rank p.value \
            --output-file {output.values} --plot &&
        awk -v min={params.min_score} '$5 >= min' {output.values} > {output.optimal} &&
        echo "$(wc -l < {output.optimal}) peaks pass IDR < {params.threshold}") > {log} 2>&1
        """


rule merge_replicates:
    input:
        bams=group_final_bams,
    output:
        bam=temp("results/bigwig/{group}_{build}.merged.bam"),
        bai=temp("results/bigwig/{group}_{build}.merged.bam.bai"),
    threads: 8
    resources:
        mem_mb=8000,
        runtime=480,
    log:
        "logs/merge_replicates/{group}_{build}.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        (samtools merge -f -@ {threads} {output.bam} {input.bams} &&
        samtools index {output.bam}) > {log} 2>&1
        """


rule pooled_bigwig:
    input:
        bam=rules.merge_replicates.output.bam,
        bai=rules.merge_replicates.output.bai,
    output:
        bw="results/bigwig/{group}_{build}_pooled_" + config["bigwig"]["normalize"] + ".bw",
    params:
        args=rules.bigwig.params.args,
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/bigwig/{group}_{build}_pooled.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage --numberOfProcessors {threads} {params.args} \
            --bam {input.bam} -o {output.bw} > {log} 2>&1
        """
