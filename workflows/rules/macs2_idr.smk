# Helper to get all replicates for a specific condition
def get_replicates_bam(wildcards):
    return expand(
        "results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam",
        sample=config["conditions"][wildcards.condition],
        build=wildcards.build
    )

def get_replicates_peaks(wildcards):
    return expand(
        "results/macs2/{sample}_{build}_narrow_peaks.narrowPeak",
        sample=config["conditions"][wildcards.condition],
        build=wildcards.build
    )

def get_replicates_sorted_peaks(wildcards):
    # This now asks for the SORTED peaks, not the raw ones
    return expand(
        "results/macs2/{sample}_{build}_narrow_peaks.sorted.narrowPeak",
        sample=config["conditions"][wildcards.condition],
        build=wildcards.build
    )

# 1. Call Relaxed Peaks (Individual Replicates)
rule macs2_relaxed:
    input:
        blacklistfiltered_bam=rules.filter_blacklist.output.blacklistfiltered_bam,  
    output:
        peaks="results/macs2/{sample}_{build}_narrow_peaks.narrowPeak",
        summits="results/macs2/{sample}_{build}_narrow_summits.bed",
    params:
        extra="--format BAMPE --nomodel --shift -100 --extsize 200 -g 2701495711 --keep-dup all -q 0.05",
        name="{sample}_{build}_narrow",
    threads: 5
    log:
        "logs/macs2/{sample}_{build}_relaxed.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak {params.extra} \
            -t {input.blacklistfiltered_bam} \
            -n {params.name} \
            --outdir results/macs2/ \
            > {log} 2>&1 
        """

# 2. Sort Peaks by P-value
# HBC Requirement: "narrowPeak files have to be sorted by the -log10(p-value) column"
rule sort_peaks_for_idr:
    input:
        peaks="results/macs2/{sample}_{build}_narrow_peaks.narrowPeak"
    output:
        sorted_peaks="results/macs2/{sample}_{build}_narrow_peaks.sorted.narrowPeak"
    shell:
        """
        # Column 8 is -log10(p-value) in narrowPeak format
        # -k8,8nr means: Sort key is col 8, numeric, reverse (descending)
        sort -k8,8nr {input.peaks} > {output.sorted_peaks}
        """

# 2. Run IDR Analysis
rule run_idr:
    input:
        peaks=get_replicates_sorted_peaks
    output:
        idr_out="results/idr/{condition}_{build}.idr_values.txt",
        idr_peaks="results/idr/{condition}_{build}.idr.optimal.bed",
    log:
        "logs/idr/{condition}_{build}_relaxed.log"
    params:
        threshold=0.05, # The final IDR cutoff
        outdir="results/idr"
    conda:
        "../envs/idr.yaml"
    shell:
        """
        mkdir -p {params.outdir} &&
        # Run IDR on the first two replicates
        idr --samples {input.peaks[0]} {input.peaks[1]} \
            --input-file-type narrowPeak \
            --rank p.value \
            --output-file {output.idr_out} \
            --plot 

        # Filter for peaks passing the IDR threshold (Col 5 is scaled IDR score)
        # IDR score of 540 = -125*log2(0.05) approx
        # Standard approach: Filter rows where global IDR <= 0.05
        awk '$5 >= 540 {{print $0}}' {output.idr_out} > {output.idr_peaks}
        """

# 3. Merge Replicate BAMs
rule merge_replicates:
    input:
        bams=get_replicates_bam
    output:
        merged_bam=temp("results/merged/{condition}_{build}.merged.bam"),
        merged_idx=temp("results/merged/{condition}_{build}.merged.bam.bai")
    threads: 8
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        samtools merge -@ {threads} {output.merged_bam} {input.bams}
        samtools index {output.merged_bam}
        """

# 4. Generate Pooled BigWig
rule pooled_bam2bw:
    input:
        bam=rules.merge_replicates.output.merged_bam,
        idx=rules.merge_replicates.output.merged_idx
    output:
        bw="results/merged/{condition}_{build}_pooled_coverage.bw"
    threads: 8
    params:
        bam2bw="--binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2701495711"
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage --numberOfProcessors {threads} {params.bam2bw} \
        --bam {input.bam} -o {output.bw}
        """