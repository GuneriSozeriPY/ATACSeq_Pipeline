################### Differential accessibility with DiffBind (per cell line) ###################


def cell_line_idr_peaks(wildcards):
    return [
        f"results/peaks/idr/{g}_{wildcards.build}.idr.optimal.narrowPeak"
        for g in groups_of_cell_line(wildcards.cell_line)
        if len(GROUPS[g]) >= 2
    ]


def cell_line_bams(wildcards):
    return [
        f"results/alignment/{s}/{s}_{wildcards.build}.final.bam"
        for s in samples_of_cell_line(wildcards.cell_line)
    ]


rule consensus_peaks:
    """Union of all IDR peak sets of one cell line, overlapping peaks merged."""
    input:
        cell_line_idr_peaks,
    output:
        bed="results/differential/{cell_line}_{build}/consensus_peaks.bed",
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "logs/differential/{cell_line}_{build}_consensus.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        (cat {input} | cut -f1-3 | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > {output.bed} &&
        echo "$(wc -l < {output.bed}) consensus peaks") > {log} 2>&1
        """


def cell_line_peaks(wildcards):
    return [
        f"results/peaks/macs3/{s}_{wildcards.build}_peaks.idr_input.narrowPeak"
        for s in samples_of_cell_line(wildcards.cell_line)
    ]


rule diffbind:
    """DiffBind (DESeq2 inside): counts fragments in consensus peaks, ~ Replicate + Condition."""
    input:
        consensus=rules.consensus_peaks.output.bed,
        bams=cell_line_bams,
        bais=lambda w: [b + ".bai" for b in cell_line_bams(w)],
        # Per-sample peaks only register samples in DiffBind; counting uses the consensus.
        peaks=cell_line_peaks,
    output:
        summary="results/differential/{cell_line}_{build}/summary.tsv",
        results="results/differential/{cell_line}_{build}/all_contrasts.tsv",
        normalized="results/differential/{cell_line}_{build}/normalized_counts.tsv",
        plots="results/differential/{cell_line}_{build}/plots.pdf",
        rds="results/differential/{cell_line}_{build}/diffbind.rds",
    params:
        cell_line="{cell_line}",
        samples=lambda w: samples_of_cell_line(w.cell_line),
        conditions=lambda w: [SAMPLES[s]["condition"] for s in samples_of_cell_line(w.cell_line)],
        replicates=lambda w: [SAMPLES[s]["replicate"] for s in samples_of_cell_line(w.cell_line)],
        reference=config["differential"]["reference"],
        fdr=config["differential"]["fdr"],
        normalization=config["differential"]["normalization"],
        summit_width=config["differential"]["summit_width"],
        outdir="results/differential/{cell_line}_{build}",
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    log:
        "logs/differential/{cell_line}_{build}_diffbind.log",
    conda:
        "../envs/diffbind.yaml"
    script:
        "../scripts/diffbind.R"
