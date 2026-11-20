################### ATAC-seq QC: TSS enrichment, FRiP, summary, MultiQC ###################

import math


def frip_inputs(wildcards):
    inputs = {
        "bam": f"results/alignment/{wildcards.sample}/{wildcards.sample}_{wildcards.build}.final.bam",
        "bai": f"results/alignment/{wildcards.sample}/{wildcards.sample}_{wildcards.build}.final.bam.bai",
        "macs3": f"results/peaks/macs3/{wildcards.sample}_{wildcards.build}_peaks.idr_input.narrowPeak",
    }
    group = group_of(wildcards.sample)
    if len(GROUPS[group]) >= 2:
        inputs["idr"] = f"results/peaks/idr/{group}_{wildcards.build}.idr.optimal.narrowPeak"
    return inputs


def per_sample(pattern):
    return expand(pattern, sample=list(SAMPLES), allow_missing=True)


# Runs on the machine running Snakemake: compute nodes may have no internet.
localrules: gtf_download


rule gtf_download:
    output:
        "resources/annotation/{build}/annotation_{build}.gtf.gz",
    params:
        link=config["genome"]["gtf_link"],
    log:
        "logs/gtf_download/{build}.log",
    shell:
        """
        (echo "`date -R`: Downloading annotation..." &&
        wget -q {params.link} -O {output} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule tss_bed:
    input:
        rules.gtf_download.output,
    output:
        "resources/annotation/{build}/tss_{build}.bed",
    log:
        "logs/tss_bed/{build}.log",
    shell:
        """
        # One line per unique transcript TSS (0-based BED, strand-aware).
        (zcat {input} |
        awk -F'\\t' 'BEGIN{{OFS="\\t"}} $3 == "transcript" {{
            match($9, /gene_name "[^"]+"/); name = substr($9, RSTART + 11, RLENGTH - 12);
            tss = ($7 == "+") ? $4 - 1 : $5 - 1;
            print $1, tss, tss + 1, name, ".", $7 }}' |
        LC_ALL=C sort -k1,1 -k2,2n -k6,6 -u > {output} &&
        echo "`date -R`: $(wc -l < {output}) unique TSSs written." ) > {log} 2>&1
        """


rule tss_enrichment:
    input:
        bam=rules.final_bam.output.bam,
        bai=rules.final_bam.output.bai,
        tss=rules.tss_bed.output[0],
    output:
        score="results/qc/{sample}_{build}.tss_enrichment.txt",
        profile="results/qc/{sample}_{build}.tss_profile.tsv",
        png="results/qc/{sample}_{build}.tss_profile.png",
    params:
        sample="{sample}",
    threads: 1
    log:
        "logs/qc/{sample}_{build}_tss_enrichment.log",
    benchmark:
        "logs/qc/{sample}_{build}_tss_enrichment.benchmark.txt",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/tss_enrichment.py"


rule frip:
    input:
        unpack(frip_inputs),
    output:
        frip="results/qc/{sample}_{build}.frip.txt",
    params:
        idr=lambda wildcards, input: input.idr if hasattr(input, "idr") else "",
        # MACS3 peaks are called at a relaxed p-value for IDR; FRiP uses q < frip_qvalue.
        min_q=-math.log10(config["peaks"]["frip_qvalue"]),
        tmp="results/qc/{sample}_{build}.frip_tmp",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=120,
    log:
        "logs/qc/{sample}_{build}_frip.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        (echo "`date -R`: FRiP for {wildcards.sample}..." &&
        total=$(samtools view -c -@ {threads} {input.bam}) &&
        # Peaks with q < frip_qvalue, overlapping ones merged (merged regions = peak count).
        # -M counts each read once even if it overlaps two peaks.
        awk -v q={params.min_q} '$9 >= q' {input.macs3} | cut -f1-3 |
            LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > {params.tmp}.macs3.bed &&
        in_macs3=$(samtools view -c -M -@ {threads} -L {params.tmp}.macs3.bed {input.bam}) &&
        n_macs3=$(wc -l < {params.tmp}.macs3.bed) &&
        if [ -n "{params.idr}" ]; then
            cut -f1-3 {params.idr} | LC_ALL=C sort -k1,1 -k2,2n | bedtools merge -i - > {params.tmp}.idr.bed &&
            in_idr=$(samtools view -c -M -@ {threads} -L {params.tmp}.idr.bed {input.bam})
        else
            in_idr=NA
        fi &&
        awk -v s={wildcards.sample} -v t=$total -v m3=$in_macs3 -v p3=$n_macs3 -v i=$in_idr 'BEGIN {{
            OFS = "\\t";
            print "sample", "reads", "peaks_macs3", "reads_in_macs3", "frip_macs3_pct", "reads_in_idr", "frip_idr_pct";
            fi = (i == "NA") ? "NA" : sprintf("%.2f", 100 * i / t);
            print s, t, p3, m3, sprintf("%.2f", 100 * m3 / t), i, fi }}' > {output.frip} &&
        rm -f {params.tmp}.macs3.bed {params.tmp}.idr.bed &&
        cat {output.frip} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule qc_summary:
    input:
        idxstats=per_sample(rules.filter_bam.output.idxstats),
        dup=per_sample(rules.filter_bam.output.dup_txt),
        bams=per_sample(rules.final_bam.output.bam),
        bais=per_sample(rules.final_bam.output.bai),
        tss=per_sample(rules.tss_enrichment.output.score),
        profiles=per_sample(rules.tss_enrichment.output.profile),
        frip=per_sample(rules.frip.output.frip),
    output:
        table="results/qc/{build}_atac_qc_summary.tsv",
        mqc_table="results/qc/{build}_atac_qc_summary_mqc.tsv",
        mqc_profile="results/qc/{build}_tss_profiles_mqc.yaml",
    params:
        samples=list(SAMPLES),
        conditions={s: group_of(s) for s in SAMPLES},
    resources:
        mem_mb=4000,
        runtime=60,
    log:
        "logs/qc/{build}_qc_summary.log",
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/qc_summary.py"


rule multiqc:
    input:
        summary=rules.qc_summary.output.mqc_table,
        profiles=rules.qc_summary.output.mqc_profile,
        fastqc_raw=expand(rules.fastqc.output.zip, sample=list(SAMPLES), read=[1, 2]),
        fastqc_trimmed=expand(rules.fastqc_trimmed.output.zip, sample=list(SAMPLES), read=[1, 2]),
        fastp=expand(rules.fastp.output.json, sample=list(SAMPLES)),
        flagstat=per_sample(rules.filter_bam.output.flagstat),
        idxstats=per_sample(rules.filter_bam.output.idxstats),
        dup_metrics=per_sample(rules.filter_bam.output.dup_metrics),
        fraglen=per_sample(rules.fragment_length.output.stats),
        macs3=per_sample(rules.macs3_callpeak.output.xls),
    output:
        html="results/multiqc/{build}_multiqc_report.html",
    params:
        bowtie2_logs=expand("logs/bowtie2/{sample}_{{build}}.log", sample=list(SAMPLES)),
        filelist="results/multiqc/{build}_multiqc_files.txt",
        # Runs FastQC twice (raw / trimmed) so the two don't overwrite each other.
        mqc_config="config/multiqc_config.yaml",
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        "logs/qc/{build}_multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        mkdir -p results/multiqc
        printf '%s\\n' {input} {params.bowtie2_logs} > {params.filelist}
        multiqc -f --file-list {params.filelist} -c {params.mqc_config} \\
            -o results/multiqc -n $(basename {output.html}) \\
            --cl-config "extra_fn_clean_exts: ['_{wildcards.build}']" > {log} 2>&1
        """
