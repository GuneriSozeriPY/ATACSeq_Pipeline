################### Alignment, filtering, Tn5 shift, coverage ###################


rule bowtie2:
    input:
        reads=[rules.fastp.output.trimmed1, rules.fastp.output.trimmed2],
        index="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        bam=temp("results/alignment/{sample}/{sample}_{build}.sorted.bam"),
        bai=temp("results/alignment/{sample}/{sample}_{build}.sorted.bam.bai"),
    params:
        index="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra=config["alignment"]["bowtie2"],
        tmp="results/alignment/{sample}/{sample}_{build}.sorttmp",
    threads: 16
    resources:
        mem_mb=32000,
        runtime=1440,
    log:
        "logs/bowtie2/{sample}_{build}.log",
    benchmark:
        "logs/bowtie2/{sample}_{build}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        # bowtie2's alignment summary goes to the log (MultiQC reads it from there).
        # Read group is set here, so no separate reheader step is needed.
        """
        (echo "`date -R`: Aligning and sorting..." &&
        bowtie2 -p {threads} {params.extra} -x {params.index} \
            --rg-id {wildcards.sample} --rg SM:{wildcards.sample} \
            --rg LB:{wildcards.sample} --rg PL:ILLUMINA \
            -1 {input.reads[0]} -2 {input.reads[1]} |
        samtools sort -@ 4 -m 2G -T {params.tmp} -o {output.bam} - &&
        samtools index {output.bam} &&
        echo "`date -R`: Success! Alignment is done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule filter_bam:
    """Keep canonical chromosomes, mark duplicates, keep properly paired MAPQ>=min reads."""
    input:
        bam=rules.bowtie2.output.bam,
        bai=rules.bowtie2.output.bai,
    output:
        idxstats="results/alignment/{sample}/{sample}_{build}.idxstats",
        flagstat="results/alignment/{sample}/{sample}_{build}.flagstat",
        kept=temp("results/alignment/{sample}/{sample}_{build}.chroms.bam"),
        marked=temp("results/alignment/{sample}/{sample}_{build}.marked.bam"),
        dup_metrics="results/alignment/{sample}/{sample}_{build}.dup.metrics",
        dup_txt="results/alignment/{sample}/{sample}_{build}.dup.txt",
        filtered=temp("results/alignment/{sample}/{sample}_{build}.filtered.bam"),
    params:
        keep_chroms=config["genome"]["keep_chroms"],
        min_mapq=config["alignment"]["min_mapq"],
        tmpdir="results/alignment/{sample}",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=1440,
    log:
        "logs/filter_bam/{sample}_{build}.log",
    benchmark:
        "logs/filter_bam/{sample}_{build}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        # Region queries on the indexed BAM keep coordinate order, so no re-sorting.
        # -F 1804 = unmapped, mate unmapped, secondary, QC-fail, duplicate.
        """
        (echo "`date -R`: filter_bam starts..." &&
        samtools idxstats {input.bam} > {output.idxstats} &&
        samtools flagstat {input.bam} > {output.flagstat} &&
        echo "keeping chromosomes matching {params.keep_chroms}" &&
        samtools view -@ {threads} -b -o {output.kept} {input.bam} \
            $(awk -v re='{params.keep_chroms}' '$1 ~ re && $3 > 0 {{print $1}}' {output.idxstats}) &&
        echo "marking duplicates" &&
        picard -Xmx24g MarkDuplicates --INPUT {output.kept} --OUTPUT {output.marked} \
            --METRICS_FILE {output.dup_metrics} --REMOVE_DUPLICATES false \
            --VALIDATION_STRINGENCY LENIENT --TMP_DIR {params.tmpdir} --QUIET true &&
        awk -F'\\t' '/^LIBRARY/ {{for (i = 1; i <= NF; i++) col[$i] = i; getline;
            print "READ_PAIR_DUPLICATES\\tPERCENT_DUPLICATION";
            print $col["READ_PAIR_DUPLICATES"] "\\t" $col["PERCENT_DUPLICATION"]; exit}}' \
            {output.dup_metrics} > {output.dup_txt} &&
        echo "filtering: proper pairs, MAPQ >= {params.min_mapq}, no duplicates" &&
        samtools view -@ {threads} -b -f 2 -F 1804 -q {params.min_mapq} -o {output.filtered} {output.marked} &&
        echo "`date -R`: Success! Filtering bam is done." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule final_bam:
    """Remove blacklisted reads, then drop orphan mates left by any filter."""
    input:
        bam=rules.filter_bam.output.filtered,
        blacklist=config["genome"]["blacklist"],
    output:
        bam="results/alignment/{sample}/{sample}_{build}.final.bam",
        bai="results/alignment/{sample}/{sample}_{build}.final.bam.bai",
    params:
        tmp="results/alignment/{sample}/{sample}_{build}.finaltmp",
    threads: 8
    resources:
        mem_mb=24000,
        runtime=720,
    log:
        "logs/final_bam/{sample}_{build}.log",
    benchmark:
        "logs/final_bam/{sample}_{build}.benchmark.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        # Name-sort -> fixmate flags reads whose mate was removed -> keep proper pairs only.
        """
        (echo "`date -R`: Removing blacklisted reads and orphan mates..." &&
        bedtools intersect -nonamecheck -v -abam {input.bam} -b {input.blacklist} |
        samtools sort -n -@ 4 -m 2G -T {params.tmp}.n - |
        samtools fixmate -r - - |
        samtools view -u -f 2 -F 1804 - |
        samtools sort -@ 4 -m 2G -T {params.tmp}.c -o {output.bam} - &&
        samtools index {output.bam} &&
        echo "`date -R`: Success! Final BAM written." ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule tn5_shift:
    """Shift reads to Tn5 insertion sites (+4/-5) for peak calling."""
    input:
        bam=rules.final_bam.output.bam,
        bai=rules.final_bam.output.bai,
    output:
        bam=temp("results/alignment/{sample}/{sample}_{build}.shifted.bam"),
        bai=temp("results/alignment/{sample}/{sample}_{build}.shifted.bam.bai"),
    params:
        unsorted="results/alignment/{sample}/{sample}_{build}.shifted.unsorted.bam",
        tmp="results/alignment/{sample}/{sample}_{build}.shifttmp",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=720,
    log:
        "logs/tn5_shift/{sample}_{build}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        (echo "`date -R`: Tn5 shift..." &&
        alignmentSieve --ATACshift -p {threads} -b {input.bam} -o {params.unsorted} &&
        samtools sort -@ 4 -m 2G -T {params.tmp} -o {output.bam} {params.unsorted} &&
        samtools index {output.bam} &&
        rm -f {params.unsorted} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule bigwig:
    input:
        bam=rules.final_bam.output.bam,
        bai=rules.final_bam.output.bai,
    output:
        bw="results/bigwig/{sample}_{build}_" + config["bigwig"]["normalize"] + ".bw",
    params:
        args=(
            f"--binSize {config['bigwig']['bin_size']} "
            f"--normalizeUsing {config['bigwig']['normalize']} "
            f"--effectiveGenomeSize {config['genome']['effective_genome_size']}"
        ),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/bigwig/{sample}_{build}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        (echo "`date -R`: Generating bigwig..." &&
        bamCoverage --numberOfProcessors {threads} {params.args} \
            --bam {input.bam} -o {output.bw} &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """


rule fragment_length:
    input:
        bam=rules.final_bam.output.bam,
    output:
        stats="results/alignment/{sample}/{sample}_{build}.fraglen.stats",
        pdf="results/alignment/{sample}/{sample}_{build}.fraglen.pdf",
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        "logs/fragment_length/{sample}_{build}.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        (echo "`date -R`: Fragment length distribution..." &&
        picard -Xmx12g CollectInsertSizeMetrics --INPUT {input.bam} \
            --OUTPUT {output.stats} --Histogram_FILE {output.pdf} --MINIMUM_PCT 0.5 &&
        echo "`date -R`: Success!" ||
        {{ echo "`date -R`: Process failed..."; exit 1; }} ) > {log} 2>&1
        """
