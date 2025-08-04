rule bowtie2:
    input:
        fastp_sample=[rules.fastp.output.trimmed1, rules.fastp.output.trimmed2],
        bowtie2="resources/ref_genomes/{build}/Bowtie2/genome_{build}.1.bt2",
    output:
        bam=temp("results/{sample}/{sample}_{build}.bam"),
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 --seed 42",
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}.benchmark.txt",
    conda:
        "../envs/alignment.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        bowtie2 \
        -p {threads} {params.extra} -x {params.ref_genome} \
        -1 {input.fastp_sample[0]} -2 {input.fastp_sample[1]} | samtools view -bS - > {output.bam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule post_alignment:
    input:
        bowtie2=rules.bowtie2.output.bam,
    output:
        sorted_bam=temp("results/{sample}/{sample}_{build}_sorted.bam"),
    params:
        ref_genome="resources/ref_genomes/{build}/Bowtie2/genome_{build}",
        extra="--local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 --seed 42",
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}_post_alignment.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}_post_alignment.benchmark.txt",
    conda:
        "../envs/alignment.yaml"
    shell:  
        """
        (echo "`date -R`: Aligning fastq files..." &&
        samtools sort {input.bowtie2} -o {output.sorted_bam} &&
        samtools index {output.sorted_bam} &&
        echo "`date -R`: Success! Alignment is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule filter_bam:
    input:
        sorted_bam=rules.post_alignment.output.sorted_bam,
    output:
        idxstats="results/{sample}/{sample}_{build}_sorted.idxstats",
        flagstat="results/{sample}/{sample}_{build}_sorted.flagstat",
        rmChrM_bam=temp("results/{sample}/{sample}_{build}_sorted.rmChrM.bam"),
        rmChrM_reheaded_bam=temp("results/{sample}/{sample}_{build}_sorted.rmChrM.reheaded.bam"),
        marked_bam=temp("results/{sample}/{sample}_{build}_sorted.marked.bam"),
        marked_bai=temp("results/{sample}/{sample}_{build}_sorted.marked.bai"),
        dup_metrics=temp("results/{sample}/{sample}_{build}_sorted.dup.metrics"),
        dup_txt="results/{sample}/{sample}_{build}_sorted.dup.txt",
        filtered_bam=temp("results/{sample}/{sample}_{build}_sorted.filtered.bam"),
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}_filter_bam.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}_filter_bam.benchmark.txt",
    conda:
        "../envs/alignment.yaml"
    shell:  
        """
        (echo "`date -R`: filter_bam runs..." &&
        samtools idxstats {input.sorted_bam} > {output.idxstats} &&
        samtools flagstat {input.sorted_bam} > {output.flagstat} &&
        samtools view -h {input.sorted_bam} | grep -v chrM | samtools sort -O bam -o {output.rmChrM_bam} -T . &&
        echo "reheading the bam file" &&
        samtools addreplacerg -@ 16 -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o {output.rmChrM_reheaded_bam} {output.rmChrM_bam} &&
        echo "marking duplicates" &&
        picard MarkDuplicates QUIET=true INPUT={output.rmChrM_reheaded_bam} OUTPUT={output.marked_bam} \
        METRICS_FILE={output.dup_metrics} REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=. &&
        head -n 8 {output.dup_metrics} | cut -f 7,9 | grep -v ^# | tail -n 2 > {output.dup_txt} &&
        samtools view -h -b -f 2 -F 1548 -q 30 {output.marked_bam} | samtools sort -o {output.filtered_bam} &&
        samtools index {output.filtered_bam} &&
        echo "`date -R`: Success! Filtering bam is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule filter_blacklist:
    input:
        filtered_bam=rules.filter_bam.output.filtered_bam,
        blacklist="resources/blacklists/hg38-blacklist.v2.bed"
    output:
        tmp_bam=temp("results/{sample}/{sample}_{build}_sorted.tmp.bam"),
        blacklistfiltered_bam="results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam",
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}_filter_blacklist.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}_filter_blacklist.benchmark.txt",
    conda:
        "../envs/alignment.yaml"
    shell:  
        """
        (echo "`date -R`: Filtering blscklist regions runs..." &&
        bedtools intersect -nonamecheck -v -abam {input.filtered_bam} -b {input.blacklist} > {output.tmp_bam} &&
        samtools sort -O bam -o {output.blacklistfiltered_bam} {output.tmp_bam} &&
        samtools index {output.blacklistfiltered_bam} &&
        echo "`date -R`: Success! Filtering blscklist regions is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule bam2bw:
    input:
        blacklistfiltered_bam=rules.filter_blacklist.output.blacklistfiltered_bam,
    output:
        bw="results/{sample}/{sample}_{build}_coverage_RPKM.bw",
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}_bam2bw.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}_bam2bw.benchmark.txt",
    params:
        bam2bw="--binSize 10 --normalizeUsing RPKM --effectiveGenomeSize 2701495711",
    conda:
        "../envs/deeptools.yaml"
    shell:  
        """
        (echo "`date -R`: Generating bw file..." &&
        bamCoverage --numberOfProcessors {threads} {params.bam2bw} \
        --bam {input.blacklistfiltered_bam} -o {output.bw} &&
        echo "`date -R`: Success! Generating bw file done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """

rule legthdist:
    input:
        blacklistfiltered_bam=rules.filter_blacklist.output.blacklistfiltered_bam,
    output:
        fraglen_stats="results/{sample}/{sample}_{build}_blacklist-filtered.fraglen.stats",
        fraglen_pdf="results/processed_files/{sample}_{build}_blacklist-filtered.fraglen.pdf",
    threads: 16  
    log:
        "logs/rule/bowtie2/{sample}_{build}_legthdist.log",
    benchmark:
        "logs/rule/bowtie2/{sample}_{build}_legthdist.benchmark.txt",
    conda:
        "../envs/picard.yaml"
    shell:  
        """
        (echo "`date -R`: Generating length distribution of counts..." &&
        picard CollectInsertSizeMetrics \
        I={input.blacklistfiltered_bam} \
        O={output.fraglen_stats} \
        H={output.fraglen_pdf} M=0.5 &&
        echo "`date -R`: Success! Generating length distribution of counts is done." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  )  > {log} 2>&1
        """