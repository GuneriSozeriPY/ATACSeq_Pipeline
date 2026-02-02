#!/bin/env python

import warnings
import os
import subprocess

################### Helper Functions ###########################################

def getPaired(sample, read, sample_dir):

    """
    Finds the suffix of a given sample name.

    Used rules: rename_raw, rename_raw_input 
    """

    pairedR1 = f"{sample_dir}{sample}_R1.fastq.gz"
    paired1 = f"{sample_dir}{sample}_1.fastq.gz"
    
    if os.path.isfile(pairedR1) and read == "forward":
        return f"{sample_dir}{sample}_R1.fastq.gz"

    elif os.path.isfile(pairedR1) and read == "reverse":
        return f"{sample_dir}{sample}_R2.fastq.gz"

    elif os.path.isfile(paired1) and read == "forward":
        return f"{sample_dir}{sample}_1.fastq.gz"

    elif os.path.isfile(paired1) and read == "reverse":
        return f"{sample_dir}{sample}_2.fastq.gz"
    else:
        return ""


def allInput(build, metadata, conditions=None):
    """
    Generates target files for Rule All.
    Args:
        build: Genome build string
        metadata: Dictionary of individual samples
        conditions: Dictionary of IDR conditions (optional)
    """
    inputlist = []

    # --- 1. Per-Sample Outputs ---
    for sample in metadata.keys():
        sdir = "results/processed_files"
        # QC & Reads
        inputlist.append(f"{sdir}/{sample}_1_fastqc.html")
        inputlist.append(f"{sdir}/{sample}_2_fastqc.html")
        inputlist.append(f"{sdir}/{sample}_fastp_dedup_adap_1_fastqc.html")
        inputlist.append(f"{sdir}/{sample}_fastp_dedup_adap_2_fastqc.html")
        
        # Alignment
        inputlist.append(f"results/{sample}/{sample}_{build}_sorted.dup.txt")
        inputlist.append(f"results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam")

        # Individual BigWig & Metrics
        inputlist.append(f"results/{sample}/{sample}_{build}_coverage_RPKM.bw")
        inputlist.append(f"results/macs2/{sample}_{build}_narrow_summits.bed")
        inputlist.append(f"results/macs2/{sample}_{build}_narrow_peaks.narrowPeak")
        inputlist.append(f"results/macs2/{sample}_{build}_narrow_peaks.sorted.narrowPeak")
        inputlist.append(f"{sdir}/{sample}_{build}_blacklist-filtered.fraglen.pdf")
        # inputlist.append(f"results/macs3_hmmratac/{sample}_{build}_hmmratac_cutoff_analysis.tsv")
        # inputlist.append(f"results/macs3_hmmratac/{sample}_{build}_hmmratac_final_accessible_regions.gappedPeak")
    # --- 2. IDR & Pooled Outputs (New) ---
    if conditions:
        for condition in conditions.keys():
            # The final high-confidence peak set
            inputlist.append(f"results/idr/{condition}_{build}.idr.optimal.bed")
            inputlist.append(f"results/idr/{condition}_{build}.idr_values.txt")
            # The pooled BigWig file (merged replicates)
            inputlist.append(f"results/merged/{condition}_{build}_pooled_coverage.bw")
    # Genome Index
    inputlist.append(f"resources/ref_genomes/{build}/genome_{build}.fa.fai")

    return inputlist

