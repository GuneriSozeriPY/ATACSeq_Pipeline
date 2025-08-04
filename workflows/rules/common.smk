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



def allInput(build, metadata):
    inputlist = []
    for sample in metadata.keys():
        sdir = "results/processed_files"
        sprefix = f"{sample}_{build}"
        inputlist.append(f"{sdir}/{sample}_1_fastqc.html")
        inputlist.append(f"{sdir}/{sample}_2_fastqc.html")
        inputlist.append(f"resources/ref_genomes/{build}/genome_{build}.fa.fai")
        inputlist.append(f"{sdir}/{sample}_fastp_dedup_adap_1_fastqc.html")
        inputlist.append(f"{sdir}/{sample}_fastp_dedup_adap_2_fastqc.html")
        inputlist.append(f"results/{sample}/{sample}_{build}_sorted.dup.txt")
        inputlist.append(f"results/{sample}/{sample}_{build}_sorted.blacklist-filtered.bam")
        inputlist.append(f"results/{sample}/{sample}_{build}_coverage_RPKM.bw")
        inputlist.append(f"{sdir}/{sample}_{build}_blacklist-filtered.fraglen.pdf")
        inputlist.append(f"{sdir}/{sample}_{build}_narrow_summits.bed")

    return inputlist
