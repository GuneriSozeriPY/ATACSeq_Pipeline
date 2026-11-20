# ATAC-seq pipeline (Snakemake)

Paired-end ATAC-seq from FASTQ to peaks, QC and differential accessibility.

## Setup

```
conda create -c conda-forge -c bioconda -n atac_snakemake snakemake=8.11.3
conda activate atac_snakemake
pip install snakemake-executor-plugin-slurm
```

## Inputs

- `config/samples.tsv` — one row per sample: `sample`, `cell_line`, `condition`, `replicate`,
  `fastq_1`, `fastq_2`. Replicates of one `cell_line` + `condition` form a group
  (IDR, pooled bigwig). Differential analysis runs per `cell_line`.
- `config/config.yaml` — genome, blacklist, trimming, alignment, peak-calling and DESeq2 settings.

## Run

```
cd ATACSeq_Pipeline
# From the VS Code terminal (itself a SLURM job) the SLURM variables must be unset:
env -u SLURM_JOB_ID -u SLURM_JOBID snakemake --profile config/slurm/
```

Add `-n` for a dry run. The profile (`config/slurm/config.yaml`) submits to the `genomics`
partition with the `adelab` account; per-rule memory and runtime are set in the rules.

## Steps

1. FastQC on raw reads; fastp adapter + poly-G trimming (no deduplication); FastQC on trimmed reads
2. bowtie2 (`--local --very-sensitive -X 2000`), sorted during alignment
3. Keep chr1-22, X, Y; Picard MarkDuplicates; keep proper pairs, MAPQ >= 30, no duplicates
4. Remove blacklisted reads, then remove orphan mates (`samtools fixmate`) -> `*.final.bam`
5. Tn5 shift (`alignmentSieve --ATACshift`) -> MACS3 callpeak (`-f BAM --shift -75 --extsize 150 -p 0.01 --call-summits`)
6. Blacklisted peaks removed, top 300k peaks per replicate -> IDR (threshold 0.05)
7. Bigwigs per sample and pooled per group (RPKM)
8. QC: TSS enrichment, FRiP, fragment lengths, summary table, MultiQC
9. Differential: DiffBind per cell line on the consensus of IDR peaks (summit-centred, fragments counted),
   library-size normalization, DESeq2 `~ Replicate + Condition`, each condition vs `noUV`

## Main outputs

| Path | Content |
|---|---|
| `results/alignment/{sample}/{sample}_hg38.final.bam` | filtered alignments |
| `results/bigwig/` | per-sample and pooled coverage |
| `results/peaks/macs3/`, `results/peaks/idr/` | relaxed per-sample peaks; IDR peak sets per group |
| `results/qc/hg38_atac_qc_summary.tsv` | per-sample QC table |
| `results/multiqc/hg38_multiqc_report.html` | MultiQC report |
| `results/differential/{cell_line}_hg38/` | consensus peaks, DiffBind tables per contrast, summary, plots (PCA, heatmap, MA), DBA object (`diffbind.rds`) |
| `results/differential/{cell_line}_hg38/deseq2/` | plain DESeq2 on featureCounts (full-width consensus peaks): counts, tables per contrast, summary, plots (PCA, MA), `dds.rds` |
