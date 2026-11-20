"""
Collect per-sample ATAC-seq QC into one table (run via Snakemake `script:`).

Writes
- a plain TSV for you to open in R/Excel,
- the same table as MultiQC custom content (*_mqc.tsv),
- every sample's TSS profile as a MultiQC line graph (*_mqc.yaml).
"""

import pysam


def read_idxstats(path):
    total, chrm = 0, 0
    with open(path) as fh:
        for line in fh:
            chrom, _length, mapped, _unmapped = line.rstrip("\n").split("\t")
            if chrom == "*":
                continue
            total += int(mapped)
            if chrom == "chrM":
                chrm += int(mapped)
    return total, chrm


def read_two_line_table(path):
    """Header line + one value line -> dict."""
    with open(path) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return dict(zip(rows[0], rows[1]))


def read_profile(path, step):
    xs, ys = [], []
    with open(path) as fh:
        next(fh)
        for line in fh:
            pos, _raw, _norm, smoothed = line.rstrip("\n").split("\t")
            if int(pos) % step == 0:
                xs.append(pos)
                ys.append(smoothed)
    return xs, ys


def pct(numerator, denominator):
    return f"{100 * numerator / denominator:.2f}" if denominator else "NA"


def main(sm):
    samples = sm.params.samples
    conditions = sm.params.conditions
    columns = [
        "Sample",
        "Group",
        "Mapped_reads",
        "chrM_pct",
        "Picard_dup_pct",
        "Final_read_pairs",
        "TSS_enrichment",
        "Peaks_MACS3",
        "FRiP_MACS3_pct",
        "FRiP_IDR_pct",
    ]
    rows = []
    for i, sample in enumerate(samples):
        mapped, chrm = read_idxstats(sm.input.idxstats[i])
        dup = read_two_line_table(sm.input.dup[i]).get("PERCENT_DUPLICATION", "NA")
        with pysam.AlignmentFile(sm.input.bams[i]) as bam:
            final_pairs = bam.mapped // 2
        tss = read_two_line_table(sm.input.tss[i])["tss_enrichment"]
        frip = read_two_line_table(sm.input.frip[i])
        rows.append(
            [
                sample,
                conditions.get(sample, "NA"),
                str(mapped),
                pct(chrm, mapped),
                f"{100 * float(dup):.2f}" if dup != "NA" else "NA",
                str(final_pairs),
                tss,
                frip["peaks_macs3"],
                frip["frip_macs3_pct"],
                frip["frip_idr_pct"],
            ]
        )

    with open(sm.output.table, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    with open(sm.output.mqc_table, "w") as fh:
        fh.write(
            "# id: 'atac_qc_summary'\n"
            "# section_name: 'ATAC-seq QC summary'\n"
            "# description: 'Mapped reads and chrM% are before filtering. "
            "Picard duplicate % is on canonical chromosomes. Final read pairs "
            "are after MAPQ, duplicate, non-canonical/chrM, blacklist and orphan-mate "
            "filtering. TSS enrichment uses Tn5 insertion sites at GENCODE TSSs "
            "(+/-2 kb). Peaks = merged MACS3 peaks with q &lt; 0.05. FRiP = % of "
            "final reads in those peaks, or in the condition&#39;s IDR peaks.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            "#     id: 'atac_qc_summary_table'\n"
            "#     namespace: 'ATAC QC'\n"
        )
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")

    # TSS profiles as a MultiQC line graph (YAML custom content: numeric x keys).
    step = 10
    with open(sm.output.mqc_profile, "w") as fh:
        fh.write(
            "id: 'tss_profile'\n"
            "section_name: 'TSS enrichment profile'\n"
            "description: 'Aggregate Tn5 insertions around GENCODE TSSs, normalized "
            "to the outer 100 bp of the +/-2 kb window (11 bp rolling mean, every "
            "10th bp shown).'\n"
            "plot_type: 'linegraph'\n"
            "pconfig:\n"
            "  id: 'tss_profile_plot'\n"
            "  title: 'TSS enrichment profile'\n"
            "  xlab: 'Distance from TSS (bp)'\n"
            "  ylab: 'Normalized insertions'\n"
            "data:\n"
        )
        for i, sample in enumerate(samples):
            xs, ys = read_profile(sm.input.profiles[i], step)
            fh.write(f"  '{sample}':\n")
            for x, y in zip(xs, ys):
                fh.write(f"    {x}: {y}\n")


main(snakemake)  # noqa: F821
