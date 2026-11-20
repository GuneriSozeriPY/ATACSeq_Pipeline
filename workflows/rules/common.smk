################### Sample sheet and helper functions ###################

import csv


def load_samples(path):
    """Read the sample sheet into {sample: row}, keeping file order."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    required = {"sample", "cell_line", "condition", "replicate", "fastq_1", "fastq_2"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Sample sheet {path} is missing columns: {sorted(missing)}")
    samples = {}
    for row in rows:
        if row["sample"] in samples:
            raise ValueError(f"Duplicate sample name in {path}: {row['sample']}")
        samples[row["sample"]] = row
    return samples


SAMPLES = load_samples(config["samples"])
BUILD = config["genome"]["build"]


def group_of(sample):
    """Replicate group (used for IDR and pooling), e.g. 'Hela_1h'."""
    row = SAMPLES[sample]
    return f"{row['cell_line']}_{row['condition']}"


GROUPS = {}
for _s in SAMPLES:
    GROUPS.setdefault(group_of(_s), []).append(_s)

CELL_LINES = sorted({row["cell_line"] for row in SAMPLES.values()})


def samples_of_cell_line(cell_line):
    return [s for s, row in SAMPLES.items() if row["cell_line"] == cell_line]


def groups_of_cell_line(cell_line):
    return [g for g, members in GROUPS.items() if SAMPLES[members[0]]["cell_line"] == cell_line]


def contrasts_of_cell_line(cell_line):
    """Conditions compared against the reference condition in one cell line."""
    ref = config["differential"]["reference"]
    conditions = []
    for s in samples_of_cell_line(cell_line):
        c = SAMPLES[s]["condition"]
        if c != ref and c not in conditions:
            conditions.append(c)
    return conditions


wildcard_constraints:
    sample="|".join(SAMPLES),
    group="|".join(GROUPS),
    cell_line="|".join(CELL_LINES),
    build=BUILD,
    read="1|2",


################### Targets for rule all ###################


def allInput():
    build = BUILD
    targets = []

    for sample in SAMPLES:
        for read in (1, 2):
            targets.append(f"results/fastqc/{sample}_R{read}_fastqc.html")
            targets.append(f"results/fastqc/{sample}_R{read}.trimmed_fastqc.html")
        targets.append(f"results/alignment/{sample}/{sample}_{build}.final.bam")
        targets.append(f"results/alignment/{sample}/{sample}_{build}.dup.txt")
        targets.append(f"results/alignment/{sample}/{sample}_{build}.fraglen.pdf")
        targets.append(f"results/bigwig/{sample}_{build}_{config['bigwig']['normalize']}.bw")
        targets.append(f"results/peaks/macs3/{sample}_{build}_peaks.narrowPeak")
        targets.append(f"results/qc/{sample}_{build}.tss_enrichment.txt")
        targets.append(f"results/qc/{sample}_{build}.frip.txt")

    for group, members in GROUPS.items():
        if len(members) >= 2:
            targets.append(f"results/peaks/idr/{group}_{build}.idr.optimal.narrowPeak")
        targets.append(f"results/bigwig/{group}_{build}_pooled_{config['bigwig']['normalize']}.bw")

    for cell_line in CELL_LINES:
        if contrasts_of_cell_line(cell_line):
            targets.append(f"results/differential/{cell_line}_{build}/summary.tsv")

    targets.append(f"results/qc/{build}_atac_qc_summary.tsv")
    targets.append(f"results/multiqc/{build}_multiqc_report.html")
    targets.append(f"resources/ref_genomes/{build}/genome_{build}.fa.fai")
    return targets
