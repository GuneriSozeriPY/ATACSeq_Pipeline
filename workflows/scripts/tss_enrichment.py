"""
TSS enrichment score for one ATAC-seq BAM (run via Snakemake `script:`).

Method
- Every aligned read gives one Tn5 insertion site at its 5' end, shifted
  +4 bp on the forward strand and -5 bp on the reverse strand.
- Insertions within +/-2000 bp of every TSS are summed into one aggregate,
  strand-aware profile (minus-strand TSSs are flipped).
- The profile is divided by the mean of its outer 100 bp on each side
  (background), then smoothed with an 11 bp rolling mean.
- TSS enrichment = maximum of the smoothed profile within +/-100 bp of the TSS.

Insertion-based scores run somewhat higher than ENCODE's read-coverage
score, so compare samples with each other rather than against ENCODE
thresholds exactly.
"""

from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pysam  # noqa: E402

FLANK = 2000  # window half-width (bp)
EDGE = 100  # background: outer EDGE bp at each end of the window
SMOOTH = 11  # rolling-mean width (bp)
CENTER = 100  # score = max within +/-CENTER bp of the TSS


def load_tss(path):
    """Return {chrom: (positions array, is_plus array)} from a 6-column BED."""
    starts = defaultdict(list)
    plus = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            starts[fields[0]].append(int(fields[1]))
            plus[fields[0]].append(fields[5] != "-")
    return {c: (np.asarray(starts[c]), np.asarray(plus[c])) for c in starts}


def insertion_sites(bam, chrom):
    """Sorted Tn5 insertion positions (0-based) of all reads on one chromosome."""
    cuts = []
    append = cuts.append
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.is_reverse:
            append(read.reference_end - 1 - 5)
        else:
            append(read.reference_start + 4)
    return np.sort(np.asarray(cuts, dtype=np.int64))


def aggregate_profile(bam_path, tss):
    width = 2 * FLANK + 1
    profile = np.zeros(width, dtype=np.int64)
    n_tss = 0
    with pysam.AlignmentFile(bam_path) as bam:
        present = set(bam.references)
        for chrom, (pos, is_plus) in tss.items():
            if chrom not in present:
                continue
            cuts = insertion_sites(bam, chrom)
            n_tss += len(pos)
            if cuts.size == 0:
                continue
            lo = np.searchsorted(cuts, pos - FLANK, side="left")
            hi = np.searchsorted(cuts, pos + FLANK, side="right")
            offsets = []
            for t, p, a, b in zip(pos, is_plus, lo, hi):
                if b > a:
                    off = cuts[a:b] - t
                    offsets.append(off if p else -off)
            if offsets:
                profile += np.bincount(
                    np.concatenate(offsets) + FLANK, minlength=width
                )
    return profile, n_tss


def main(sm):
    sample = sm.params.sample
    tss = load_tss(sm.input.tss)
    profile, n_tss = aggregate_profile(sm.input.bam, tss)

    background = np.concatenate([profile[:EDGE], profile[-EDGE:]]).mean()
    normalized = profile / background if background > 0 else np.zeros(profile.shape)
    # Pad with edge values so the rolling mean doesn't dip at the window ends.
    padded = np.pad(normalized, SMOOTH // 2, mode="edge")
    smoothed = np.convolve(padded, np.ones(SMOOTH) / SMOOTH, mode="valid")
    score = float(smoothed[FLANK - CENTER : FLANK + CENTER + 1].max())

    with open(sm.output.score, "w") as fh:
        fh.write("sample\ttss_enrichment\tn_tss\tinsertions_in_windows\n")
        fh.write(f"{sample}\t{score:.3f}\t{n_tss}\t{int(profile.sum())}\n")

    positions = np.arange(-FLANK, FLANK + 1)
    with open(sm.output.profile, "w") as fh:
        fh.write("position\tinsertions\tnormalized\tsmoothed\n")
        for x, raw, norm, sm_val in zip(positions, profile, normalized, smoothed):
            fh.write(f"{x}\t{raw}\t{norm:.4f}\t{sm_val:.4f}\n")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(positions, smoothed, color="#1f77b4", lw=1.2)
    ax.axvline(0, color="grey", lw=0.6, ls="--")
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Normalized insertions")
    ax.set_title(f"{sample}\nTSS enrichment = {score:.2f}", fontsize=9)
    fig.tight_layout()
    fig.savefig(sm.output.png, dpi=150)


main(snakemake)  # noqa: F821
