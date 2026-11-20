# Differential accessibility with DiffBind for one cell line (run via Snakemake `script:`).
#
# - Peaks: the cell line's consensus of IDR peak sets, re-centred on summits (+/- summit_width bp).
# - Counts: fragments (read pairs) from the final BAMs, counted by DiffBind.
# - Normalization: config differential.normalization
#     "lib" = total library size (DiffBind default; robust when many peaks change)
#     "RLE" = DESeq2 median-of-ratios on reads in peaks
# - Model: DESeq2, design ~ Replicate + Condition; every condition vs the reference.
#
# Outputs (in the cell line's results/differential folder)
#   <condition>_vs_<reference>.tsv   all peaks for one contrast, sorted by FDR
#   all_contrasts.tsv                the same tables stacked, with a `contrast` column
#                                    (Conc_<condition>/Conc_<reference> -> Conc_group/Conc_reference)
#   summary.tsv                      peaks with FDR < fdr: total / more / less accessible
#   normalized_counts.tsv            normalized counts per peak and sample
#   plots.pdf                        PCA, sample correlation heatmap, MA plot per contrast
#   diffbind.rds                     the analysed DBA object, for further work in R

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(DiffBind)
    library(BiocParallel)
})

p <- snakemake@params
samples <- unlist(p[["samples"]])
reference <- p[["reference"]]
fdr <- as.numeric(p[["fdr"]])
outdir <- p[["outdir"]]
threads <- as.integer(snakemake@threads)

sheet <- data.frame(
    SampleID = samples,
    Tissue = p[["cell_line"]],
    Condition = unlist(p[["conditions"]]),
    # DiffBind wants a numeric replicate; R1 / R3 -> 1 / 3 (two levels, so same model as a factor)
    Replicate = as.integer(sub("^R", "", unlist(p[["replicates"]]))),
    bamReads = unlist(snakemake@input[["bams"]]),
    Peaks = unlist(snakemake@input[["peaks"]]),
    PeakCaller = "narrow"
)
print(sheet[, c("SampleID", "Condition", "Replicate")])

dba_obj <- dba(sampleSheet = sheet)
dba_obj$config$singleEnd <- FALSE        # count fragments, not reads
dba_obj$config$cores <- threads
dba_obj$config$doBlacklist <- FALSE      # reads and peaks are already blacklist-filtered
dba_obj$config$doGreylist <- FALSE       # no input controls
register(MulticoreParam(threads))

consensus <- read.table(snakemake@input[["consensus"]], col.names = c("chr", "start", "end"))
consensus$start <- consensus$start + 1   # BED (0-based) -> GRanges (1-based)
consensus <- GenomicRanges::makeGRangesFromDataFrame(consensus)
message(length(consensus), " consensus peaks")

dba_obj <- dba.count(dba_obj, peaks = consensus, summits = p[["summit_width"]],
                     minOverlap = 1, bParallel = TRUE)
print(dba.show(dba_obj))

# RLE factors are computed with DESeq2 on the raw counts and passed to DiffBind as user factors:
# DiffBind's own DBA_NORM_RLE / DBA_LIBSIZE_PEAKREADS path calls edgeR, which this env does not have.
norm <- switch(p[["normalization"]],
    lib = DBA_NORM_LIB,
    RLE = DESeq2::estimateSizeFactorsForMatrix(sapply(dba_obj$peaks, function(pk) pk$Reads)),
    stop("differential.normalization must be 'lib' or 'RLE'")
)
dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2, normalize = norm, library = DBA_LIBSIZE_FULL)
print(dba.normalize(dba_obj, method = DBA_DESEQ2, bRetrieve = TRUE))

conditions <- unique(sheet$Condition)
contrasts <- setdiff(conditions, reference)
dba_obj <- dba.contrast(dba_obj, design = "~Replicate + Condition",
                        reorderMeta = list(Condition = reference))
for (cond in contrasts) {
    dba_obj <- dba.contrast(dba_obj, contrast = c("Condition", cond, reference))
}
dba_obj <- dba.analyze(dba_obj, method = DBA_DESEQ2, bParallel = FALSE)
print(dba.show(dba_obj, bContrasts = TRUE))

# Contrast numbers in the DBA object, matched by name (they follow the order added above).
contrast_info <- dba.show(dba_obj, bContrasts = TRUE)

all_tables <- list()
summary_rows <- list()
pdf(file.path(outdir, "plots.pdf"), width = 8, height = 7)
dba.plotPCA(dba_obj, attributes = DBA_CONDITION, label = DBA_REPLICATE, score = DBA_SCORE_NORMALIZED)
dba.plotHeatmap(dba_obj, score = DBA_SCORE_NORMALIZED)
for (cond in contrasts) {
    i <- which(contrast_info$Group == cond & contrast_info$Group2 == reference)
    if (length(i) != 1) stop("Could not find contrast ", cond, " vs ", reference)
    rep <- dba.report(dba_obj, contrast = i, th = 1, bCounts = FALSE)
    tab <- as.data.frame(rep)
    tab <- tab[, setdiff(colnames(tab), "strand")]
    names(tab)[names(tab) == "seqnames"] <- "chr"
    tab <- tab[order(tab$FDR), ]
    name <- paste0(cond, "_vs_", reference)
    write.table(tab, file.path(outdir, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    # Stacked table: per-contrast column names (Conc_4h, Conc_noUV, ...) -> generic names.
    stacked <- tab
    conc_cols <- grep("^Conc_", names(stacked))
    names(stacked)[conc_cols] <- c("Conc_group", "Conc_reference")
    all_tables[[name]] <- cbind(contrast = name, stacked)
    sig <- tab$FDR < fdr
    summary_rows[[name]] <- data.frame(
        contrast = name,
        peaks_tested = nrow(tab),
        significant = sum(sig),
        more_accessible = sum(sig & tab$Fold > 0),
        less_accessible = sum(sig & tab$Fold < 0)
    )
    message(name, ": ", sum(sig), " peaks with FDR < ", fdr)
    dba.plotMA(dba_obj, contrast = i, th = fdr)
    title(sub = name)
}
invisible(dev.off())

write.table(do.call(rbind, all_tables), snakemake@output[["results"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, summary_rows), snakemake@output[["summary"]], sep = "\t", quote = FALSE, row.names = FALSE)

normalized <- as.data.frame(dba.peakset(dba_obj, bRetrieve = TRUE, score = DBA_SCORE_NORMALIZED, DataType = DBA_DATA_FRAME))
write.table(normalized, snakemake@output[["normalized"]], sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(dba_obj, snakemake@output[["rds"]])
