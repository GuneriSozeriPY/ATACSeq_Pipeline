# Differential accessibility with plain DESeq2 for one cell line (run via Snakemake `script:`).
# Runs next to DiffBind as an independent check: different counting (featureCounts, full-width
# consensus peaks) and DESeq2's own median-of-ratios size factors.
#
# Design: ~ replicate + condition (replicate series differ, e.g. HeLa R3 has lower TSS enrichment).
# Every condition is compared to the reference condition (config: differential.reference).
#
# Outputs (in the cell line's results/differential/<cell>_<build>/deseq2 folder)
#   <condition>_vs_<reference>.tsv  all tested peaks for one contrast, sorted by padj;
#                                   log2FoldChange = MLE, log2FC_shrunk = apeglm-shrunken (for ranking/plots)
#   all_contrasts.tsv               the same tables stacked, with a `contrast` column
#   summary.tsv                     peaks with padj < fdr: total / more / less accessible
#   normalized_counts.tsv           size-factor-normalized counts
#   plots.pdf                       PCA (VST; raw and replicate effect removed), MA plot per contrast
#   dds.rds                         the fitted DESeqDataSet

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
    library(ggplot2)
    library(BiocParallel)
})

p <- snakemake@params
samples <- unlist(p[["samples"]])
conditions <- unlist(p[["conditions"]])
replicates <- unlist(p[["replicates"]])
reference <- p[["reference"]]
fdr <- as.numeric(p[["fdr"]])
min_count <- as.integer(p[["min_count"]])
outdir <- p[["outdir"]]
register(MulticoreParam(as.integer(snakemake@threads)))

# featureCounts: 6 annotation columns, then one column per BAM (same order as `samples`).
fc <- read.delim(snakemake@input[["counts"]], comment.char = "#", check.names = FALSE)
peaks <- fc[, c("Geneid", "Chr", "Start", "End")]
colnames(peaks) <- c("peak", "chr", "start", "end")
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid
colnames(counts) <- samples

coldata <- data.frame(
    row.names = samples,
    # Reference first, then sample-sheet order (keeps contrasts in time order, not alphabetical).
    condition = factor(conditions, levels = c(reference, setdiff(unique(conditions), reference))),
    replicate = factor(replicates)
)
print(coldata)

design <- if (nlevels(coldata$replicate) > 1) ~ replicate + condition else ~condition
message("Design: ", deparse(design))

dds <- DESeqDataSetFromMatrix(counts, coldata, design)
dds <- dds[rowSums(counts(dds) >= min_count) >= 2, ]
message(nrow(dds), " of ", nrow(counts), " peaks kept (>= ", min_count, " fragments in >= 2 samples)")
dds <- DESeq(dds, parallel = TRUE)
print(round(sizeFactors(dds), 3))

norm <- counts(dds, normalized = TRUE)
write.table(cbind(peak = rownames(norm), round(norm, 2)), snakemake@output[["normalized"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(snakemake@output[["plots"]], width = 7, height = 5)
vsd <- vst(dds, blind = FALSE)
pca_plot <- function(v, title) {
    d <- plotPCA(v, intgroup = c("condition", "replicate"), returnData = TRUE)
    pv <- round(100 * attr(d, "percentVar"))
    ggplot(d, aes(PC1, PC2, colour = condition, shape = replicate)) +
        geom_point(size = 3) +
        labs(title = title, x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
        theme_bw()
}
print(pca_plot(vsd, "VST counts (top 500 peaks)"))
if (nlevels(coldata$replicate) > 1) {
    vsd_corrected <- vsd
    assay(vsd_corrected) <- removeBatchEffect(assay(vsd), batch = vsd$replicate,
                                              design = model.matrix(~condition, colData(vsd)))
    print(pca_plot(vsd_corrected, "VST counts, replicate effect removed"))
}

contrasts <- setdiff(levels(coldata$condition), reference)
all_tables <- list()
summary_rows <- list()
for (cond in contrasts) {
    name <- paste0(cond, "_vs_", reference)
    res <- results(dds, contrast = c("condition", cond, reference), alpha = fdr)
    shrunk <- lfcShrink(dds, coef = paste0("condition_", name), type = "apeglm", quiet = TRUE)
    res$log2FC_shrunk <- shrunk$log2FoldChange
    tab <- merge(peaks, as.data.frame(res), by.x = "peak", by.y = "row.names")
    tab <- tab[order(tab$padj, na.last = TRUE), ]
    write.table(tab, file.path(outdir, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    all_tables[[name]] <- cbind(contrast = name, tab)
    sig <- !is.na(tab$padj) & tab$padj < fdr
    summary_rows[[name]] <- data.frame(
        contrast = name,
        peaks_tested = sum(!is.na(tab$padj)),
        significant = sum(sig),
        more_accessible = sum(sig & tab$log2FoldChange > 0),
        less_accessible = sum(sig & tab$log2FoldChange < 0)
    )
    message(name, ": ", sum(sig), " peaks with padj < ", fdr)
    plotMA(shrunk, alpha = fdr, main = paste(name, "(apeglm-shrunken LFC)"))
}
invisible(dev.off())

write.table(do.call(rbind, all_tables), snakemake@output[["results"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, summary_rows), snakemake@output[["summary"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(dds, snakemake@output[["rds"]])
