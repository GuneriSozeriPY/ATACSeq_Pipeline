#### Libraries ####
suppressPackageStartupMessages(library(ATACseqQC))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(Rsamtools))

library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(Rsamtools)
set.seed(1) 
args <- commandArgs(trailingOnly=T)

#### Rearrange ####

sample <- paste0("",args[1],"")

bamFile=paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/alignment/",sample,".blacklist-filtered.bam")  ### at AkD1n's dekstop
bamFileLabels <- sample

# bam_qc=bamQC(bamFile, outPath = NULL)
# bam_qc[1:10]

## files will be saved into outPath respective to the working directory
outPath <- paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample,"/")
dir.create(outPath)

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

bamTop100 <- scanBam(BamFile(bamFile, yieldSize = 100),
                     param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]

gal <- readBamFile(bamFile, tag=tags, asMates=TRUE, bigFile=TRUE)

shiftedBamFile <- file.path(outPath, paste0("",sample,".shifted.bam"))
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamFile)



 ### save the GRanges object for future use
saveRDS(gal1, file = paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample,"/",sample,"_gal1.rds"), ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
genome <- Hsapiens
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)
saveRDS(objs, file = paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample,"/",sample,"_objs.rds"), ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)

tsse <- TSSEscore(gal1, txs)
print(tsse$TSSEscore)

pdf(paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample,"/",sample,"_TSSEscore.pdf"))
plot(100*(-9:10-.5), tsse$values, type="b", col = "#264653",
     xlab="distance to TSS",
     ylab="aggregate TSS score",
     main=paste0("TSSE score for ",sample,""))
dev.off()


###Signal in NFR and Mononucleosome Fractions
bamFiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))

TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

librarySize <- estLibSize(bamFiles)

NTILE <- 101
dws <- ups <- 1010
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")],
                          TSS=TSS, seqlev = paste0("chr", c(1:22, "X", "Y")),
                          librarySize=librarySize,
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))


pdf(paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample,"/",sample,"_Heatmap_splitbam.pdf"))
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE) 
dev.off()

#Signal at TSS

out <- featureAlignedDistribution(sigs,
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l",
                                  ylab="Averaged coverage")

## rescale the nucleosome-free and nucleosome signals to 0~1 for plotting
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)

pdf(paste0("/cta/users/ygsuzeri/ATAC/NHF1_rep2/ATACseqQC/",sample, "/",sample, "_TSSprofile_splitbam.pdf"))
matplot(out, type="l", xaxt="n",
        xlab="Position (bp), TSS centered",
        ylab="Fraction of signal",
        main = paste0("",bamFileLabels,""))
axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()
