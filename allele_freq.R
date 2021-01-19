library(argparser)
library(SeqVarTools)

argp <- arg_parser("allele frequency")
argp <- add_argument(argp, "gds.file", help="GDS file")
argp <- add_argument(argp, "out.file", help="out file")
argp <- add_argument(argp, "--sample.file", help="sample file")
argp <- add_argument(argp, "--cpu", help="Number of CPUs to utilize for parallel processing", default=1)
argv <- parse_args(argp)

gds <- seqOpen(argv$gds.file)

if (!is.na(argv$sample.file)) {
    samp <- readRDS(argv$sample.file)
    seqSetFilter(gds, sample.id=samp)
}

ref.freq <- seqAlleleFreq(gds, parallel=argv$cpu)
alt.freq <- 1-ref.freq

var.info <- variantInfo(gds, alleles=FALSE)

res <- cbind(var.info, alt.freq)
saveRDS(res, file=argv$out.file)

seqClose(gds)
