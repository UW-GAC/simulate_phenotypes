library(argparser)

argp <- arg_parser("LD prune to select variants in multiple groups")
argp <- add_argument(argp, "--gds_file", help="GDS file")
argp <- add_argument(argp, "--strata_file", help="RDS file with list of sample.id")
argp <- add_argument(argp, "--variant_file", help="RDS file with list of variant.id")
argp <- add_argument(argp, "--maf", help="MAF threshold for pruning in pooled sample", default=0.01)
argp <- add_argument(argp, "--ld_threshold", "LD pruning |r| threshold", default=sqrt(0.1))
argp <- add_argument(argp, "--out_file", help="out file")
argv <- parse_args(argp)
print(argv)

library(SeqArray)
#remotes::install_github("UW-GAC/simulate_phenotypes/simphen", upgrade=FALSE)
library(simphen)
sessionInfo()


strata <- readRDS(argv$strata_file)

variant.id <- readRDS(argv$variant_file)

gds <- seqOpen(argv$gds_file)

var.id <- ldprune_strata(gds, strata, variant.id,
                         maf=argv$maf, ld.threshold=argv$ld_threshold,
                         verbose=TRUE)

saveRDS(var.id, file=argv$out_file)

seqClose(gds)
