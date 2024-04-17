library(argparser)

argp <- arg_parser("test variant effects")
argp <- add_argument(argp, "--gds_file", help="GDS file")
argp <- add_argument(argp, "--covmat_file", help="RDS file with block-diagonal covariance matrix")
argp <- add_argument(argp, "--variant_file", help="RDS file with list of variant.id")
argp <- add_argument(argp, "--h2", help="heritability")
argp <- add_argument(argp, "--beta", help="beta")
argp <- add_argument(argp, "--varComp1", help="variance component 1", type="integer")
argp <- add_argument(argp, "--varComp2", help="variance component 2", type="integer")
argp <- add_argument(argp, "--seed", help="seed for simulating outcome")
argp <- add_argument(argp, "--out_file", help="out file")
argv <- parse_args(argp)
print(argv)

library(foreach)
library(SeqArray)
library(simphen)
sessionInfo()

# load covmat
covmat <- readRDS(argv$covmat_file)

# load variant.id and h2 or beta
variant.id <- readRDS(argv$variant_file)
if ("h2" %in% names(variant.id)) {
    h2 <- variant.id$h2
    variant.id <- variant.id$variant.id
    input <- "h2"
} else if ("beta" %in% names(variant.id)) {
    beta <- variant.id$beta
    variant.id <- variant.id$variant.id
    input <- "beta"
} else if (!is.na(argv$h2)) {
    stopifnot(!is.vector(variant.id))
    h2 <- rep(as.numeric(h2), length(variant.id))
    input <- "h2"
} else if (!is.na(argv$beta)) {
    stopifnot(!is.vector(variant.id))
    beta <- rep(as.numeric(beta), length(variant.id))
    input <- "beta"
} else {
    stop("must specify either h2 or beta")
}
nvar <- length(variant.id)

varComp <- c(argv$varComp1, argv$varComp2)

# simulate outcome
message("Simulating outcome")
seed <- if (is.na(argv$seed)) Sys.time() else as.integer(argv$seed)
set.seed(seed)
outcome <- outcome_from_covMat_blocks(covmat, varComp)

# read genotypes
message("Reading genotypes")
gds <- seqOpen(argv$gds_file)
geno <- variant_genotypes(gds, variant.id=variant.id)
stopifnot(setequal(colnames(geno), variant.id))
seqClose(gds)

# order outcome to match GDS (was block-diagonal order)
sample.id <- rownames(geno)
outcome <- outcome[sample.id]

message("Calculating variant effects")
eff <- foreach::foreach(i = 1:nvar, .combine='+') %do% {
    var_ind <- as.character(variant.id[i])
    if (input == "h2") {
        tmp <- variant_effect(G=geno[,var_ind], h2=h2[i], varComp=varComp)
    } else {
        tmp <- variant_effect(G=geno[,var_ind], beta=beta[i], varComp=varComp)
    }
    tmp$Gbeta
}
outcome <- outcome + eff


# save output
saveRDS(outcome, file=argv$out_file)

