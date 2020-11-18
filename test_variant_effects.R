library(argparser)

argp <- arg_parser("test variant effects")
argp <- add_argument(argp, "--gds_file", help="GDS file")
argp <- add_argument(argp, "--covmat_file", help="RData file with objects covmat and blocks")
argp <- add_argument(argp, "--covar_file", help="RDS file with AnnotatedDataFrame containing sample.id and covariates")
argp <- add_argument(argp, "--outcome_file", help="RDS file with random outcomes")
argp <- add_argument(argp, "--strata_file", help="RDS file with list of sample.id")
argp <- add_argument(argp, "--variant_file", help="RDS file with list of variant.id")
argp <- add_argument(argp, "--h2", help="heritability for each stratum", nargs=Inf)
argp <- add_argument(argp, "--beta", help="beta for each stratum", nargs=Inf)
argp <- add_argument(argp, "--varComp1", help="variance component 1", type="integer")
argp <- add_argument(argp, "--varComp2", help="variance component 2", type="integer")
argp <- add_argument(argp, "--num_variants", help="number of variants to test", type="integer")
argp <- add_argument(argp, "--variant_block_size", help="number of variants to add to each outcome at the same time", type="integer", default=10)
argp <- add_argument(argp, "--out_file", help="out file")
argv <- parse_args(argp)
print(argv)

library(doParallel)
library(parallel)
library(foreach)
library(SeqArray)
library(Biobase)
remotes::install_github("UW-GAC/simulate_phenotypes/simphen", upgrade=FALSE)
library(simphen)
sessionInfo()

# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores())

# Activate cluster for foreach library
doParallel::registerDoParallel(cl)

# load outcomes
outcomes <- readRDS(argv$outcome_file)

# load covmat
load(argv$covmat_file)

# load covariates
dat <- readRDS(argv$covar_file)
covars <- setdiff(varLabels(dat), "sample.id")

# load strata
strata <- readRDS(argv$strata_file)

# load variant ids and sample
variant.id <- readRDS(argv$variant_file)
if (!is.na(argv$num_variants)) {
    variant.id <- variant.id[sort(sample(1:argv$num_variants))]
}
nvar <- length(variant.id)

varComp <- c(argv$varComp1, argv$varComp2)

# match beta/h2 to strata
if (!is.na(argv$h2)) {
    h2 <- setNames(as.numeric(argv$h2), names(strata))
} else if (!is.null(argv$beta)) {
    beta <- setNames(as.numeric(argv$beta), names(strata))
} else {
    stop("must specify either h2 or beta")
}

# read genotypes
gds <- seqOpen(argv$gds_file)
geno <- variant_genotypes(gds, variant.id=variant.id, sample.id=unlist(strata))
seqClose(gds)

n_iter <- ceiling(nvar / argv$variant_block_size)
var_blocks <- unname(split(variant.id, cut(1:nvar, n_iter)))

eff <- foreach::foreach(i = 1:n_iter) %dopar% {
    dat$outcome <- outcomes[[sample(length(outcomes), 1)]]
    var_ind <- as.character(var_blocks[[i]])
    if (!is.na(argv$h2)) {
        variant_assoc(geno[,var_ind,drop=FALSE], h2=h2, varComp=varComp,
                    dat=dat, outcome="outcome", cov.mat=covmat, 
                    strata=strata, covars=covars)
    } else {
        variant_assoc(geno[,var_ind,drop=FALSE], beta=beta, varComp=varComp,
                      dat=dat, outcome="outcome", cov.mat=covmat, 
                      strata=strata, covars=covars)
    }
}

eff.df <- as.data.frame(data.table::rbindlist(eff))

# save output
saveRDS(eff.df, file=argv$out_file)

# Stop cluster to free up resources
parallel::stopCluster(cl)
