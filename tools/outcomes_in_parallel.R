library(doParallel)
library(parallel)
library(foreach)
library(argparser)
#remotes::install_github("UW-GAC/simulate_phenotypes/simphen", upgrade=FALSE)

argp <- arg_parser("generate random outcomes")
argp <- add_argument(argp, "covmat_file", help="RData file with objects covmat and blocks")
argp <- add_argument(argp, "out_file", help="out file")
argp <- add_argument(argp, "--varComp1", help="variance component 1", type="integer")
argp <- add_argument(argp, "--varComp2", help="variance component 2", type="integer")
argp <- add_argument(argp, "--num_outcomes", help="number of outcomes to generate", type="integer")
argv <- parse_args(argp)

# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores())

# Activate cluster for foreach library
doParallel::registerDoParallel(cl)

# load covmat and blocks
load(argv$covmat_file)

varComp <- c(argv$varComp1, argv$varComp2)

outcomes <- foreach::foreach(i = 1:argv$num_outcomes) %dopar% {
    simphen:::outcome_from_covMat_block_indices(covmat, blocks, varComp)
}

# save output
saveRDS(outcomes, file=argv$out_file)

# Stop cluster to free up resources
parallel::stopCluster(cl)
