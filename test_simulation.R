library(GENESIS)
source("simulation_functions.R")

## test matrix (this is a Matrix object)
covMat <- GWASTools::getobj("/Users/stephanie/GCC_code/analysis_pipeline/testdata/pcrelate_Matrix.RData")
 
varComp <- c(40, 100) # (genetic, error) - from height?
outcome <- outcome_from_covMat(covMat, varComp)

## fit the null model and see if variance is what I expect
dat <- data.frame(y=outcome, row.names=rownames(covMat))
nullmod <- fitNullModel(dat, outcome="y", cov.mat=covMat)
nullmod$varComp


## test with grm
grm <- GWASTools::getobj("/Users/stephanie/GCC_code/analysis_pipeline/testdata/grm.RData")
covMat <- grm$grm
dimnames(covMat) <- list(grm$sample.id, grm$sample.id)

varComp <- c(40, 100) # (genetic, error) - from height?
outcome <- outcome_from_covMat(covMat, varComp)

dat <- data.frame(y=outcome, row.names=rownames(covMat))
nullmod <- fitNullModel(dat, outcome="y", cov.mat=covMat)
nullmod$varComp
