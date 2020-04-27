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


## test with dense vs sparse matrix
dense <- GWASTools::getobj("/sbgenomics/project-files/1KG_pcrelate_dense_Matrix.RData")
sparse <- GWASTools::getobj("/sbgenomics/project-files/1KG_pcrelate_Matrix.RData")

lapply(list(dense=dense, sparse=sparse), function(covMat) {
  varComp <- c(40, 100)
  set.seed(100)
  outcome <- outcome_from_covMat(covMat, varComp)

  dat <- data.frame(y=outcome, row.names=rownames(covMat))
  nullmod <- fitNullModel(dat, outcome="y", cov.mat=covMat)
  nullmod$varComp
})


## test with TOPMed freeze 8
covMat <- GWASTools::getobj("/sbgenomics/project-files/pcrelate_kinshipMatrix_sparseDeg4_v2.RData")

## svd gives "Cholmod error 'problem too large'"
cov25k <- covMat[1:25000,1:25000]

varComp <- c(40, 100)
outcome <- outcome_from_covMat(cov25k, varComp)

dat <- data.frame(y=outcome, row.names=rownames(covMat))
nullmod <- fitNullModel(dat, outcome="y", cov.mat=covMat)
nullmod$varComp
