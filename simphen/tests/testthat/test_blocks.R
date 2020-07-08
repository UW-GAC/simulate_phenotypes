context("block tests")

test_that("block indices", {
    x <- diag(nrow=10, ncol=10)
    x[2:4,2:4] <- 0.5
    x[7:9,7:9] <- 0.5
    diag(x) <- 1
    bdMatrix <- Matrix(x)
    blocks <- block_indices(bdMatrix)
    expect_true(all(bdMatrix == do.call(bdiag, lapply(blocks, function(ind) bdMatrix[ind,ind]))))
})


test_that("outcome from blocks matches original", {
    data("1KG_pcrelate_blockDiag_Matrix")
    covMat <- blockDiagMatrix1KG[1:500,1:500]
    varComp <- c(40, 100)
    set.seed(80)
    outcome <- rnorm(nrow(covMat))
    outcome.all <- outcome_from_covMat(covMat, varComp, outcome=outcome)
    outcome.blocks <- outcome_from_covMat_blocks(covMat, varComp, outcome=outcome)
    expect_equal(outcome.all, outcome.blocks)
})
