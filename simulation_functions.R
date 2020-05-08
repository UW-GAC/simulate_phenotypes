library(Matrix)

#' generate random normal outcome from covariance matrix
#'
#' @param covMat a covariance matrix
#' @param varComp 2-element vector with (genetic, error) variance components
#' @return numeric vector of outcomes
#' 
#' @import Matrix
outcome_from_covMat <- function(covMat, varComp, outcome=NULL) {
    
    if (is.null(outcome)) {
        outcome <- rnorm(nrow(covMat))
    } else {
        stopifnot(length(outcome) == nrow(covMat))
    }
    
    ## could add heterogeneous residual variance here
    covMat <- covMat*varComp[1] + Diagonal(nrow(covMat))*varComp[2]

    ## singular value decomposition
    svd.covMat <- svd(covMat)
    sqrt.covMat <- svd.covMat$v %*% diag(sqrt(svd.covMat$d)) %*% t(svd.covMat$v)

    ## svd.covMat$v - eigenvectors
    ## svd.covMat$d - eigenvalues (diagonals)

    ## sanity check:
    ## svd.covMat$v %*% diag(svd.covMat$d) %*% t(svd.covMat$v) should give original matrix
    ## chk <- svd.covMat$v %*% diag(svd.covMat$d) %*% t(svd.covMat$v)
    ## diff <- abs(chk - covMat)
    ## c(min(diff), max(diff))

    ## random normal outcome with covariance of the matrix
    outcome <- as.vector(sqrt.covMat %*% outcome)

    return(outcome)
}



#' Return list of blocks in a block-diagonal matrix
#' 
#' https://stackoverflow.com/questions/54472962/find-the-indices-for-the-sub-matrices-in-block-matrix
#' 
#' @param x Block diagonal matrix
#' @return list of indices of matrix blocks
block_indices <- function(x) {
    split(seq_len(nrow(x)), max.col(abs(x) > 0, "first"))
}


outcome_from_covMat_blocks <- function(covMat, varComp, outcome=NULL) {
    
    if (is.null(outcome)) {
        outcome <- rnorm(nrow(covMat))
    } else {
        stopifnot(length(outcome) == nrow(covMat))
    }
    
    outcome <- rnorm(nrow(covMat))
    outcome.list <- lapply(block_indices(covMat), function(ind) {
        if (length(ind) == 1) return(outcome[ind])
        outcome_from_covMat(covMat[ind,ind], varComp, outcome=outcome[ind])
    })
    return(unlist(outcome.list))
}
