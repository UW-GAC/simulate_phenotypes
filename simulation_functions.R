
#' generate random normal outcome from covariance matrix
#'
#' @param covMat a covariance matrix
#' @param varComp 2-element vector with (genetic, error) variance components
#' @return numeric vector of outcomes
outcome_from_covMat <- function(covMat, varComp) {
    
    ## could add heterogeneous residual variance here
    covMat <- covMat*varComp[1] + diag(nrow(covMat))*varComp[2]

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
    outcome <- rnorm(nrow(covMat))
    outcome <- as.vector(sqrt.covMat %*% outcome)

    return(outcome)
}
