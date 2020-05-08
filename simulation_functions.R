library(Matrix)
library(data.table)

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
#' only works if all block elements are nonzero
#' https://stackoverflow.com/questions/54472962/find-the-indices-for-the-sub-matrices-in-block-matrix
#' 
#' @param x Block diagonal matrix
#' @return list of indices of matrix blocks
block_indices_nonzero <- function(x) {
    # from stackoverflow
    #ind <- split(seq_len(nrow(x)), max.col(abs(x) > 0, "first"))
    
    #using sparse Matrix representation
    row.counts <- table(x@i + 1)
    r <- 1
    b <- 1
    ind <- list()
    while (r < length(row.counts)) {
        n <- r + row.counts[r] - 1
        ind[[b]] <- r:n
        r <- n + 1
        b <- b + 1
    }
    return(ind)
}

#' Return list of blocks in a block-diagonal matrix
#' 
#' Code taken from GENESIS::makeSparseMatrix
#' 
#' @param x Block diagonal matrix
#' @return list of indices of matrix blocks
#' @import data.table
#' @import igraph
block_indices <- function(x) {
    
    # get the table of all related pairs
    rel <- apply(x, MARGIN = 1, FUN = function(v){ which(v > 0) })
    rel <- lapply(seq_along(rel), function(i) { data.table('ID1' = i, 'ID2' = rel[[i]]) })
    rel <- do.call(rbind, rel)
    setkeyv(rel, c('ID1', 'ID2'))
    
    # create graph of relatives
    g <- igraph::graph_from_data_frame(rel[ID1 != ID2])
    # extract cluster membership
    clu <- igraph::components(g)
    mem <- clu$membership
    
    blocks <- list()
    for(i in 1:clu$no){
        # samples in the cluster
        blocks[[i]] <- unname(which(mem == i))
    }
    
    return(blocks)
}


outcome_from_covMat_blocks <- function(covMat, varComp, outcome=NULL) {
    
    if (is.null(outcome)) {
        outcome <- rnorm(nrow(covMat))
    } else {
        stopifnot(length(outcome) == nrow(covMat))
    }
    
    outcome.list <- lapply(block_indices(covMat), function(ind) {
        if (length(ind) == 1) return(outcome[ind])
        outcome_from_covMat(covMat[ind,ind], varComp, outcome=outcome[ind])
    })
    return(unlist(outcome.list))
}
