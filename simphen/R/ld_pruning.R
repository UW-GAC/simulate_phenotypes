#' LD prune variants in all strata
#' 
#' @param gdsobj \code{\link[SeqArray]{SeqVarGDSClass}} object
#' @param strata named list of sample.id in strata
#' @param variant.id ids of variants to consider
#' @param maf MAF threshold for pruning in pooled sample
#' @param ld.threshold LD pruning threshold: single value or named list matching \code{strata}
#' @param seed seed for LD pruning: single value or named list matching \code{strata}
#' @param verbose logical for whether to print messages
#' @return vector with variant.id
#' @export
ldprune_strata <- function(gdsobj, strata, variant.id=NULL,
                           maf=0.01, ld.threshold=sqrt(0.1),
                           seed=NULL, verbose=FALSE) {
    if (!requireNamespace("SNPRelate")) {
        stop("must install SNPRelate to use ldprune_strata")
    }
    
    n.grps <- length(strata)
    grps <- names(strata)
    if (!is.list(ld.threshold)) {
        ld.threshold = as.list(rep(ld.threshold, n.grps+1))
        names(ld.threshold) <- c("all", grps)
    }
    if (!is.null(seed) && !is.list(seed)) {
        seed = as.list(rep(seed, n.grps+1))
        names(seed) <- c("all", grps)
    }
    
    if (!is.null(seed)) set.seed(seed[["all"]])
    pruned.all <- SNPRelate::snpgdsLDpruning(gdsobj, sample.id=unlist(strata), snp.id=variant.id,
                                  maf=maf, method="corr", ld.threshold=ld.threshold[["all"]], 
                                  verbose=verbose)
    pruned <- unlist(pruned.all, use.names=FALSE)
    
    for (grp in names(strata)) {
        samp <- as.character(strata[[grp]])
        
        if (!is.null(seed)) set.seed(seed[[grp]])
        pruned.grp <- SNPRelate::snpgdsLDpruning(gdsobj, sample.id=strata[[grp]], snp.id=pruned,
                                      remove.monosnp=FALSE, method="corr", 
                                      ld.threshold=ld.threshold[[grp]], 
                                      verbose=verbose)
        pruned <- unlist(pruned.grp, use.names=FALSE)
    }
    return(pruned)
}