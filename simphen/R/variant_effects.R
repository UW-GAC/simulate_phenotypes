#' Add effect of a variant to an outcome
#'
#' @param G genotype vector
#' @param h2 heritability
#' @param beta effect size
#' @param varComp 2-element vector with (genetic, error) variance components
#' @return A list with the variant effect and related parameters
#' \itemize{
#'   \item h2 - heritability
#'   \item beta - effect size
#'   \item Gbeta - vector of G * beta
#' }
#' @importFrom stats var
#' @export
variant_effect <- function(G, h2, beta, varComp) {
    if (!missing(h2) & !missing(beta)) {
         stop("only one of h2 and beta may be specified")
    }
    
    if (missing(beta)) {
        htot2 <- (varComp[1] + h2*varComp[2])/sum(varComp)
        beta <- sqrt((h2*varComp[2])/((1-htot2)*var(G)))
        beta[beta == Inf] <- 0
    } else if (missing(h2)) {
        h2 <- (var(G)*beta^2) / (var(G)*beta^2 + sum(varComp))
    } else {
        stop("h2 or beta must be specified")
    }
    
    Gbeta <- G * beta
    return(list(h2=h2, beta=beta, Gbeta=Gbeta))
}


#' Calculate power
#' 
#' @param N number of samples
#' @param h2 heritability
#' @param pval p-value threshold for significance
#' @return power
#' @references https://github.com/kaustubhad/gwas-power/blob/master/power_calc_functions.R
#' @importFrom stats qchisq pchisq
#' @export
power <- function(N, h2, pval=5e-9) {
    # Significance threshold for chi-square, corresponding to P-value threshold
    th <- qchisq(pval, df=1, lower.tail=FALSE)

    # non-centrality parameter
    ncp <- N*h2/(1-h2)

    # power
    pow <- pchisq(th, df=1, lower.tail=FALSE, ncp=ncp)

    return(pow)
}


#' Run association test for an outcome and a variant
#'
#' Run association test for an outcome and a variant
#'
#' Either h2 or beta must be specified.
#'
#' Order of samples in \code{outcome} and \code{covars} must match \code{cov.mat}. The rownames of \code{cov.mat} must correspond to the \code{sample.id} node in \code{gdsobj}.
#'
#' @param variant.sel index of variant in unfilted gdsobj
#' @param h2 heritability
#' @param beta effect size for variant
#' @param varComp 2-element vector with (genetic, error) variance components
#' @param gdsobj SeqVarGDSClass object
#' @param dat AnnotatedDataFrame with sample.id, outcome, and covariates
#' @param outcome character string specifying the name of the outcome variable in \code{dat}
#' @param cov.mat covariance matrix
#' @param covars A vector of character strings specifying the names of the fixed effect covariates in \code{dat}
#' @return association test results for outcome and variant
#' @export
variant_assoc <- function(variant.sel, h2=NULL, beta=NULL, varComp, gdsobj, dat, outcome, cov.mat, covars=NULL) {
    if (!requireNamespace("SeqArray") | !requireNamespace("GENESIS")) {
        stop("must install SeqArray and GENESIS to use variant_assoc")
    }
    
    if (is.null(h2) & is.null(beta)) stop("one of h2 or beta must be specified")
    
    # should we use variant.id, or variant.sel for speed?
    stopifnot(length(variant.sel) == 1)

    SeqArray::seqSetFilter(gdsobj, variant.sel=variant.sel, verbose=FALSE)
    
    geno <- SeqArray::seqGetData(gdsobj, "$dosage_alt")
    rownames(geno) <- SeqArray::seqGetData(gdsobj, "sample.id")
    geno <- geno[as.character(dat$sample.id),,drop=FALSE]
    if (is.null(h2)) {
        eff <- variant_effect(G=as.vector(geno), beta=beta, varComp=varComp)
    } else {
        eff <- variant_effect(G=as.vector(geno), h2=h2, varComp=varComp)
    }
    dat[[outcome]] <- dat[[outcome]] + eff$Gbeta

    if (nrow(dat) < nrow(cov.mat)) {
        sel <- rownames(cov.mat) %in% dat$sample.id
        cov.mat <- cov.mat[sel,sel]
    }

    nullmod <- GENESIS::fitNullModel(dat, outcome=outcome, covars=covars, cov.mat=cov.mat, start=varComp, verbose=FALSE)
    assoc <- GENESIS:::testGenoSingleVar(nullmod, geno[as.character(nullmod$sample.id),])
    assoc <- cbind(eff[c("beta", "h2")], assoc)

    return(assoc)
}
