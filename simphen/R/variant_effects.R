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



#' Run association test for an outcome and a set of variants
#'
#' Run association test for an outcome and a set of variants
#'
#' Either h2 or beta must be specified.
#'
#' @param G genotype matrix with sample.id as rownames
#' @param h2 heritability
#' @param beta effect size for variant
#' @param varComp 2-element vector with (genetic, error) variance components
#' @param dat AnnotatedDataFrame with sample.id, outcome, and covariates
#' @param outcome character string specifying the name of the outcome variable in \code{dat}
#' @param cov.mat covariance matrix with sample.id as rownames
#' @param covars A vector of character strings specifying the names of the fixed effect covariates in \code{dat}
#' @return association test results for outcome and variant
#' @export
variant_assoc <- function(G, h2=NULL, beta=NULL, varComp, dat, outcome, cov.mat, covars=NULL) {
    if (!requireNamespace("GENESIS")) {
        stop("must install GENESIS to use variant_assoc")
    }
    
    if (is.null(h2) & is.null(beta)) stop("one of h2 or beta must be specified")
    
    G <- G[as.character(dat$sample.id),,drop=FALSE]
    res <- list()
    for (i in 1:ncol(G)) {
        if (is.null(h2)) {
            eff <- variant_effect(G=as.vector(G[,i]), beta=beta, varComp=varComp)
        } else {
            eff <- variant_effect(G=as.vector(G[,i]), h2=h2, varComp=varComp)
        }
        dat[[outcome]] <- dat[[outcome]] + eff$Gbeta
        res[[i]] <- eff[c("beta", "h2")]
    }
    res <- as.data.frame(data.table::rbindlist(res))

    if (nrow(dat) < nrow(cov.mat)) {
        sel <- rownames(cov.mat) %in% dat$sample.id
        cov.mat <- cov.mat[sel,sel]
    }

    nullmod <- GENESIS::fitNullModel(dat, outcome=outcome, covars=covars, cov.mat=cov.mat, start=varComp, verbose=FALSE)
    assoc <- GENESIS:::testGenoSingleVar(nullmod, G[as.character(nullmod$sample.id),])
    assoc <- cbind(res, assoc)

    return(assoc)
}



#' Return genotypes for set of samples and variants
#' 
variant_genotypes <- function(gdsobj, variant.id, sample.id=NULL) {
    if (!requireNamespace("SeqArray")) {
        stop("must install SeqArray and GENESIS to use variant_genotypes")
    }
    
    #SeqArray::seqSetFilter(gdsobj, variant.sel=variant.sel, sample.id=sample.id, verbose=FALSE)
    SeqArray::seqSetFilter(gdsobj, variant.id=variant.id, sample.id=sample.id, verbose=FALSE)
    
    geno <- SeqArray::seqGetData(gdsobj, "$dosage_alt")
    rownames(geno) <- SeqArray::seqGetData(gdsobj, "sample.id")
    return(geno)
}
