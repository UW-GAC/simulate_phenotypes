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
        h2 <- heritability(G, beta, varComp)
    } else {
        stop("h2 or beta must be specified")
    }
    
    Gbeta <- G * beta
    return(list(h2=h2, beta=beta, Gbeta=Gbeta))
}



#'  Calculate heritability from beta
#'
#' @param G genotype vector
#' @param beta effect size
#' @param varComp 2-element vector with (genetic, error) variance components
#' @return h2 (heritability)
#' @importFrom stats var
#' @export
heritability <- function(G, beta, varComp) {
    h2 <- (var(G)*beta^2) / (var(G)*beta^2 + sum(varComp))
    return(h2)
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
#' @param G genotype matrix with sample.id as rownames and variant.id as colnames
#' @param h2 heritability
#' @param beta effect size for variant
#' @param varComp 2-element vector with (genetic, error) variance components
#' @param dat AnnotatedDataFrame with sample.id, outcome, and covariates
#' @param outcome character string specifying the name of the outcome variable in \code{dat}
#' @param cov.mat covariance matrix with sample.id as rownames
#' @param strata named list of sample.id in strata
#' @param covars A vector of character strings specifying the names of the fixed effect covariates in \code{dat}
#' @param power.signif p-value threshold for significance in power calculations
#' @return association test results for outcome and variant
#' @export
variant_assoc <- function(G, h2=NULL, beta=NULL, varComp, dat, outcome, cov.mat, strata=NULL, covars=NULL, power.signif=5e-9) {
    if (!requireNamespace("GENESIS")) {
        stop("must install GENESIS to use variant_assoc")
    }
    
    if (is.null(h2) & is.null(beta)) {
        stop("one of h2 or beta must be specified")
    }
    
    if (is.null(strata)) {
        strata <- rownames(G)
    }
    if (!is.list(strata)) {
        strata <- list(all=strata)
        if (!is.null(h2)) h2 <- list(all=h2)
        if (!is.null(beta)) beta <- list(all=beta)
    }
    if (!(setequal(names(strata), names(h2)) | setequal(names(strata), names(beta)))) {
        stop("h2 or beta must match strata")
    }
    
    dat$sample.id <- as.character(dat$sample.id)
    rownames(dat) <- dat$sample.id
    
    if (nrow(dat) < nrow(cov.mat)) {
        sel <- rownames(cov.mat) %in% dat$sample.id
        cov.mat <- cov.mat[sel,sel]
    }
    
    variant.id <- colnames(G)
    if (is.null(variant.id)) variant.id <- 1:ncol(G)
    res <- list()
    for (grp in names(strata)) {
        samp <- as.character(strata[[grp]])
        N <- length(samp)
        
        G.grp <- G[as.character(samp),,drop=FALSE]
        freq <- 0.5*colMeans(G.grp, na.rm=TRUE)
        res.grp <- list()
        for (i in 1:ncol(G.grp)) {
            if (is.null(h2)) {
                eff <- variant_effect(G=as.vector(G.grp[,i]), beta=beta[[grp]], varComp=varComp)
            } else {
                eff <- variant_effect(G=as.vector(G.grp[,i]), h2=h2[[grp]], varComp=varComp)
            }
            
            tmp <- pData(dat)
            tmp[samp, outcome] <- tmp[samp, outcome] + eff$Gbeta
            pData(dat) <- tmp
            res.grp[[i]] <- list(group=grp, N=N, beta=eff$beta, h2=eff$h2, power=power(N, eff$h2, pval=power.signif))
        }
        res.grp <- as.data.frame(data.table::rbindlist(res.grp))
        
        nullmod <- GENESIS::fitNullModel(dat, outcome=outcome, covars=covars, cov.mat=cov.mat, 
                                         sample.id=samp, start=varComp, verbose=FALSE)
        assoc <- GENESIS:::testGenoSingleVar(nullmod, G[as.character(nullmod$sample.id),])
        res[[grp]] <- cbind(variant.id, res.grp, freq, assoc, stringsAsFactors=FALSE)
    }
    
    if (length(strata) > 1) {
        grp <- "pooled"
        samp <- as.character(unlist(strata))
        N <- length(samp)
        G.grp <- G[as.character(samp),,drop=FALSE]
        freq <- 0.5*colMeans(G.grp, na.rm=TRUE)
        beta.pooled <- sum(sapply(names(strata), function(grp) length(strata[[grp]])*beta[[grp]])) / N
        h2.pooled <- apply(G.grp, 2, heritability, beta.pooled, varComp)
        power.pooled <- power(N, h2.pooled, pval=power.signif)
        res.grp <- list(group=grp, N=N, beta=beta.pooled, h2=h2.pooled, power=power.pooled)
        
        nullmod <- GENESIS::fitNullModel(dat, outcome=outcome, covars=covars, cov.mat=cov.mat, 
                                         sample.id=samp, start=varComp, verbose=FALSE)
        assoc <- GENESIS:::testGenoSingleVar(nullmod, G[as.character(nullmod$sample.id),])
        res[[grp]] <- cbind(variant.id, res.grp, freq, assoc, stringsAsFactors=FALSE)
    }
    
    res <- as.data.frame(data.table::rbindlist(res))
    return(res)
}



#' Return genotypes for set of samples and variants
#' 
#' @param gdsobj \code{\link{SeqVarGDSClass}} object
#' @param variant.sel indices of variants to return
#' @param variant.id ids of variants to return
#' @param sample.id ids of samples to return
#' @return matrix of genotypes
#' @export
variant_genotypes <- function(gdsobj, variant.sel, variant.id, sample.id=NULL) {
    if (!requireNamespace("SeqArray")) {
        stop("must install SeqArray and GENESIS to use variant_genotypes")
    }
    
    if (missing(variant.sel) & missing(variant.id)) {
        stop("one of variant.sel or variant.id must be specified")
    }
    
    if (!missing(variant.sel)) {
        SeqArray::seqSetFilter(gdsobj, variant.sel=variant.sel, sample.id=sample.id, verbose=FALSE)
    } else {
        SeqArray::seqSetFilter(gdsobj, variant.id=variant.id, sample.id=sample.id, verbose=FALSE)
    }
    
    geno <- SeqArray::seqGetData(gdsobj, "$dosage_alt")
    rownames(geno) <- SeqArray::seqGetData(gdsobj, "sample.id")
    colnames(geno) <- SeqArray::seqGetData(gdsobj, "variant.id")
    return(geno)
}
