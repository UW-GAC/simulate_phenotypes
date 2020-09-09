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
        h2 <- (1 - varComp[1]/sum(varComp)) / (varComp[2]/(var(G)*beta^2) + varComp[2]/sum(varComp))
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
#' @param beta effect size
#' @param pval p-value threshold for significance
#' @return power
#' @references https://github.com/kaustubhad/gwas-power/blob/master/power_calc_functions.R
#' @export
power <- function(N, h2, beta, pval=5e-8) {
    # Significance threshold for chi-square, corresponding to P-value threshold
    th <- qchisq(pval, df=1, lower.tail=FALSE)

    # non-centrality parameter
    ncp <- N*h2/(1-h2)

    # power
    pow <- pchisq(th, df=1, lower.tail=FALSE, ncp=ncp)

    return(pow)
}
