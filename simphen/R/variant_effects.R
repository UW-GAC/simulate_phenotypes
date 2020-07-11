#' Add effect of a variant to an outcome
#'
#' @param G genotype vector
#' @param h2 heritability
#' @param varComp 2-element vector with (genetic, error) variance components
#' @return variant effect
#' @importFrom stats var
#' @export
variant_effect <- function(G, h2, varComp) {
    htot2 <- varComp[1]/sum(varComp)
    alpha <- sqrt((h2*varComp[2])/((1-htot2)*var(G)))
    alpha[alpha == Inf] <- 0
    snpeff <- G * alpha
    return(snpeff)
}
