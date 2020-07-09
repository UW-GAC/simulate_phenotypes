#' Add effect of a variant to an outcome
#'
#' @param outcome vector
#' @param G genotype vector
#' @param h2 heritability
#' @return outcome with variant effect
#' @importFrom stats var
#' @export
add_variant_effect <- function(outcome, G, h2) {
    Ve <- var(outcome)
    alpha <- sqrt((h2*Ve)/((1-h2)*var(G)))
    alpha[alpha == Inf] <- 0
    snpeff <- G * alpha
    return(outcome + snpeff)
}
