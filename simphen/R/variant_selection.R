#' Select variants based on allele frequency
#'
#' \code{var_multiple_strata} selects variants present
#' in groups a and b with frequency \code{>min.freq}
#'
#' \code{var_single_stratum} selects variants
#' with frequency \code{>min.freq} in group a and
#' frequency 0 in group b
#'
#' @param a vector of variant frequency in group a
#' @param b vector of variant frequency in group b
#' @param min.freq miniumum frequency to be included in results
#' @return index of variant meeting criteria
#' @export
var_multiple_strata <- function(a, b, min.freq=0.1) {
    which(a > min.freq & b > min.freq)
}

#' @rdname var_multiple_strata
#' @export
var_single_stratum <- function(a, b, min.freq=0.1) {
  which(a > min.freq & b == 0)
}
