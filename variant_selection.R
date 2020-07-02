var_multiple_strata <- function(a, b, min.freq=0.1) {
    which(a > min.freq & b > min.freq)
}

var_single_stratum <- function(a, b, min.freq=0.1) {
  which(a > min.freq & b == 0)
}
