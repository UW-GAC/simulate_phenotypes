% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_effects.R
\name{variant_effect}
\alias{variant_effect}
\title{Add effect of a variant to an outcome}
\usage{
variant_effect(G, h2, beta, varComp)
}
\arguments{
\item{G}{genotype vector}

\item{h2}{heritability}

\item{beta}{effect size}

\item{varComp}{2-element vector with (genetic, error) variance components}
}
\value{
A list with the variant effect and related parameters
\itemize{
  \item h2 - heritability
  \item beta - effect size
  \item Gbeta - vector of G * beta
}
}
\description{
Add effect of a variant to an outcome
}
