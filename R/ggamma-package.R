
#' Generalized Gamma Probability Distribution
#'
#' Density, distribution function, quantile function and random
#' generation for the Generalized Gamma lifetime distributions.
#' 
#' @details
#' 
#' This package follows naming convention that is consistent with base R,
#' where density (or probability mass) functions, distribution functions,
#' quantile functions and random generation functions names are followed by
#' \code{d}, \code{p}, \code{q}, and \code{r} prefixes.
#' 
#' Behaviour of the functions is consistent with base R, where for
#' not valid parameters values \code{NaN}'s are returned, while
#' for values beyond function support \code{0}'s are returned
#' (e.g. for non-integers in discrete distributions, or for
#' negative values in functions with non-negative support).
#'
#' C++ was not used, as the R code proved itself most efficient.
#' See the package website page for more details.
#' 
#' @docType package
#' @name ggamma
#'
#' @importFrom stats dgamma pgamma qgamma rgamma
#' @importFrom stats runif
NULL
