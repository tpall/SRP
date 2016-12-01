
#' Stable retrospective power
#'
#' \code{srp} uses p-values and false discovery rate (FDR) to calculate stable retrospective power.
#'
#' @param pvalues A vector of raw p-values.
#' @param FDR A level at which to control the false discovery rate (FDR). Deafults to 0.05.
#' @param ... Additional arguments passed to \code{\link{qvalue}}.
#' @return A numeric that is the estimated stable retrospective power.
#' @seealso \code{\link{qvalue}} how q-values and pi0 are calculated and for arguments.
#' @export
#' @examples
#' # import data
#' library(qvalue)
#' data("hedenfalk")
#' pvalues <- hedenfalk$p
#'
#' # calculate SRP
#' pw <- srp(pvalues, FDR = 0.05)
#' pw
#'

srp <- function(pvalues, FDR = 0.05, ...){

  qobj <- qvalue::qvalue(pvalues, ...)
  qvalues <- qobj$qvalues
  pi0 <- qobj$pi0

  qvalues <- qvalues[!is.na(qvalues)]
  q <- sum(qvalues <= FDR)
  n_tests <- length(qvalues)

  if(q == 0){
    warning("No significant effects", call. = FALSE)
    pw <- 0
    return(pw)
  }

  pw <- ((1 - FDR) * q) / ((1 - pi0) * n_tests)

  if(pw > 1){
    warning("Infinite power", call. = FALSE)
    pw <- 1
    return(pw)
  }

  return(pw)
}


