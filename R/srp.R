
#' Stable retrospective power.
#'
#' \code{srp} uses q-values, the proportion of true null p-values (pi0) and false discovery rate (FDR) to calculate stable retrospective power.
#'
#' @param qvalues A vector of q-values.
#' @param pi0 The proportion of true null p-values.
#' @param FDR A level at which to control the false discovery rate (FDR). Deafults to 0.05.
#' @return A numeric that is the estimated stable retrospective power.
#' @seealso \code{\link{qvalue}} to calculate q-values and pi0
#' @export
#' @examples
#' # import data
#' library(qvalue)
#' data("hedenfalk")
#' qobj <- qvalue(hedenfalk$p)
#' qvalues <- qobj$qvalues
#' pi0 <- qobj$pi0
#'
#' # calculate SRP
#' pw <- srp(qvalues, pi0, FDR = 0.05)
#' pw
#'

srp <- function(qvalues, pi0, FDR = 0.05){

  qvalues <- qvalues[!is.na(qvalues)]
  q <- sum(qvalues <= FDR)
  n_tests <- length(qvalues)

  pw <- ((1 - FDR) * q) / ((1 - pi0) * n_tests)

  if(pw > 1){
    warning("Infinite power", call. = FALSE)
    pw <- 1
    return(pw)
  }

  return(pw)
}


