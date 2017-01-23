
#' Stable retrospective power
#'
#' \code{srp} uses p-values and false discovery rate (FDR) to calculate stable retrospective power.
#'
#' @param pvalues A vector of raw p-values.
#' @param FDR A level at which to control the false discovery rate (FDR). Deafults to 0.05.
#' @param ... Additional arguments passed to \code{\link{qvalue}} and \code{\link{pi0est}}.
#' @return A numeric that is the estimated stable retrospective power.
#' @seealso \code{\link{qvalue}} and \code{\link{pi0est}} how q-values and pi0 are calculated and for arguments.
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

srp <- function (pvalues, FDR = 0.05, ...)
{
  qobj <- try(qvalue::qvalue(pvalues, FDR, ...), silent = T)

  if(inherits(qobj, "try-error")){
    msg <- gsub("\\n","", qobj[1])
    stop(msg)
  }

  pi0 <- qobj$pi0

  if(pi0==1){
    stop("The estimated pi0 == 1: no effects. Power calculation is not meaningful.")
  }

  ntests <- length(pvalues)
  qsig <- sum(qobj$significant)

  th1 <- (1 - pi0) * ntests
  dh1 <- (1 - FDR) * qsig
  SRP <- dh1/th1

  data.frame(SRP, pi0)
}


