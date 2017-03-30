
#' Stable retrospective power
#'
#' \code{srp} uses p-values and false discovery rate (FDR) to calculate stable retrospective power.
#'
#' @param pvalues A vector of raw p-values.
#' @param FDR A level at which to control the false discovery rate (FDR). Deafults to 0.05.
#' @param ... Additional arguments passed to \code{\link{qvalue}} and \code{\link{pi0est}}.
#' @return A data.frame containing:
#'  SRP that is the estimated stable retrospective power.
#'  pi0 propotion of the true non-null effects.
#'  fp the number of the statistically significant effects that are likely to be false positives.
#'  rs the number of statistically significant true effects seen in this study that will most likely come up as significant in a replication study.
#'  ud the number undiscovered results awaiting to be uncovered by replication studies.
#'
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

  k <- length(pvalues)
  d <- sum(qobj$significant)

  tnn <- (1 - pi0) * k
  td <- (1 - FDR) * d

  if(td>tnn){
    stop("Power calculation is not reliable.")
  }

  SRP <- td/tnn
  fp <- FDR * d
  rs <- td * SRP
  ud <- tnn - td

  data.frame(SRP, pi0, fp, rs, ud)
}


