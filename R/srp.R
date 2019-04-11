
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
#' @seealso \code{\link{limma}} and \code{\link{propTrueNull}} how pi0 is calculated and for arguments.
#' @import limma
#' @export
#' @examples
#' z <- rnorm(200)
#' z[1:40] <- z[1:40]+2
#' 160/200
#' p <- 2*pnorm(-abs(z))
#'
#' # Calculate SRP
#' pw <- srp(pvalues, FDR = 0.05)
#' pw
#'

srp <- function (pvalues, method = "lfdr", FDR = 0.05, ...)
{
  pi0 <- try(limma::propTrueNull(pvalues, method, ...), silent = TRUE)

  if (inherits(pi0, "try-error")) {
    msg <- gsub("\\n","", pi0[1])
    stop(msg)
  }

  if (dplyr::near(pi0, 1)) {
    stop("The estimated pi0 == 1: no effects. Power calculation is not meaningful.")
  }

  k <- length(pvalues)
  d <- sum(pvalues < FDR)

  tnn <- (1 - pi0) * k
  td <- (1 - FDR) * d

  if (td > tnn) {
    stop("Power calculation is not reliable.")
  }

  SRP <- td / tnn
  fp <- FDR * d
  rs <- td * SRP
  ud <- tnn - td

  data.frame(SRP, pi0, fp, rs, ud)
}

