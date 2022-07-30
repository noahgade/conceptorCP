#' Change Point Detection with Conceptors
#'
#' @description Performs the conceptor change point algorithm to determine the location and significance of the most likely change point in a dependent, multivariate time series. Requires specification of a training length, or a combined washout and training length.
#'
#' @details Provides an estimate of the most likely change point location in a multivariate time series. Fits a series of conceptor matrices to a representative training window of data, and compares the evolution of the RNN reservoir states to the original computed conceptor spaces. Method assumes that the training window is at least wide-sense cyclostationary, or there is not a long run trend present. The training window should be representative in the sense that it captures a full range of dynamics of the system. Change points are identified from a Kolmogorov-Smirnov like statistic based on a univariate sequence of derived cosine similarity measures. Significance estimates are obtained from a moving block bootstrap of the orginal data.
#'
#' @param data A T\code{x}d data set with variables as columns.
#' @param washL_plus_trainL Number of time points used for both reservoir washout and conceptor training.
#' @param trainL Number of time points used for conceptor training.
#' @param washoutL Number of time points used for reservoir washout.
#' @param tol Error tolerance for conceptor fit to the data.
#' @param nboots Number of bootstraps to estimate statistic null distribution.
#'
#' @return List of output:
#' \describe{
#' \item{\code{estimate}}{Estimated change point location.}
#' \item{\code{statistic}}{Statistic from method.}
#' \item{\code{MBBsig}}{Significance estimate from MBB.}
#' \item{\code{MBBnull}}{Simulated null distribution from MBB.}
#' \item{\code{statSeries}}{Time-ordered series of statistics.}
#' \item{\code{angles}}{Time-ordered cosine similarities between reservoir states and conceptor space.}
#' \item{\code{netParams}}{List of RNN parameters:}
#' \itemize{
#' \item{\code{rscale}} {Scale of reservoir matrix.}
#' \item{\code{iscale}} {Scale of input matrix.}
#' \item{\code{bscale}} {Scale of bias matrix.}
#' \item{\code{N}} {RNN reservoir size.}
#' \item{\code{aperture} }{Parameter for conceptor computation.}
#' \item{\code{washoutL}} {Number of time points used for reservoir washout.}
#' \item{\code{trainL}} {Number of time points used for conceptor training.}
#' \item{\code{error}} {Error of conceptor fit to data (<= tolerance).}
#' \item{\code{C}} {Conceptor matrices from each RNN fit.}
#' \item{\code{W}} {RNN matrices for each fit.}
#' \item{\code{ResStates}} {Time-ordered reservoir states for each RNN fit.}
#' }
#' }
#' @export
#' @importFrom dplyr %>%
#'
#' @examples
#' \donttest{
#' ccp(test_data, washL_plus_trainL = 150)
#' ccp(test_data, trainL = 100)
#' ccp(test_data, trainL = 100, tol = 0.08, nboots = 500)
#' }
ccp <- function(data, washL_plus_trainL = "", trainL = "", washoutL = "", tol = 0.04, nboots = 200, plot.it = TRUE) {

  if(check_integer(washL_plus_trainL) == FALSE & check_integer(trainL) == FALSE) {
    stop(cat("Error: Please specify a number of time points for conceptor training.\nSpecify one of the following parameters: washL_plus_trainL, trainL"))
  } else if(check_integer(washL_plus_trainL) == TRUE & check_integer(washoutL) == TRUE & check_integer(trainL) == FALSE) {
    trainL <- washL_plus_trainL - washoutL
  } else if(check_integer(washL_plus_trainL) == FALSE & check_integer(washoutL) == TRUE & check_integer(trainL) == TRUE) {
    washL_plus_trainL <- washoutL + trainL
  } else if(check_integer(washL_plus_trainL) == TRUE & check_integer(washoutL) == FALSE & check_integer(trainL) == TRUE) {
    washoutL <- washL_plus_trainL - trainL
  } else if(check_integer(washL_plus_trainL) == TRUE & check_integer(washoutL) == TRUE & check_integer(trainL) == TRUE) {
    if(washL_plus_trainL != (washoutL + trainL)) {
      stop(cat("Error: Please check specifications of washout and training lengths.\nSpecify parameters such that washL_plus_trainL = washoutL + trainL"))
    }
  }

  L <- nrow(data)
  CRNNFit <- fitCRNN(data, washL_plus_trainL, trainL, washoutL, tol)
  KSseries <- KSstatCalc(CRNNFit$output$angles[(CRNNFit$params$washoutL + CRNNFit$params$trainL + 1):L])
  statistic <- max(KSseries)
  estimate <- which.max(KSseries) + CRNNFit$params$washoutL + CRNNFit$params$trainL

  MBBblockL <- ceiling(L^(1/3))
  binput <- replicate(nboots, bootdata(data, CRNNFit$params$washoutL, CRNNFit$params$trainL, MBBblockL))
  MBBnull <- CRNNBootstrap(binput, CRNNFit$output$W, CRNNFit$output$C, 0, CRNNFit$params$washoutL,
                           CRNNFit$params$trainL, CRNNFit$params$bscale, CRNNFit$params$iscale)
  MBBsig <- sum(statistic <= MBBnull) / nboots

  output <- list("estimate" = estimate,
                 "statistic" = statistic,
                 "MBBsig" = MBBsig,
                 "MBBnull" = MBBnull,
                 "statSeries" = KSseries,
                 "angles" = CRNNFit$output$angles,
                 "netParams" = append(CRNNFit$params, list(C = CRNNFit$output$C, W = CRNNFit$output$W, ResStates = CRNNFit$output$ResStates)))

  if(plot.it == TRUE) {
    CPplot <- plotCP(output)
    print(CPplot)
  }

  return(output)
}
