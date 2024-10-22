#' Change Point Detection with Conceptors
#'
#' @description Performs the conceptor change point algorithm to determine the location and significance of the most likely change point in a dependent, multivariate time series. Requires specification of a training length, or a combined washout and training length.
#'
#' @details Provides an estimate of the most likely change point location in a multivariate time series. Fits a series of conceptor matrices to a representative training window of data, and compares the evolution of the RNN reservoir states to the original computed conceptor spaces. Method assumes that the training window is at least wide-sense cyclostationary, or there is not a long run trend present. The training window should be representative in the sense that it captures a full range of dynamics of the system. Change points are identified from a CUSUM like statistic based on a univariate sequence of derived cosine similarity measures. Quantile estimates are obtained from a moving block bootstrap of the original data.
#'
#' @param data A T\code{x}d data set with variables as columns.
#' @param washoutL_plus_trainL Number of time points used for both reservoir washout and conceptor training.
#' @param trainL Number of time points used for conceptor training.
#' @param washoutL Number of time points used for reservoir washout.
#' @param tol Error tolerance for conceptor fit to the data.
#' @param nboots Number of bootstraps to estimate statistic null distribution.
#' @param MBBblockL Block length for moving block bootstrap.
#' @param plot.it Indicator to toggle export of plot.
#' @param kappa Small constant that sets the change point domain (set as <= 0.01).
#'
#' @return List of output:
#' \describe{
#' \item{\code{estimate}}{Estimated change point location.}
#' \item{\code{statistic}}{Statistic from method.}
#' \item{\code{MBBquant}}{Estimated quantile of statistic from MBB.}
#' \item{\code{MBBnull}}{Simulated null distribution from MBB.}
#' \item{\code{MBBblockL}}{Block length used for the MBB.}
#' \item{\code{statSeries}}{Time-ordered series of statistics.}
#' \item{\code{angles}}{Time-ordered cosine similarities between reservoir states and conceptor space.}
#' \item{\code{CLimits}}{Lower and upper limits of time points for a consistent estimate of the change point.}
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
#' ccp(test_data, washoutL_plus_trainL = 150)
#' ccp(test_data, trainL = 100)
#' ccp(test_data, trainL = 100, tol = 0.08, nboots = 500)
#' }
ccp <- function(data, washoutL_plus_trainL = "", trainL = "", washoutL = "", tol = 0.04, nboots = 240, kappa = 0.01, MBBblockL = NULL, plot.it = TRUE) {

  if(check_integer(washoutL_plus_trainL) == FALSE & check_integer(trainL) == FALSE) {
    stop("Please specify a number of time points for conceptor training. Specify one of the following parameters: washoutL_plus_trainL, trainL")
  } else if(check_integer(washoutL_plus_trainL) == TRUE & check_integer(washoutL) == TRUE & check_integer(trainL) == FALSE) {
    trainL <- washoutL_plus_trainL - washoutL
  } else if(check_integer(washoutL_plus_trainL) == FALSE & check_integer(washoutL) == TRUE & check_integer(trainL) == TRUE) {
    washoutL_plus_trainL <- washoutL + trainL
  } else if(check_integer(washoutL_plus_trainL) == TRUE & check_integer(washoutL) == FALSE & check_integer(trainL) == TRUE) {
    washoutL <- washoutL_plus_trainL - trainL
  } else if(check_integer(washoutL_plus_trainL) == TRUE & check_integer(washoutL) == TRUE & check_integer(trainL) == TRUE) {
    if(washoutL_plus_trainL != (washoutL + trainL)) {
      stop("Please check specifications of washout and training lengths. Specify parameters such that washoutL_plus_trainL = washoutL + trainL")
    }
  }

  if(length(dim(data)) != 2) {
    stop("Please enter data in format Txd with columns as variables.")
  }

  if(tol < 1e-3 | tol > 0.50) {
    stop("Please specify a relevant error tolerance.")
  }

  if(nboots <= 0) {
    stop("Please specify a positive number of bootstrap iterations.")
  } else if(nboots <= 100) {
    warning("Consider increasing the number of bootstrap iterations.")
  } else if(nboots > 1000) {
    warning("Computation time may be hindered by a large number of bootstrap iterations.")
  }

  if(is.logical(plot.it) == FALSE) {
    stop("Please specify plot.it to include or suppress plot output.")
  }

  if(!is.null(MBBblockL)) {
    if(!check_integer(MBBblockL) | MBBblockL >= nrow(data)) {
      stop("Please specify a valid block length for the MBB.")
    }
  }

  if(kappa > 0.01) {
    warning("Consider changing kappa to <= 0.01 for consistent change point estimation.")
  } else if(kappa <= 0) {
    stop("kappa must be a small constant > 0.")
  }

  L <- nrow(data)
  CRNNFit <- fitCRNN(data, washoutL_plus_trainL, trainL, washoutL, tol)
  KSseries <- KSstatCalc(CRNNFit$output$angles[(CRNNFit$params$washoutL + CRNNFit$params$trainL + 1):L], kappa)
  statistic <- max(KSseries)
  estimate <- which.max(KSseries) + CRNNFit$params$washoutL + CRNNFit$params$trainL

  if(is.null(MBBblockL)) {
    MBBblockL <- chooseMBBblockL(data, CRNNFit, kappa)
  }
  binput <- replicate(nboots, bootdata(data, CRNNFit$params$washoutL, CRNNFit$params$trainL, MBBblockL))
  MBBnull <- CRNNBootstrap(binput, CRNNFit$output$W, CRNNFit$output$C, 0, CRNNFit$params$washoutL,
                           CRNNFit$params$trainL, CRNNFit$params$bscale, CRNNFit$params$iscale, kappa)
  MBBquant <- sum(statistic <= MBBnull) / nboots

  CLimits <- CRNNFit$params$washoutL + CRNNFit$params$trainL + (L - CRNNFit$params$washoutL - CRNNFit$params$trainL)*(0.5 + c(-0.5*sqrt(1 - 4*kappa^2), 0.5*sqrt(1 - 4*kappa^2)))

  if(!dplyr::between(estimate, CLimits[1], CLimits[2])) {
    warning("Consider lowering kappa, CP estimate may not be a consistent estimate.")
  }

  output <- list("estimate" = estimate,
                 "statistic" = statistic,
                 "MBBquant" = MBBquant,
                 "MBBnull" = MBBnull,
                 "MBBblockL" = MBBblockL,
                 "statSeries" = KSseries,
                 "angles" = CRNNFit$output$angles,
                 "CLimits" = CLimits,
                 "netParams" = append(CRNNFit$params, list(C = CRNNFit$output$C, W = CRNNFit$output$W, ResStates = CRNNFit$output$ResStates)))

  if(plot.it == TRUE) {
    CPplot <- plotCP(output)
    print(CPplot)
  }

  return(output)
}
