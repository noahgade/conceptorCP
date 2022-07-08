#' Scale columns of data to range (-1,1)
#'
#' @description Scales columns of data such that the middle 95% of the values are between -1 and 1 without changing the relative center.
#'
#' @param input_data A data set with columns as variables.
#'
#' @return The scaled data set with the same dimension as \code{input_data}.
#' @importFrom dplyr %>%
scalein <- function(input_data) {
  out <- tidyr::as_tibble(input_data, .name_repair = "unique") %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), function(x) 2 * (x - stats::quantile(x, 0.025)) / (stats::quantile(x, 0.975) - stats::quantile(x, 0.025)) - 1))
  return(out)
}

#' Randomly initialize RNN matrices
#' @description Randomly initializes matrices Wbias, Winput, and Wres for use in a recurrent neural network. Wbias and Winput are dense matrices with random normal entries and Wres is a sparse matrix with a spectral radius of 0.8.
#'
#' @param N RNN reservoir size.
#' @param dimen Dimension of original data, d.
#'
#' @return A N\code{x}(N+d+1) matrix of stacked RNN matrices in order: \cr
#' \[Wbias (N\code{x}1), Winput (N\code{x}d), Wres (N\code{x}N)\].
initRNN <- function(N, dimen) {
  rscale <- 0.8
  W <- matrix(stats::rnorm(N*dimen + N), nrow = N)
  repeat{
    Wres_sparse <- Matrix::rsparsematrix(N, N, nnz = 10*N, rand.x = function(x) stats::rnorm(x))
    maxeigval <- max(abs(eigen(Wres_sparse, only.values = T)$values))
    Wres <- rscale / maxeigval * Wres_sparse
    if(sum(is.infinite(Wres)) == 0){break}
  }
  W <- cbind(W, as.matrix(Wres))
  return(W)
}

#' Fit RNNs and conceptor matrices to data
#' @description Fits 100 RNNs to the data and computes the corresponding conceptor matrices. RNN parameters, initialized matrices, reservoir states, conceptor matrices, and cosine similarity measures between the reservoir states and conceptor spaces are stored.
#'
#' @param data A T\code{x}d data set with columns as variables.
#' @param trainL Number of time points used for conceptor training.
#' @param washoutL Number of time points used for reservoir washout.
#' @param tol Error tolerance for conceptor fit to the data.
#'
#' @return List of RNN parameters and output from conceptor fit process.
fitCRNN <- function(data, trainL, washoutL = "", tol = 0.04) {
  L <- nrow(data)
  dimen <- ncol(data)
  input <- as.matrix(scalein(data))

  if(washoutL == ""){
    TEMPwashoutL <- 50
  }else{
    TEMPwashoutL <- washoutL
  }

  N <- 10 * dimen
  regular <- 1e-4
  bscales <- seq(0.1, 0.5, 0.2)
  iscales <- seq(0.2, 1.4, 0.4)
  ninit <- 10
  rscale <- 0.8

  W1 <- replicate(ninit, initRNN(N, dimen))
  trainErrors <- RNNParamFit(input, W1, 0, TEMPwashoutL, trainL, bscales, iscales)
  error1 <- min(trainErrors)

  bscale <- bscales[which(trainErrors == error1, arr.ind = T)[1]]
  iscale <- iscales[which(trainErrors == error1, arr.ind = T)[2]]

  error <- 1
  aperture <- N
  while(error > tol) {
    W2 <- replicate(ninit, initRNN(N, dimen))
    rnn0 <- runRNN(input, W2, 0, bscale, iscale)

    if(washoutL == "") {
      rnn1 <- runRNN(input, W2, 1, bscale, iscale)
      washoutL <- 10 * ceiling((min(which(apply(rnn1 - rnn0, 3, max) < 1e-6)) - 1) / 10)
    }

    C <- conceptorCalc(rnn0[,,(washoutL + 2):(washoutL + trainL + 1)], aperture)
    Wout <- WoutCalc(input[(washoutL + 1):(washoutL + trainL),], rnn0[,,(washoutL + 2):(washoutL + trainL + 1)], regular)
    crnn0 <- runCRNN(input, W2, C, 0, washoutL, bscale, iscale)
    output <- outputCalc(Wout, crnn0[,,(L + 1 + washoutL + 2):(L + 1 + washoutL + trainL + 1)])
    error <- mean(apply(output, 3, NRMSE, input[(washoutL + 1):(washoutL + trainL),]))

    if(error > tol) {
      if(aperture > 100 * N) {
        N <- N * dimen
        aperture <- N
      } else {
        aperture <- sqrt(10) * aperture
      }
    }
    if(N > 1000) {
      print('Unable to meet error tolerance - consider raising threshold.')
      break
    }
  }

  nfits <- 100
  W <- replicate(nfits, initRNN(N, dimen))
  rnn <- runRNN(input, W, 0, bscale, iscale)
  C <- conceptorCalc(rnn[,,(washoutL + 2):(washoutL + trainL + 1)], aperture)
  crnn <- runCRNN(input, W, C, 0, washoutL, bscale, iscale)

  angles <- apply(angleCalc(crnn[,,2:(L + 1)], crnn[,,(L + 3):(2*L + 2)]), 1, mean)

  RNNparams <- list("rscale" = rscale,
                    "iscale" = iscale,
                    "bscale" = bscale,
                    "N" = N,
                    "aperture" = aperture,
                    "washoutL" = washoutL,
                    "trainL" = trainL,
                    "error" = error)

  RNNoutput <- list("angles" = angles,
                    "C" = C,
                    "W" = W,
                    "ResStates" = crnn)

  out = list("params" = RNNparams, "output" = RNNoutput)
  return(out)
}

#' Generate bootstrapped data
#' @description Creates bootstrapped data for MBB procedure of the same dimension as the original data. Randomly appends data sections of a specified block length after a washout and training period.
#'
#' @param data A T\code{x}d data set with columns as variables.
#' @param washoutL Number of time points used for reservoir washout.
#' @param trainL Number of time points used for conceptor training.
#' @param blockL Number of time points for each block in bootstrap procedure.
#'
#' @return Bootstrapped data for MBB procedure of the same dimension as \code{data}.
bootdata <- function(data, washoutL, trainL, blockL) {
  L <- nrow(data)
  input <- as.matrix(scalein(data))
  wrapinput <- rbind(input, input[(washoutL + trainL + 1):L,])
  sindex <- sample((washoutL + trainL + 1):L, ceiling((L - washoutL - trainL) / blockL))
  bootinput <- as.matrix(rbind(as.matrix(input[1:(washoutL + trainL),]), do.call(rbind, lapply(sindex, function(x, dat, blockL) as.matrix(dat[x:(x + blockL - 1),]), wrapinput, blockL)))[1:L,])
  return(bootinput)
}

#' Adjust ggplot legend
#'
#' @param data Data series to plot.
#' @param params Parameters of plot.
#' @param size Size of key.
#'
#' @return Adjusted plot settings.
draw_key_cust <- function(data, params, size) {
  if (data$colour == "blue") {
    ggplot2::draw_key_vpath(data, params, size)
  } else {
    ggplot2::draw_key_path(data, params, size)
  }
}
