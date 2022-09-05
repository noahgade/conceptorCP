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
  nconnect <- ifelse(N < 10, N, 10)
  repeat{
    Wres_sparse <- Matrix::rsparsematrix(N, N, nnz = nconnect*N, rand.x = function(x) stats::rnorm(x))
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
#' @param washoutL_plus_trainL Number of time points used for both reservoir washout and conceptor training.
#' @param trainL Number of time points used for conceptor training.
#' @param washoutL Number of time points used for reservoir washout.
#' @param tol Error tolerance for conceptor fit to the data.
#'
#' @return List of RNN parameters and output from conceptor fit process.
fitCRNN <- function(data, washoutL_plus_trainL, trainL, washoutL, tol) {
  L <- nrow(data)
  dimen <- ncol(data)
  input <- as.matrix(scalein(data))

  if(check_integer(washoutL) == FALSE) {
    TEMPwashoutL <- 50
    if(check_integer(trainL) == FALSE) {
      trainL <- washoutL_plus_trainL - TEMPwashoutL
      scenario <- 1
    } else {
      washoutL_plus_trainL <- trainL + TEMPwashoutL
      scenario <- 2
    }
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

    if(check_integer(washoutL) == FALSE) {
      rnn1 <- runRNN(input, W2, 1, bscale, iscale)
      washoutL <- 10 * ceiling((min(which(apply(rnn1 - rnn0, 3, max) < 1e-6)) - 1) / 10)

      if(scenario == 1) {
        trainL <- washoutL_plus_trainL - washoutL
      } else if (scenario == 2) {
        washoutL_plus_trainL <- trainL + washoutL
      }
    }

    C <- conceptorCalc(rnn0[,,(washoutL + 2):(washoutL + trainL + 1)], aperture)
    Wout <- WoutCalc(as.matrix(input[(washoutL + 1):(washoutL + trainL),]), rnn0[,,(washoutL + 2):(washoutL + trainL + 1)], regular)
    crnn0 <- runCRNN(input, W2, C, 0, washoutL, bscale, iscale)
    output <- outputCalc(Wout, crnn0[,,(L + 1 + washoutL + 2):(L + 1 + washoutL + trainL + 1)])
    error <- mean(apply(output, 3, NRMSE, as.matrix(input[(washoutL + 1):(washoutL + trainL),])))

    if(error > tol) {
      if(aperture > 100 * N) {
        N <- N * dimen
        aperture <- N
      } else {
        aperture <- sqrt(10) * aperture
      }
    }
    if(N > 1000) {
      warning('Unable to meet error tolerance - consider raising threshold.')
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
  wrapinput <- rbind(input, as.matrix(input[(washoutL + trainL + 1):L,]))
  sindex <- sample((washoutL + trainL + 1):L, ceiling((L - washoutL - trainL) / blockL))
  bootinput <- as.matrix(rbind(as.matrix(input[1:(washoutL + trainL),]), do.call(rbind, lapply(sindex, function(x, dat, blockL) as.matrix(dat[x:(x + blockL - 1),]), wrapinput, blockL)))[1:L,])
  return(bootinput)
}

#' Select optimal MBB block length
#' @description Selects optimal MBB block length with the Hall, Horowitz, and Jing (1995) algorithm.
#'
#' @param data A T\code{x}d data set with columns as variables.
#' @param CRNNFit Output from fitCRNN function.
#'
#' @return Estimated optimal block length for MBB per Hall, Horowitz, and Jing (1995).
chooseMBBblockL <- function(data, CRNNFit) {
  Lstar <- CRNNFit$params$washoutL

  M <- 40
  B <- 5
  L <- nrow(data)
  Lb <- ceiling(seq(L^(1/5), L^(1/2), length.out = B))
  split_input <- array(dim = c(L - M + 1, ncol(data), M*B))
  for(b in 1:B){
    for(m in 1:M){
      split_input[,,(b-1)*M + m] <- bootdata(rbind(data[1:(CRNNFit$params$washoutL + CRNNFit$params$trainL),],
                                              data[(CRNNFit$params$washoutL + CRNNFit$params$trainL + m):(L - M + m),]),
                                            CRNNFit$params$washoutL, CRNNFit$params$trainL, Lb[b])
    }
  }

  base_input <- replicate(M, bootdata(data, CRNNFit$params$washoutL, CRNNFit$params$trainL, Lstar))

  MBBests <- CRNNBootstrap(split_input, CRNNFit$output$W, CRNNFit$output$C, 0, CRNNFit$params$washoutL,
                             CRNNFit$params$trainL, CRNNFit$params$bscale, CRNNFit$params$iscale)

  MBBbase <- CRNNBootstrap(base_input, CRNNFit$output$W, CRNNFit$output$C, 0, CRNNFit$params$washoutL,
                             CRNNFit$params$trainL, CRNNFit$params$bscale, CRNNFit$params$iscale)

  Lstar <- ceiling((L / (L - M + 1))^(1/3) * Lb[which.min(colSums(matrix((MBBests - mean(MBBbase))^2, nrow = M)))])
  return(Lstar)
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

#' Check if value is a positive integer
#'
#' @param value Numeric value to check if positive integer.
#'
#' @return Boolean value (T/F) if value is a positive integer.
check_integer <- function(value) {
  if(is.numeric(value) == FALSE) {
    output <- FALSE
  } else if(value <= 0) {
    output <- FALSE
  } else {
    output <- (abs(value - round(value)) < .Machine$double.eps^0.5)
  }
  return(output)
}

#' Vizualize conceptor CP method
#' @description Plots estimate and internal dynamics of the conceptor change point method.
#'
#' @details Plots the time-ordered series of Kolmogorov-Smiornov like statistics from the conceptor change point method along with quantiles of the moving block bootstrap null distribution. Provides a visual aid of the relationship between the computed conceptor spaces and the propagating reservoir states over time with a cosine similarity measure. A comparison of the empirical disribution functions for windows of the cosine similarities are also included. Shading is relative and not on the same scale for all plots. Red shading represents time points where the dynamics are further away from the original conceptor space, and blue shading represents dynamics closer to the training window and the conceptor space.
#'
#' @param ccp_output Output from the conceptorCP function.
#' @param nbreaks Number of windows to divide series for visual.
#'
#' @return Plots in the following order:
#' \itemize{
#' \item \emph{Top}: Time-ordered series of statistics with most likely change point location and MBB results.
#' \item \emph{Middle}: Time-ordered cosine similarities between reservoir states and conceptor space.
#' \item \emph{Bottom}: Comparison of cosine similarity empirical CDFs (each window against the full time series).
#' }
#' @export
#' @importFrom dplyr %>%
#' @importFrom grDevices rgb
#' @importFrom scales label_number
#'
#' @examples
#' \donttest{
#' ccp_output <- ccp(test_data, washoutL_plus_trainL = 150)
#' plotCP(ccp_output)
#' }
plotCP <- function(ccp_output, nbreaks = 10) {
  if(check_integer(nbreaks) == FALSE) {
    stop("Enter an integer number of break points for plot.")
  }

  if(nbreaks < 2 | nbreaks > 40) {
    warning("Plot may be improved by adjusting number of break points.")
  }

  Time <- Window <- Values <- Reference <- WindowLength <- RCDF <- NULL
  Angles <- dplyr::tibble(Angles = ccp_output$angles)
  L <- nrow(Angles)
  Angles <- Angles %>% dplyr::mutate(Time = seq(1, L))
  PWAngles <- Angles[(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1):L,]
  AngleMin <- min(PWAngles) - (1 - min(PWAngles)) / 100
  PWAngles <- PWAngles %>% dplyr::mutate(PWRanks = rank(PWAngles$Angles) / nrow(PWAngles))
  EndPts <- floor(seq(0, nrow(PWAngles), length.out = nbreaks + 1)) + ccp_output$netParams$washoutL + ccp_output$netParams$trainL
  PWAngles <- PWAngles %>% dplyr::rowwise() %>% dplyr::mutate(Window = sum(EndPts < Time))
  PWAngles$Values <- paste("Value", unlist(sapply(diff(EndPts), seq, simplify = F)))

  plotM <- ggplot2::ggplot(PWAngles, ggplot2::aes_string(x = "Time", y = "Angles")) +
    ggplot2::geom_segment(data = PWAngles, ggplot2::aes_string(x = "Time", xend = "Time", y = "AngleMin", yend = 1, color = "PWRanks")) +
    ggplot2::scale_color_gradientn(name = "",
                                   colors = c(rgb(1, 0, 0), rgb(1, 1, 1, 0), rgb(0, 0, 1)),
                                   limits = c(0, 1), na.value = "white",
                                   breaks = c(0, 0.5, 1),
                                   labels = c("Away from \nConceptor \nSpace", "\nMiddle: Percentiles of Cosine Similarities \n\nBottom: Relative ECDF Difference", "Towards \nConceptor \nSpace")) +
    ggplot2::geom_point(ggplot2::aes_string(x = "Time", y = "Angles"), shape = 20, size = 1) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(name = "Time",
                                expand = c(0, 0),
                                breaks = EndPts + c(rep(1, nbreaks), 0)) +
    ggplot2::scale_y_continuous(name = paste("Cosine Similaritiy\n Smaller  \u2194  Larger", sep = ''),
                                labels = label_number(accuracy = 0.1),
                                limits = c(AngleMin, 1),
                                breaks = 1,
                                expand = c(0, 0)) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(color = "white", size = 2),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "bottom",
                   legend.text = ggplot2::element_text(size = 10),
                   axis.title.x = ggplot2::element_text(size = 10),
                   axis.text.x = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 10),
                   legend.justification = "top",
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   legend.key.height = ggplot2::unit(0.3, "cm"),
                   legend.key.width = ggplot2::unit(2, "cm"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  CDF <- dplyr::tibble(Reference = rep(sort(PWAngles$Angles), nbreaks), RCDF = rep(seq(nrow(PWAngles)) / nrow(PWAngles), nbreaks), Window = rep(seq(1, nbreaks), each = nrow(PWAngles)), WindowLength = rep(diff(EndPts), each = nrow(PWAngles)))
  WCDF <- dplyr::select(PWAngles, Angles, Window, Values) %>% tidyr::pivot_wider(names_from = Values, values_from = Angles)
  WCDF <- dplyr::left_join(CDF, WCDF, by = "Window")
  WCDF <- dplyr::transmute(WCDF, dplyr::across(5:max(diff(EndPts)), function(X) X <= Reference))
  CDF$WCDF <- dplyr::rowwise(WCDF) %>% dplyr::transmute(WCDF = sum(dplyr::c_across(cols = dplyr::everything())))
  CDF <- dplyr::mutate(CDF, WCDF = WCDF$WCDF / WindowLength, Shading = RCDF - WCDF)

  plotB <- ggplot2::ggplot(data = CDF, ggplot2::aes_string(x = "RCDF", y = "WCDF")) +
    ggplot2::geom_segment(ggplot2::aes_string(x = "RCDF", xend = "RCDF", y = 0, yend = 1, color = "Shading")) +
    ggplot2::scale_color_gradientn(name = "",
                                   colors = c(rgb(1, 0, 0), rgb(1, 1, 1, 0.5), rgb(0, 0, 1)),
                                   limits = c(-1, 1), na.value = "white") +
    ggplot2::facet_wrap(. ~ Window, nrow = 1) +
    ggplot2::geom_point(data = CDF, ggplot2::aes_string(x = "RCDF", y = "WCDF"), col = "black", shape = 20, size = 0.2) +
    ggplot2::scale_x_continuous(name = "ECDF (Full Time Series)", limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = "ECDF (Time Window)", limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::theme(strip.text.x = ggplot2::element_blank(),
                   strip.text.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(color = "black", fill = "white"),
                   panel.spacing.y = ggplot2::unit(0.05, "cm"),
                   panel.spacing.x = ggplot2::unit(0, "cm"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(0.3, "cm"),
                   legend.key.height = ggplot2::unit(0.5, "cm"),
                   axis.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 10),
                   legend.justification = "bottom",
                   axis.text.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = 10),
                   legend.position = "none",
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  SSeries <- dplyr::tibble(Time = seq(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1, L), Stat = ccp_output$statSeries)
  upper.limit <- max(SSeries$Stat, stats::quantile(ccp_output$MBBnull, 0.99)[[1]]) + 0.25

  plotT <- ggplot2::ggplot(SSeries) +
    ggplot2::geom_line(ggplot2::aes_string(x = "Time", y = "Stat")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.90)[[1]], linetype = "Upper 10% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.95)[[1]], linetype = "Upper 5% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(ccp_output$MBBnull, 0.99)[[1]], linetype = "Upper 1% \nMBB Null Dist."), color = "red", key_glyph = "cust") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = ccp_output$estimate, linetype = "Most Likely \nChange Point"), color = "blue", key_glyph = "cust") +
    ggplot2::scale_linetype_manual("", values = c("Most Likely \nChange Point" = "solid", "Upper 10% \nMBB Null Dist." = "dotted", "Upper 5% \nMBB Null Dist." = "dashed", "Upper 1% \nMBB Null Dist." = "longdash"),
                                   guide = ggplot2::guide_legend(override.aes = list(colour = c("blue", "red", "red", "red")))) +
    ggplot2::scale_x_continuous(name = "Time", limits = c(ccp_output$netParams$washoutL + ccp_output$netParams$trainL + 1, L),
                                expand = c(0, 0), EndPts + c(rep(1, nbreaks), 0)) +
    ggplot2::scale_y_continuous(name = "Statistic", limits = c(0, upper.limit), labels = label_number(accuracy = 0.1),
                                expand = c(0, 0)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10),
                   legend.position = "top",
                   legend.justification = "right",
                   axis.title.x = ggplot2::element_text(size = 10),
                   axis.text.x = ggplot2::element_text(size = 10),
                   axis.title.y = ggplot2::element_text(size = 10),
                   plot.margin = ggplot2::unit(c(0.2, 0.4, 0.2, 0.2), "cm"),
                   plot.caption = ggplot2::element_text(size = 12),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  plotS <- cowplot::axis_canvas(plotT, axis = "y", coord_flip = TRUE) +
    ggplot2::geom_histogram(ggplot2::aes(x = ccp_output$MBBnull), binwidth = 0.4, color = "black", fill = rgb(1, 0, 0, 0.2)) +
    ggplot2::coord_flip()
  suppressWarnings({
    plot1 <- cowplot::insert_yaxis_grob(plotT, plotS, grid::unit(0.1, "null"), position = "right")
    plot2 <- cowplot::plot_grid(plotM + ggplot2::theme(legend.position = "none"), plotB, ncol = 1, align = "hv", axis = "lr", rel_heights = c(1.1, 1))
    plot3 <- cowplot::insert_yaxis_grob(plot2, grid::nullGrob(), grid::unit(0.1, "null"), position = "right")
    plot4 <- cowplot::plot_grid(plot1, plot3, ncol = 1)
    plot5 <- cowplot::plot_grid(plot4, cowplot::get_legend(plotM), rel_heights = c(1, 0.15), ncol = 1, align = "hv", axis = "r")})
  return(plot5)
}

