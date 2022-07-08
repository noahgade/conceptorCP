#' Change Point Detection with Conceptors
#'
#' @description Performs the conceptor change point algorithm to determine the location and significance of the most likely change point in a dependent, multivariate time series.
#'
#' @details Provides an estimate of the most likely change point location in a multivariate time series. Fits a series of conceptor matrices to a representative training window of data, and compares the evolution of the RNN reservoir states to the original computed conceptor spaces. Method assumes that the training window is at least wide-sense cyclostationary, or there is not a long run trend present. The training window should be representative in the sense that it captures a full range of dynamics of the system. Change points are identified from a Kolmogorov-Smirnov like statistic based on a univariate sequence of derived cosine similarity measures. Significance estimates are obtained from a moving block bootstrap of the orginal data.
#'
#' @param data A T\code{x}d data set with variables as columns.
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
#' ccp(test_data, trainL = 100, washoutL = 50)
ccp <- function(data, trainL = 100, washoutL = "", tol = 0.04, nboots = 200) {
  L <- nrow(data)
  CRNNFit <- fitCRNN(data, trainL, washoutL, tol)
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
  return(output)
}

#' Viaualize conceptor CP method
#' @description Plots estimate and internal dynamics of the conceptor change point method.
#'
#' @details Plots the time-ordered series of Kolmogorov-Smiornov like statistics from the conceptor change point method along with quantiles of the moving block bootstrap null distribution. Provides a visual aid of the relationship between the computed conceptor spaces and the propagating reservoir states over time with a cosine similarity measure. A comparison of the empirical disribution functions for windows of the cosine similarities are also included. Shading is relative and not on the same scale for all plots. Red shading represents time points where the dynamics are further away from the original conceptor space, and blue shading represents dynamics closer to the training window and the conceptor space.
#'
#' @param conceptorCPoutput Output from the conceptorCP function.
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
#' ccp_output <- ccp(test_data)
#' plotCP(ccp_output, nbreaks = 10)
plotCP <- function(ccp_output, nbreaks = 10) {
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
    ggplot2::geom_histogram(ggplot2::aes(x = conceptorCPoutput$MBBnull), binwidth = 0.5, color = "black", fill = rgb(1, 0, 0, 0.2)) +
    ggplot2::coord_flip()
  suppressWarnings({
  plot1 <- cowplot::insert_yaxis_grob(plotT, plotS, grid::unit(0.1, "null"), position = "right")
  plot2 <- cowplot::plot_grid(plotM + ggplot2::theme(legend.position = "none"), plotB, ncol = 1, align = "hv", axis = "lr", rel_heights = c(1.1, 1))
  plot3 <- cowplot::insert_yaxis_grob(plot2, grid::nullGrob(), grid::unit(0.1, "null"), position = "right")
  plot4 <- cowplot::plot_grid(plot1, plot3, ncol = 1)
  plot5 <- cowplot::plot_grid(plot4, cowplot::get_legend(plotM), rel_heights = c(1, 0.15), ncol = 1, align = "hv", axis = "r")})
  return(plot5)
}
