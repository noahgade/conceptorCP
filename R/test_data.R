#' @title Test data for use in examples
#' @description 3-Dimensional VAR1 data set with chagne point at T = 500
#' @format Tibble with 1000 rows and 3 variables (Series1, Series2, Series3)
#'
#' @details
#' \code{set.seed(138)} \cr
#' \code{CoeffMatrix1 <- matrix(c(0.6, -0.4, 0, -0.2, 0.8, -0.3, 0, 0.1, 0.3), nrow = 3)} \cr
#' \code{CoeffMatrix2 <- matrix(c(0.2, -0.7, -0.2, 0, 0.9, -0.1, 0.2, 0.5, 0.1), nrow = 3)} \cr
#' \code{RandInnovs <- matrix(rnorm(3000), ncol = 3)} \cr
#' \code{test_data <- dplyr::tibble(Series1 = rnorm(1), Series2 = rnorm(1), Series3 = rnorm(1))} \cr
#' \code{for(j in 1:500) {} \cr
#' \code{   test_data[j + 1,] <- as.matrix(test_data[j,]) %*% CoeffMatrix1 + RandInnovs[j,]} \cr
#' \code{}} \cr
#' \code{for(j in 501:1000) {} \cr
#' \code{   test_data[j + 1,] <- as.matrix(test_data[j,]) %*% CoeffMatrix2 + RandInnovs[j,]} \cr
#' \code{}} \cr
#' \code{test_data <- test_data[-1,]} \cr
"test_data"
