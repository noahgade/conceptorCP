#' @title Test data for use in examples
#' @description 3-Dimensional VAR1 data set with change point at T = 500
#' @format Tibble with 1000 rows and 3 variables (Series1, Series2, Series3)
#'
#' @details Details of data creation:
#' \itemize{
#' \item Random seed set to 138 with \code{set.seed(138)}
#' \item VAR1 coefficient matrix for time points 1 - 500, \cr
#' \eqn{[0.6, -0.4, 0 // -0.2, 0.8, -0.3 // 0, 0.1, 0.3]}
#' \item VAR1 coefficient matrix for time points 501 - 1000, \cr
#' \eqn{[0.2, -0.7, -0.2 // 0, 0.9, -0.1 // 0.2, 0.5, 0.1]}
#' \item Random innovations generated with \code{matrix(rnorm(3000), ncol = 3)}
#' \item Each of the three variables initialized with \code{rnorm(1)}
#' }
"test_data"
