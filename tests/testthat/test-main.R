testthat::test_that("ccp Function", {
  testthat::expect_error(
    ccp(test_data)
  )

  testthat::expect_error(
    ccp(test_data, washoutL_plus_trainL = 100, trainL = 100, washoutL = 50)
  )

  testthat::expect_error(
    ccp(1:1000, trainL = 100)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, tol = 0)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, nboots = 0)
  )

  testthat::expect_warning(
    ccp(test_data, trainL = 100, nboots = 5, plot.it = FALSE)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, plot.it = 0)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, MBBblockL = 0)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, MBBblockL = nrow(test_data) + 1)
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, kappa = 0, plot.it = FALSE)
  )

  testthat::expect_warning(
    ccp(test_data, trainL = 100, kappa = 0.2, plot.it = FALSE)
  )

  testthat::expect_silent(
    ccp(test_data, trainL = 100, plot.it = FALSE, MBBblock = 20)
  )

  testthat::expect_silent(
    ccp(test_data, washoutL_plus_trainL = 150, plot.it = FALSE)
  )

  testthat::expect_silent(
    output <- ccp(test_data, trainL = 100, plot.it = FALSE)
  )

  testthat::expect_type(
    output,
    "list"
  )

  p <- plotCP(output)
  testthat::expect_true(
    ggplot2::is.ggplot(p)
  )

  testthat::expect_error(
    plotCP(1:1000)
  )

  testthat::expect_error(
    plotCP(output, nbreaks = 0)
  )

  testthat::expect_warning(
    plotCP(output, nbreaks = 50)
  )

  testthat::expect_true(
    output$netParams$error < 0.04
  )

  testthat::expect_true(
    min(output$MBBnull) >= 0 & max(output$MBBnull) <= 0.5
  )

  testthat::expect_true(
    output$statistic >= 0 & output$statistic <= 0.5
  )
}
)
