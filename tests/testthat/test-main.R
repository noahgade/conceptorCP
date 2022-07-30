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

  testthat::expect_type(
    ccp(test_data, trainL = 100, plot.it = FALSE),
    "list"
  )

  testthat::expect_silent(
    ccp(test_data, washoutL_plus_trainL = 150, plot.it = FALSE)
  )

  testthat::expect_output(
    ccp(test_data, trainL = 100)
  )
}
)
