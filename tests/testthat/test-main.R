testthat::test_that("InputErrors", {
  testthat::expect_error(
    ccp(test_data),
    "Error in `ccp(test_data)`: Please specify a number of time points for conceptor training. Specify one of the following parameters: washoutL_plus_trainL, trainL"
  )

  testthat::expect_error(
    ccp(test_data, washoutL_plus_trainL = 100, trainL = 100, washoutL = 50),
    "Please check specifications of washout and training lengths. Specify parameters such that washoutL_plus_trainL = washoutL + trainL"
  )

  testthat::expect_error(
    ccp(1:1000, trainL = 100),
    "Please format data with columns as variables."
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, tol = 0),
    "Please specify a relevant error tolerance."
  )

  testthat::expect_error(
    ccp(test_data, trainL = 100, nboots = 0),
    "Please specify a positive number of bootstrap iterations."
  )

  testthat::expect_warning(
    ccp(test_data, trainL = 100, nboots = 5),
    "Consider increasing the number of bootstrap iterations."
  )

  test_that::expect_error(
    ccp(test_data, trainL = 100, plot.it = 0),
    "Please specify plot.it to include or suppress plot output."
  )
}
)
