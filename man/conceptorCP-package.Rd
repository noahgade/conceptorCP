\name{conceptorCP-package}
\alias{conceptorCP-package}
\alias{conceptorCP}
\docType{package}
\title{
\packageTitle{conceptorCP}
}
\description{
\packageDescription{conceptorCP}
}
\details{
Package to estimate the most likely change point location in a multivariate time series. Fits a series of conceptor matrices to a representative training window of data, and compares the evolution of the RNN reservoir states to the original computed conceptor spaces. Method assumes that the training window is at least wide-sense cyclostationary, or there is not a long run trend present. The training window should be representative in the sense that it captures a full range of dynamics of the system. Change points are identified from a Kolmogorov-Smirnov like statistic based on a univariate sequence of derived cosine similarity measures. Significance estimates are obtained from a moving block bootstrap of the orginal data. Plotting functions included to help with the visualization of the data and the methodology.
}
\author{
Noah D. Gade \href{ndg5e@virginia.edu}
}
\examples{
ccp_output <- ccp(test_data)
plotCP(ccp_output)
}
