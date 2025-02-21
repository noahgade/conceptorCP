% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CRNNBootstrap}
\alias{CRNNBootstrap}
\title{Compute MBB null distribution}
\usage{
CRNNBootstrap(
  bootinput,
  W,
  C,
  RSInitial,
  washoutL,
  trainL,
  bscale,
  iscale,
  kappa
)
}
\arguments{
\item{bootinput}{bootstrapped data input}

\item{W}{cube of RNN matrices}

\item{C}{cube of conceptor matrices}

\item{RSInitial}{initial reservoir state}

\item{washoutL}{length of reservoir washout}

\item{trainL}{length of conceptor training}

\item{bscale}{bias matrix scaling}

\item{iscale}{input matrix scaling}

\item{kappa}{small constant > 0}
}
\value{
bootstrap null distribution of CUSUM-like statistics
}
\description{
Compute MBB null distribution
}
