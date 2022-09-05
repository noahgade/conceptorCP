# conceptorCP

## Overview

**conceptor CP** is a package of tools to identify a change point in time-ordered multivariate and dependent data. Estimates of significance are obtained with a moving block bootstrap, and the assumption is an initial section of data at least wide-sense stationary, or it does not exhibit some long run trend.

## Example Data

Package includes one example VAR1 data set, `test_data`. Description of its generation is included in the documentation files.

## References
@article{HHJ1995,
    author = {HALL, PETER and HOROWITZ, JOEL L. and JING, BING-YI},
    title = "{On blocking rules for the bootstrap with dependent data}",
    journal = {Biometrika},
    volume = {82},
    number = {3},
    pages = {561-574},
    year = {1995},
    month = {09},
    issn = {0006-3444},
    doi = {10.1093/biomet/82.3.561}
}

## Installation

Install via CRAN:
`install.packages("conceptorCP")`
`library(conceptorCP)`

