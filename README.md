
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slr

<!-- badges: start -->
<!-- badges: end -->

The goal of slr is to apply simple linear regression

## Installation

You can install the development version of slr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jerrywangrs/slr")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(slr)
data("iris")
slr(iris$Sepal.Length,iris$Sepal.Width)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #> Beta 0:  3.418947 
    #> Beta 1:  -0.0618848 
    #> r square: 0.01382265 
    #> F-statistic:  2.074427 p-value:  0.1518983 
    #> t-statistic:  -1.440287 p-value:  0.1518983
