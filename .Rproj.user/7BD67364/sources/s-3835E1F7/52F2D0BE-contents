
<!-- README.md is generated from README.Rmd. Please edit that file -->
EdGwas
======

<!-- badges: start -->
<!-- badges: end -->
The goal of EdGwas is to help clustering outcome components (traits) that share some feature (genetic component) using polygenic risk scores (PRS).

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("abuchardt/EdGwas")
```

Example
-------

This is a basic example on simulated data:

``` r
library(EdGwas)
N <- 1000 #
q <- 10 #
p <- 5000 #
set.seed(1)
# Sample 1
x0 <- matrix(rbinom(n = N*p, size = 2, prob = 0.3), nrow=N, ncol=p)
B <- matrix(0, nrow = p, ncol = q)
B[1, 1:2] <- 2
y0 <- x0 %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
```

Compute polygenic scores and coefficients

``` r
psobj <- ps.edgwas(x0, y0)
ps <- psobj$PS
beta <- psobj$beta
```

Create new sample

``` r
x <- matrix(rbinom(n = N*p, size = 2, prob = 0.3), nrow=N, ncol=p)
y <- x %*% B + matrix(rnorm(N*q), nrow = N, ncol = q)
ps <- x %*% beta
```

Run 10-fold cross-validation for edgwas

``` r
pc <- cv.edgwas(ps, y)
```

Plot cross-validated error curve

``` r
plot(pc, 1)
```

![](README-plot1.png)

Plot estimated optimal adjacency matrix

``` r
plot(pc, 2)
```

![](README-plot2.png)
