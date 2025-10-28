---
output: github_document
editor_options: 
  markdown: 
    wrap: 72
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# spforest

<!-- badges: start -->

<!-- badges: end -->

The goal of \texttt{spforest} is to estimate the intensity of a 2D point
process with spatial random intensity forest. The package but can be
considered in its pre-alpha version and is currently under heavy
development. In particular, its documentation can 
be considered in its infancy. 

## Installation

The package \texttt{spforest} depends on the \texttt{R} package
[spatstat.](https://github.com/spatstat/spatstat) 
which you should install first.

Then, you  can install the development version of \texttt{spforest} from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("biscio/spforest")
```

## Basic example without covariate

By default, the function \texttt{spforest} will work without 
any covariates. 

```{r}
library(spforest)
forestnocov <- spforest(
  X = spatstat.data::bei,
  Ntree = 100,
  listcovariates = NULL,
  lambda = 100,
  dimyx = c(50, 50),
  test.connected = FALSE
)
```
```{r}
plot(forestnocov)
```


## Basic example with covariates


```{r example}
library(spforest)
forest <- spforest(
  X = spatstat.data::bei,
  listcovariates = spatstat.data::bei.extra, ,
  Ntree = 50,
  mtry = 1,
  minpts = 200,
  cores = 1
)
```

```{r}
plot(forest)
```

