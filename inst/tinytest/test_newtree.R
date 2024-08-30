expect_silent(intensitytree(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  mtry = 1,
  minpts = 500
))

ftree <- function(){
  treerec(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    mtry = 2/3,
    minpts = 20
  )
}

gtree <- function(){
  intensitytree(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    mtry = 2/3,
    minpts = 20
  )
}

library(microbenchmark)
microbenchmark(ftree(), gtree(), times = 10)
