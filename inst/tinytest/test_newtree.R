expect_silent(intensitytree(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  mtry = 1,
  minpts = 500
))

