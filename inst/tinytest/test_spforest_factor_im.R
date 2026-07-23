expect_silent(
  convertedim <- imfcttoreal(
    Z = cut(beisoilres[[1]],
      quantile(beisoilres[[1]]),
      include.lowest = T
    ),
    X = spatstat.data::bei
  )
)

expect_true(
  is.im(convertedim)
)

expect_true(
  convertedim$type == "real"
)

expect_silent(
  facbei <- lapply(beisoilres, FUN = function(i) {
    W <- spatstat.geom::cut.im(i,
      breaks = quantile(i),
      include.lowest = T
    )

    levels(W) <- letters[1:4]
    return(W)
  })
)

expect_silent(
  foo <- spforest(
    X = bei,
    listcovariates = facbei,
    Ntree = 2,
    mtry = 10,
    minpts = 400,
    parallel = F
  )
)
