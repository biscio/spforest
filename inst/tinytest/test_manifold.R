expect_silent(
  res <- pptomesh(X=spatstat.data::bei,,
                  elev= spatstat.data::bei.extra$elev)
)

expect_silent(
  forest <- spforest(X = res)
)

expect_silent(
  forest2 <- spforest(X = res)
)

expect_silent(
  log(forest)
)

expect_silent(
  A <- plot.spforestmesh(log(forest+exp(-8)), points=TRUE),
  rgl::close3d()
)

expect_silent(
  B <- plot.spforestmesh(log(forest+exp(-8)), points=FALSE),
  rgl::close3d()
)

expect_silent(
  dummypponmesh(n=100, mesh = res$mesh)
)
rgl::close3d()