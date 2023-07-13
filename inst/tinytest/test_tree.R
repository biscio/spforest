expect_silent(
  arbre <- treerec(
    X = spatstat.data::bei,
    threshold = 1000,
    score = "lcv2",
    listcovariates = list(
      grad = spatstat.data::bei.extra$grad,
      elev = spatstat.data::bei.extra$elev
    ),
    mtry = 1,
    tol = Inf,
    minpts = 50
  )
)

# 2 is length(listcovariates) in the example above
expect_length(arbre$namecov, 2)

expect_length(arbre$tree, 207)

expect_true(spatstat.geom::is.im(arbre$im))

expect_true(spatstat.geom::is.ppp(arbre$X))

expect_true(sum(rand_covar(
  listcovariates = list(1, 1, 1, 1, 1),
  mtry = 1 / 2
)) > 0)

expect_equal(arbre$tree[[37]]$nodeID, 37)

expect_equal(arbre$tree[[37]]$left_daughter, 60)

expect_equal(arbre$tree[[37]]$right_daughter, 61)

expect_equal(arbre$tree[[37]]$split_var, 1)

expect_equal(arbre$tree[[37]]$split_val, 0.04586901)

expect_equal(arbre$tree[[37]]$status, 1)

expect_equal(arbre$tree[[37]]$intensity_pred, 0.007962085,
  tolerance = 1e-5
)

expect_equal(arbre$tree[[37]]$scr_parent, -859.8082,
  tolerance = 1e-5
)


expect_length(arbre$tree[[37]]$improvement, 5)
