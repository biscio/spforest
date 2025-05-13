# Test for testcovtreee ----

## Part I ----
vecval0 <- lapply(spatstat.data::bei.extra, FUN = function(i) {
  c(as.matrix.im(i))
})

expect_silent(
  arbre <- tesscovtree(
    X = spatstat.data::bei,
    vecval = vecval0,
    areapixel = beisoilres[[1]]$xstep * beisoilres[[1]]$ystep,
    dimcov = beisoilres[[1]]$dim,
    covrangex = beisoilres[[1]]$xrange,
    covrangey = beisoilres[[1]]$yrange,
    listcovariates = bei.extra,
    mtry = 1,
    randmtry = TRUE,
    minpts = 1000
  )
)

### Test the object is correct

# 2 is length(listcovariates) in the example above
expect_length(arbre$namecov, 2)

expect_length(arbre$tree, 9)

expect_true(spatstat.geom::is.im(arbre$im))

expect_true(spatstat.geom::is.ppp(arbre$X))

expect_true(sum(rand_covar(
  listcovariates = list(1, 1, 1, 1, 1),
  mtry = 1 / 2,
  randmtry = TRUE
)) > 0)

expect_equal(arbre$tree[[8]]$nodeID, 8)

expect_equal(arbre$tree[[3]]$left_daughter, 6)

expect_equal(arbre$tree[[2]]$right_daughter, 5)

expect_equal(arbre$tree[[2]]$split_var, 2)

expect_equal(arbre$tree[[3]]$split_val, 141.62)

expect_equal(arbre$tree[[8]]$status, 0)

expect_equal(arbre$tree[[8]]$intensity_pred, 0.009661151,
  tolerance = 1e-5
)

### Test methods on sptree
expect_silent(print(arbre))
expect_silent(summary(arbre))
expect_equal(dim(summary(arbre))[2], 7)
# expect_equal(dim(summary(arbre, fulltree = T))[2], 11)

expect_equal(class(plot(arbre)), "im")

expect_length(predict(arbre), spatstat.geom::npoints(arbre$X))
expect_true(is.numeric(predict(arbre)))
expect_length(predict(
  object = arbre,
  newdata = spatstat.random::runifpoint(
    n = 50,
    win = arbre$X$window
  )
), 50)
expect_true(is.numeric(
  predict(
    object = arbre,
    newdata = spatstat.random::runifpoint(
      n = 50,
      win = arbre$X$window
    )
  )
))
expect_null(predict(
  object = arbre,
  newdata = spatstat.random::runifpoint(n = 0)
))

expect_true(is.numeric(predict(arbre, newdata = c(827, 319))))


expect_true(is.numeric(
  predicttree(
    object = arbre,
    newdata = spatstat.random::runifpoint(
      n = 50,
      win = arbre$X$window
    )
  )
))
expect_null(predicttree(arbre,
  newdata = spatstat.random::runifpoint(n = 0)
))
expect_equal(
  predict(arbre, newdata = c(827, 319)),
  predicttree(arbre, newdata = c(827, 319))
)


## Part II ----

areapixel0 <- beisoilres[[1]]$xstep * beisoilres[[1]]$ystep
vecval0 <- lapply(beisoilres, FUN = function(i) {
  c(spatstat.geom::as.matrix.im(i))
})
dimcov0 <- beisoilres[[1]]$dim
covrangex0 <- beisoilres[[1]]$xrange
covrangey0 <- beisoilres[[1]]$yrange

expect_silent(arbre <- tesscovtree(
  X = spatstat.data::bei,
  vecval = vecval0,
  areapixel = areapixel0,
  dimcov = dimcov0,
  covrangex = covrangex0,
  covrangey = covrangey0,
  listcovariates = beisoilres,
  mtry = 1,
  randmtry = TRUE,
  minpts = 500
))

expect_equal(
  names(arbre),
  c(
    "tree", "X", "namecov",
    "listcov", "im"
  )
)

expect_length(arbre$namecov, 15)


expect_length(arbre$tree, 23)

expect_true(spatstat.geom::is.im(arbre$im))

expect_true(spatstat.geom::is.ppp(arbre$X))


expect_equal(arbre$tree[[4]]$nodeID, 4)

expect_equal(arbre$tree[[4]]$left_daughter, 8)
expect_equal(arbre$tree[[4]]$right_daughter, 9)

expect_equal(arbre$tree[[4]]$split_var, 15)

expect_equal(arbre$tree[[4]]$split_val, 4.599173,
  tolerance = 1e-5
)

expect_equal(arbre$tree[[4]]$status, 1)

expect_equal(arbre$tree[[4]]$intensity_pred, 0.00862606,
  tolerance = 1e-5
)


A <- tesscovtree(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  vecval = vecval0,
  areapixel = areapixel0,
  dimcov = dimcov0,
  covrangex = covrangex0,
  covrangey = covrangey0,
  mtry = 1,
  minpts = 100
)
B <- format(object.size(A), units = "Mb")
expect_true(as.numeric(gsub(" Mb", "", B)) < 4)



# Test for tesstree ----

# TODO
