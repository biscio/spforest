# Test spforest call tesscovforest ----

expect_silent(
  forest <- spforest(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 2,
    threshold = 0.01,
    cores = 1,
    mtry = 1 / 3,
    minpts = 100
  )
)

expect_equal(class(forest), "spforest")


# Test spforest call tessforest ----

expect_silent(
  forest2 <- spforest(
    X = spatstat.data::bei,
    Ntree = 5,
    lambda = 50,
    dimyx = c(51, 101),
    test.connected = FALSE,
    cores = 1
  )
)

expect_equal(class(forest2), "spforest")

expect_silent(
  extractforest(forest = forest2, whichtrees = c(1, 3))
)

# Test method for spforest from tessforest ----

expect_silent(print(forest2))
expect_silent(plot(forest))

# Test method for spforest from tesscovforest ----

expect_silent(print(forest))

expect_length(predict(forest), spatstat.geom::npoints(forest$X))
expect_true(is.numeric(predict(forest)))
expect_length(predict(
  object = forest,
  newdata = spatstat.random::runifpoint(
    n = 50,
    win = forest$X$window
  )
), 50)
expect_true(is.numeric(
  predict(
    object = forest,
    newdata = spatstat.random::runifpoint(
      n = 50,
      win = forest$X$window
    )
  )
))
expect_null(predict(
  object = forest,
  newdata = spatstat.random::runifpoint(n = 0)
))
expect_true(is.numeric(predict(forest,
  newdata = c(827, 319)
)))
expect_true(is.numeric(predict(forest,
  newdata = cbind(
    c(827, 319),
    c(418, 88)
  )
)))

# Testing boxplot and vipplot methods for spforest ----

expect_silent(boxplot(forest))


expect_silent(vipplot(forest, sorted = TRUE))
expect_silent(vipplot(forest, sorted = FALSE))
expect_true(vipplot(forest, sorted = TRUE)$Importance[1] >
  vipplot(forest, sorted = TRUE)$Importance[2])

