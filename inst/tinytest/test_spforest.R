# Test spforest call tesscovforest ----

expect_silent(
  forest <- spforest(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 2,
    threshold = 0.01,
    mtry = 1 / 3,
    randmtry = TRUE,
    minpts = 100
  )
)

expect_equal(class(forest), "spforest")


## test parallel future 
future::plan("multisession", workers = 2)
expect_silent(
  forest <- spforest(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    score = "lcv",
    p = 1 / 2,
    Ntree = 4,
    threshold = 0.01,
    mtry = 8,
    randmtry = FALSE,
    minpts = 100,
    parallel = TRUE
  )
)
future::plan("sequential")


# Test spforest with tesscovforest - other parameters ----

## mtry to choose a fixed number of param each times ----

expect_silent(
  forest <- spforest(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    score = "lcv",
    p = 1 / 2,
    Ntree = 4,
    threshold = 0.01,
    mtry = 8,
    randmtry = FALSE,
    minpts = 100
  )
)


# Test spforest call tessforest ----

expect_silent(
  forest2 <- spforest(
    X = spatstat.data::bei,
    Ntree = 5,
    gamma = 50,
    dimyx = c(51, 101),
    test.connected = FALSE
  )
)

expect_equal(class(forest2), "spforest")

expect_silent(
  extractforest(forest = forest2, whichtrees = c(1, 3))
)

## Test parallelism future 

future::plan("multisession", workers = 2)
expect_silent(
  forest2 <- spforest(
    X = spatstat.data::bei,
    Ntree = 4,
    gamma = 50,
    dimyx = c(51, 101),
    test.connected = FALSE,
    parallel = TRUE
  )
)
future::plan("sequential")
expect_error(
  forest2 <- spforest(
    X = spatstat.data::bei,
    Ntree = 4,
    gamma = 50,
    dimyx = c(51, 101),
    test.connected = FALSE,
    parallel = TRUE
  )
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

# expect_silent(boxplot(forest))
# 
# 
# expect_silent(vipplot(forest, sorted = TRUE))
# expect_silent(vipplot(forest, sorted = FALSE))
# expect_true(vipplot(forest, sorted = TRUE)$Importance[1] >
#   vipplot(forest, sorted = TRUE)$Importance[2])

