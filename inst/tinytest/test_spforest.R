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
