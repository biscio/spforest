Mnsin <- 1 + sin(spforest::beisoilres[[10]] / 50)
z500sin <- 500 * Mnsin / integral(Mnsin)
Xpoi <- spatstat.random::rpoispp(lambda = z500sin, nsim = 1)


foo1 <- spforest(
  X = Xpoi,
  listcovariates = spforest::beisoilres[sample(15)],
  Ntree = 200,
  mtry = 0.5,
  minpts = 100,
  randmtry = TRUE,
  p = 0,
  score = "lcv",
  threshold = 50,
  cores = 2
)

expect_true(max(abs(as.im(foo1) - z500sin)) < 0.002)

vipfoo1 <- sapply(
  X = seq_along(foo1$listcov),
  FUN = function(i) {
    mean(importance(foo1,
                    id_cov = i,
                    viptype = 4
    ))
  }
)
names(vipfoo1) <- names(foo1$listcov)

expect_equal(names(which.max(abs(vipfoo1))), "Mn")

