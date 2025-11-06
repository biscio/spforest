# Test tesscovforest ---- 

expect_silent(
  forest <- tesscovforest(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 3,
    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
    mtry = 1 / 3,
    randmtry = TRUE,
    minpts = 50,
    parallel = FALSE
  )
)
expect_length(forest$trees, 3)

expect_silent(
  forest2 <- tesscovforest(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 2,
    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
    mtry = 1 / 3,
    randmtry = TRUE,
    minpts = 50,
    parallel = FALSE
  )
)
expect_length(forest2$trees, 2)

expect_silent(merge(x=forest, y=forest2))
expect_length(merge(x=forest, y=forest2)$trees, 5)

expect_equal(forest$X, spatstat.data::bei)

expect_equal(sum(forest$listcov[[1]] - spatstat.data::bei.extra$elev), 0)
expect_equal(sum(forest$listcov[[2]] - spatstat.data::bei.extra$grad), 0)

expect_equal(forest$p, 1 / 2)

expect_equal(forest$mtry, 1 / 3)

expect_silent(
  imtrees <- lapply(forest$trees,
                    FUN = function(i) {
                      i$im
                    }
  )
)

expect_silent(
  imforest <- Reduce("+", imtrees) / length(forest$trees) / forest$p
)

expect_true(spatstat.geom::is.im(imforest))


expect_silent(
  outputoob <- OOBscr(forest=forest)
)

expect_inherits(current=OOBscr(forest=forest), class="numeric")


A <- tesscovforest(
  X = spatstat.data::bei,
  listcovariates = lapply(beisoilres, FUN=function(i){
    as.im(i, dimyx=c(10,20))
  }),
  Ntree = 10,
  minpts = 100,
  mtry = 1,
  p = 0,
  parallel = FALSE
)
B <- format(object.size(A), units = "Mb")
expect_true(as.numeric(gsub(" Mb", "",B)) < 2.2)


expect_inherits(
  OOBoptim(X = spatstat.data::bei, 
           listcovariates = spatstat.data::bei.extra, 
           params = list(Ntree = 5,
                         mtry = 2, 
                         minpts = c(500,800)),
           parallel = FALSE),
  class="data.frame"
)

