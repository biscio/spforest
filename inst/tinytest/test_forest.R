expect_silent(library(spatstat.data))
expect_silent(library(spatstat.geom))

expect_silent(
  forest <- RforestPP(
    X = spatstat.data::bei,
    listcovariates = list(
      grad = spatstat.data::bei.extra$grad,
      elev = spatstat.data::bei.extra$elev
    ),
    score = "lcv",
    p = 1 / 2,
    Ntree = 3,
    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
    cores_trees = 1,
    mtry = 1 / 3,
    tol = Inf,
    minpts = 50,
    minsplitq = 0.5,
    maxsplitq = 0.5
  )
)

# expect_silent(validate_spforest(x=forest))

expect_silent(
    forest2 <- RforestPP(
        X = spatstat.data::bei,
        listcovariates = list(
            grad = spatstat.data::bei.extra$grad,
            elev = spatstat.data::bei.extra$elev
        ),
        score = "lcv",
        p = 1 / 2,
        Ntree = 2,
        threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
        cores_trees = 1,
        mtry = 1 / 3,
        tol = Inf,
        minpts = 50,
        minsplitq = 0.5,
        maxsplitq = 0.5
    )
)


expect_silent(merge.spforest(x=forest, y=forest2))

expect_length(forest$trees, 3)

expect_true(length(forest$trees[[3]]$tree)>0)

expect_equal(forest$X, spatstat.data::bei)

expect_equal(sum(forest$listcov[[1]] - spatstat.data::bei.extra$grad), 0)
expect_equal(sum(forest$listcov[[2]] - spatstat.data::bei.extra$elev), 0)

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
    outputoob <- OOBscr(forest=forest, cores=3)
)

expect_inherits(current=OOBscr(forest=forest, cores=3), class="numeric")







