# expect_silent(library(spatstat.data))
# expect_silent(library(spatstat.geom))

expect_silent(
  forest <- RforestPP2(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 3,
    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
    cores_trees = 1,
    mtry = 1 / 3,
    minpts = 50
  )
)
expect_length(forest$trees, 3)

expect_silent(
  forest2 <- RforestPP2(
    X = spatstat.data::bei,
    listcovariates = spatstat.data::bei.extra,
    score = "lcv",
    p = 1 / 2,
    Ntree = 2,
    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
    cores_trees = 1,
    mtry = 1 / 3,
    minpts = 50
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
  outputoob <- OOBscr(forest=forest, cores=1)
)

expect_inherits(current=OOBscr(forest=forest, cores=1), class="numeric")


A <- RforestPP2(
  X = spatstat.data::bei,
  listcovariates = lapply(beisoilres, FUN=function(i){
    as.im(i, dimyx=c(10,20))
  }),
  Ntree = 50,
  minpts = 100,
  mtry = 1,
  p = 0,
  cores_trees = 2
)
B <- format(object.size(A), units = "Mb")
expect_true(as.numeric(gsub(" Mb", "",B)) < 8)




# 
# 
# x <- beisoilres[[1]]
# 
# plot(as.im(x, dimyx=c(101,201)))
# length(c(x$v))
# y<-as.im(beisoilres[[8]], dimyx=c(25,50))
# length(c(y$v))
# plot(y)
# 
# lapply(beisoilres, FUN=function(i){
#   as.im(i, dimyx=c(25,50))
# })
# 
# timer <- proc.time()
# forestnewsmall <- RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = lapply(beisoilres, FUN=function(i){
#     as.im(i, dimyx=c(10,20))
#   }),
#   Ntree = 50,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 2
# )
# timer <- proc.time() - timer; timer
# # vipplot(forestnewsmall, sorted=T, cores=5)
# # OOBscr(forestnewsmall, cores=1)
# 
# timer <- proc.time()
# forestoldsmall <- RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = lapply(beisoilres, FUN=function(i){
#     as.im(i, dimyx=c(10,20))
#   }),
#   Ntree = 50,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 2
# )
# timer <- proc.time() - timer; timer



# forestnew <- RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = beisoilres,
#   Ntree = 200,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 1
# )
# f(forestnew)
# 
# vipplot(forestnew, sorted=T)
# boxplot(forestnew)
# 
# timer <- proc.time()
# forestold <- RforestPP(
#   X = spatstat.data::bei,
#   listcovariates = beisoilres,
#   Ntree = 10,
#   minpts = 100,
#   mtry = 1 / 3,
#   p = 0,
#   cores_trees = 1
# )
# timer <- proc.time() - timer; timer
# 
# format(object.size(forestnew), units = "Mb")
# format(object.size(forestold), units = "Mb")
# 
# gold <- function(){
#   RforestPP(
#     X = spatstat.data::bei,
#     listcovariates = beisoilres,
#     Ntree = 10,
#     minpts = 100,
#     mtry = 1 / 3,
#     p = 0,
#     cores_trees = 1
#   )
# }
# 
# gnew <- function(){
#   RforestPP2(
#     X = spatstat.data::bei,
#     listcovariates = beisoilres,
#     Ntree = 10,
#     minpts = 100,
#     mtry = 1 / 3,
#     p = 0,
#     cores_trees = 1
#   )
# }
# 
# microbenchmark(gold(), gnew(), times = 10)
# bench::mark(gold(), gnew())
# 
# 
# f <- function(x){
#   format(object.size(x), units = "Mb")
# }
# 
# length(forestnew$trees)
# length(forestnew$trees[[1]]$tree)
# sapply(forestnew, FUN=function(i) f(i))
# sapply(forestnew$trees[[1]], FUN=function(i) f(i))
# 
# sapply(forestold$trees[[1]]$tree, FUN=function(i) f(i))
# forestnew$trees[[1]]$tree[[1]]
# names(forestnew$trees[[1]]$tree[[1]])
# names(forestold$trees[[1]]$tree[[1]])
# 
# 
# sapply(forestnew$trees[[1]]$tree[[1]], FUN=function(i) f(i))
# 
# f(forestnew$trees[[1]]$tree)
# 
# names(forestnew)
# f(forestnew$trees[[1]]$tree[[1]]$
# 
# 
# f(forestnew$listcov)
# 
# f(forestold$listcov)
# 
# 
# plot(forest)
# 
# 
# timer <- proc.time()
# forest <- testf(
#   X = spatstat.data::bei,
#   listcovariates = beisoilres,
#   Ntree = 100,
#   minpts = 100,
#   mtry = 1 / 3,
#   p = 0,
#   cores_trees = 1
# )
# forest
# timer <- proc.time() - timer; timer
