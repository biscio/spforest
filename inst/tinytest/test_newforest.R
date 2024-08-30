forestnew <- RforestPP2(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  Ntree = 10,
  minpts = 100,
  mtry = 1 / 3,
  p = 0,
  cores_trees = 1
)

timer <- proc.time()
forestold <- RforestPP(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  Ntree = 10,
  minpts = 100,
  mtry = 1 / 3,
  p = 0,
  cores_trees = 1
)
timer <- proc.time() - timer; timer

format(object.size(forestnew), units = "Mb")
format(object.size(forestold), units = "Mb")

gold <- function(){
  RforestPP(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    Ntree = 10,
    minpts = 100,
    mtry = 1 / 3,
    p = 0,
    cores_trees = 1
  )
}

gnew <- function(){
  RforestPP2(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    Ntree = 10,
    minpts = 100,
    mtry = 1 / 3,
    p = 0,
    cores_trees = 1
  )
}

microbenchmark(gold(), gnew(), times = 10)
bench::mark(gold(), gnew())


f <- function(x){
  format(object.size(x), units = "Mb")
}

length(forestnew$trees)
length(forestnew$trees[[1]]$tree)
sapply(forestnew, FUN=function(i) f(i))
sapply(forestnew$trees[[1]], FUN=function(i) f(i))

sapply(forestold$trees[[1]]$tree, FUN=function(i) f(i))
forestnew$trees[[1]]$tree[[1]]
names(forestnew$trees[[1]]$tree[[1]])
names(forestold$trees[[1]]$tree[[1]])


sapply(forestnew$trees[[1]]$tree[[1]], FUN=function(i) f(i))

f(forestnew$trees[[1]]$tree)

names(forestnew)
f(forestnew$trees[[1]]$tree[[1]]$)


f(forestnew$listcov)

f(forestold$listcov)


plot(forest)


timer <- proc.time()
forest <- testf(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  Ntree = 100,
  minpts = 100,
  mtry = 1 / 3,
  p = 0,
  cores_trees = 1
)
forest
timer <- proc.time() - timer; timer
