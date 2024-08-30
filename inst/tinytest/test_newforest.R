# 
# 

x <- beisoilres[[1]]

plot(as.im(x, dimyx=c(101,201)))
length(c(x$v))
y<-as.im(beisoilres[[8]], dimyx=c(25,50))
length(c(y$v))
plot(y)

lapply(beisoilres, FUN=function(i){
  as.im(i, dimyx=c(25,50))
})

timer <- proc.time()
forestnewsmall <- RforestPP2(
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
timer <- proc.time() - timer; timer
# vipplot(forestnewsmall, sorted=T, cores=5)
# OOBscr(forestnewsmall, cores=1)

timer <- proc.time()
forestoldsmall <- RforestPP(
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
timer <- proc.time() - timer; timer



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
