# tesstree2 <- function(X,
#                       lambda = 100,
#                       dimyx = c(128, 128),
#                       test.connected = FALSE) {
#   # lambda is the nb of points (not the intensity)
#   wind <- Window(X)
#   enclose.rect <- spatstat.geom::owin(wind$xrange, wind$yrange)
# 
#   tol <- 2 / sqrt(lambda) # in order to add points outside the window
#   # Alternative to avoid a polygonal windows.
#   tess.points <- spatstat.geom::runifrect(
#     lambda,
#     owin(
#       xrange = enclose.rect$xrange + c(-tol, tol),
#       yrange = enclose.rect$yrange + c(-tol, tol)
#     )
#   )
# 
#   del <- spatstat.geom::dirichlet(tess.points) # associated tessellation
#   tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows
#   if (test.connected) {
#     tmp <- spatstat.geom::connected(tmp)
#   }
# 
#   delarea <- spatstat.geom::tile.areas(tmp) # collect the areas
# 
#   # mX <- marks(cut(X, tmp)) ## the alternative is very slightly quicker
#   mX <- tileindex(X$x, X$y, tmp)
# 
#   ptintess <- c()
#   if (test.connected) {
#     for (i in levels(tmp$image)) {
#       ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
#     }
#   } else {
#     for (i in names(tmp$tiles)) {
#       ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
#     }
#   }
#   return(as.im(tmp, values = ptintess / delarea, dimyx = dimyx))
# }
# 
# tessforest2 <- function(X,
#                         Ntree = 1,
#                         lambda = 100,
#                         dimyx = c(50, 50),
#                         test.connected = FALSE,
#                         cores = 1) {
#   if (is.null(lambda)) {
#     lambda <- floor(mean(c(
#       grDevices::nclass.FD(X$x),
#       grDevices::nclass.FD(X$y)
#     ))^2)
#   }
# 
#   if (cores > 1) {
#     listtree <- parallel::mclapply(1:Ntree, FUN = function(i) {
#       tesstree2(
#         X = X,
#         lambda = lambda,
#         dimyx = dimyx,
#         test.connected = test.connected
#       )
#     }, mc.cores = cores)
#   } else {
#     listtree <- lapply(1:Ntree, FUN = function(i) {
#       tesstree2(
#         X = X,
#         lambda = lambda,
#         dimyx = dimyx,
#         test.connected = test.connected
#       )
#     })
#   }
# 
#   return(Reduce("+", listtree) / length(listtree))
# }
# 
# 
# load("~/Downloads/anuraPP.RData")
# load("~/Downloads/anura.RData")
# plot(anura)
# 
# X <- rthin(anura, 0.1)
# 
# plot(X)
# A <- tesstree2(X = X, lambda = 100, dimyx = c(128, 128), test.connected = T)
# timer <- proc.time()
# B <- tessforest2(X = X, Ntree = 10, 
#                  lambda = 100, dimyx = c(128, 128), 
#                  test.connected = F)
# proc.time() - timer
# timer2 <- proc.time()
# Bpar <- tessforest2(X = X, Ntree = 10, 
#                     lambda = 100, dimyx = c(128, 128), 
#                     test.connected = F, cores = 5)
# proc.time() - timer2
# 
# plot(A)
# plot(B)
# plot(anura)
# 
# A <- tessforest2(X = X, lambda = 100, dimyx = c(128, 128), Ntree = 15)
# format(object.size(A), unit = "Mb")
# 
# 
# f <- function() {
#   tessforest2(X = X, lambda = 100, dimyx = c(128, 128), Ntree = 5)
# }
# floop <- function() {
#   tessforest3(X = X, lambda = 100, dimyx = c(128, 128), Ntree = 5)
# }
# 
# g <- function() {
#   tessforest(X = X, lambda = 100, Ntree = 5, test.connected = FALSE)
# }
# 
# microbenchmark::microbenchmark(f(), g())
# microbenchmark::microbenchmark(f(), floop(), times = 10)
# 
# library(bench)
# bench::mark(f())
# bench::mark(g())
# 
# profvis::profvis(f())
# 
# points(X[1], col = 1, pch = 16, cex = 0.8)
# 
# points(X, cex = 0.8, col = 4)
# points(m, cex = 1, col = 2, )
# marks(m)
# library(microbenchmark)
# 
# microbenchmark(marks(cut(X, dd)), tileindex(X$x, X$y, dd))
# 
# X <- rpoispp(function(x, y) {
#   3000 * exp(-3 * x)
# })
# 
# II <- as.im(function(x, y) {
#   3000 * exp(-3 * x)
# }, W = square(1))
# lambda <- 100
# B <- tesstree2(X = X, lambda = lambda)
# plot(II)
# plot(B)
# 
# 
# B <- lapply(1:50, FUN = function(i) {
#   tesstree2(X = X, lambda = lambda)
# })
# 
# plot(Reduce("+", B) / 50)
# library(bench)
# 
# profvis(tesstree(
#   X = X, lambda = lambda,
#   target.points = target.points, test.connected = F
# ))
# 
# microbenchmark(
#   tesstree(
#     X = X, lambda = lambda,
#     target.points = target.points, test.connected = F
#   ),
#   tesstree2(
#     X = X, lambda = lambda,
#     target.points = target.points
#   )
# )
