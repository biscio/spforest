tesstree2 <- function(X,
                      lambda) {
  # lambda is the nb of points (not the intensity)
  wind <- Window(X)
  enclose.rect <- spatstat.geom::owin(wind$xrange, wind$yrange)

  tol <- 2 / sqrt(lambda) # in order to add points outside the window
  # Alternative to avoid a polygonal windows.
  tess.points <- spatstat.geom::runifrect(
    lambda,
    owin(
      xrange = enclose.rect$xrange + c(-tol, tol),
      yrange = enclose.rect$yrange + c(-tol, tol)
    )
  )


  del <- spatstat.geom::dirichlet(tess.points) # associated tesselation
  tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows
  delarea <- spatstat.geom::tile.areas(tmp) # collect the areas

  mX <- cut(tess.points, tmp)

  ptintess <- c()
  for (i in as.numeric(names(tmp$tiles))) {
    ptintess <- c(ptintess, sum(marks(mX) == i, na.rm = TRUE))
  }


  return(as.im(tmp, values = ptintess / delarea, dimyx = c(100, 100)))
}


X <- rpoispp(function(x, y) {
  100 * exp(-3 * x)
}, 100)
lambda <- 200
B<- tesstree2(X = X, lambda = lambda)
plot(B)


B <- lapply(1:50, FUN=function(i) {
  tesstree2(X = X, lambda = lambda)})

DD <- lapply(B, FUN=function(i) {
  i$v
})
Reduce("+", DD)/50

A <- B[[1]]

A$v <- DD
plot(A)

profvis(tesstree(
  X = X, lambda = lambda,
  target.points = target.points, test.connected = F
))

microbenchmark(
  tesstree(
    X = X, lambda = lambda,
    target.points = target.points, test.connected = F
  ),
  tesstree2(
    X = X, lambda = lambda,
    target.points = target.points
  )
)


# # library(microbenchmark)
# # library(profvis)
# # A<-profvis(splitcell2(
# #      X=X,
# #      valpts = valpts,
# #      vecval = vecval,
# #      usecovariates = usecovariates,
# #      dimcov = dimcov,
# #      covrangex = covrangex,
# #      covrangey = covrangey,
# #      areapixel = areapixel,
# #      threshold = 100
# #    ))
# # A
# #
# # microbenchmark(RforestPP2(
# #   X = spatstat.data::bei,
# #   listcovariates = lapply(beisoilres, FUN=function(i){
# #     as.im(i, dimyx=c(10,20))
# #   }),
# #   Ntree = 10,
# #   minpts = 100,
# #   mtry = 1,
# #   p = 0,
# #   cores_trees = 1
# # ))
#
#
# deltirppp <- function (X)
# {
#   rw <- with(X$window, c(xrange, yrange))
#   dd <- try(deldir(X$x, X$y, rw = rw))
#   if (!inherits(dd, "try-error") && inherits(dd, "deldir"))
#     return(dd)
#   warning("deldir failed; re-trying with slight perturbation of coordinates.",
#           call. = FALSE)
#   Y <- rjitter(X, mean(nndist(X))/100)
#   dd <- try(deldir(Y$x, Y$y, rw = rw))
#   if (!inherits(dd, "try-error") && inherits(dd, "deldir"))
#     return(dd)
#   warning("deldir failed even after perturbation of coordinates.",
#           call. = FALSE)
#   return(NULL)
# }
#
#
#
# df2poly <- function (z)
# {
#   owin(poly = z[c("x", "y")])
# }
# test<-function(X){
#   stopifnot(is.ppp(X))
#   X <- unique(X, rule = "deldir", warn = TRUE)
#   nX <- npoints(X)
#   w <- X$window
#   if (nX == 0)
#     return(NULL)
#   if (nX == 1)
#     return(as.tess(w))
#   dd <- safedeldir(X)
#   if (is.null(dd))
#     return(NULL)
#   tt <- deldir::tile.list(dd)
#   pp <- lapply(tt, df2poly)
#   if (length(pp) == npoints(X))
#     names(pp) <- seq_len(npoints(X))
#   dir <- tess(tiles = pp, window = as.rectangle(w))
#   if (w$type != "rectangle")
#     dir <- intersect.tess(dir, w, keepempty = TRUE)
#   return(dir)
# }
#
# profvis(test(X))
#
# X <- rpoispp(50)
#
# library(microbenchmark)
#
#
# f <- function(){
#   deldir::deldir(tess.points$x, tess.points$y, rw= c(0, 1, 0, 1))
# }
# g <- function(){
#   del <- spatstat.geom::dirichlet(tess.points) # associated tesselation
#   tmp <- spatstat.geom::intersect.tess(del, wind)
#   tmp
# }
#
# plot(f())
# plot(square(1),add=T)
# plot(g())
#
# microbenchmark(f(),g())
#
# deldir::deldir(X$x, X$y,plot = T)
#
# del <- spatstat.geom::dirichlet(tess.points)
#
# deldir::deldir(tess.points$x, tess.points$y)
#
# plot(spatstat.geom::dirichlet(X))
#
# profvis(spatstat.geom::dirichlet(X))
#
# #
# #
# # A<-RforestPP2(
# #   X = spatstat.data::bei,
# #   listcovariates = beisoilres,
# #   Ntree = 10,
# #   minpts = 100,
# #   mtry = 1,
# #   p = 0,
# #   cores_trees = 1
# # )
# #
# # plot(A)
# #
# # microbenchmark(RforestPP2(
# #   X = spatstat.data::bei,
# #   listcovariates = lapply(beisoilres, FUN=function(i){
# #     as.im(i, dimyx=c(10,20))
# #   }),
# #   Ntree = 10,
# #   minpts = 100,
# #   mtry = 1,
# #   p = 0,
# #   cores_trees = 1
# # ),
# # RforestPP(
# #   X = spatstat.data::bei,
# #   listcovariates = lapply(beisoilres, FUN=function(i){
# #     as.im(i, dimyx=c(10,20))
# #   }),
# #   Ntree = 10,
# #   minpts = 100,
# #   mtry = 1,
# #   p = 0,
# #   cores_trees = 1
# # ))
# # #
# # # library(profvis)
# # #
# # # profvis(RforestPP2(
# # #   X = spatstat.data::bei,
# # #   listcovariates = lapply(beisoilres, FUN=function(i){
# # #     as.im(i, dimyx=c(10,20))
# # #   }),
# # #   Ntree = 1,
# # #   minpts = 100,
# # #   mtry = 1,
# # #   p = 0,
# # #   cores_trees = 1
# # # ))
# # #
# # # profvis(splitcell2(
# # #   X=X,
# # #   valpts = valpts,
# # #   vecval = vecval,
# # #   usecovariates = usecovariates,
# # #   dimcov = dimcov,
# # #   covrangex = covrangex,
# # #   covrangey = covrangey,
# # #   areapixel = areapixel,
# # #   threshold = 100
# # # ))
# # #
# # # # ftree <- function(){
# # # #   treerec(
# # # #     X = spatstat.data::bei,
# # # #     listcovariates = Rsandbox::beisoilres,
# # # #     mtry = 1,
# # # #     minpts = 500
# # # #   )
# # # # }
# # # #
# # # # gtree <- function(){
# # # #   intensitytree(
# # # #     X = spatstat.data::bei,
# # # #     listcovariates = Rsandbox::beisoilres,
# # # #     mtry = 1,
# # # #     minpts = 500
# # # #   )
# # # # }
# # # #
# # # # A<-ftree()
# # # # B <- gtree()
# # # # plot(A)
# # # # plot(B)
# # # # test<-A$im-B$im ### Small difference why ??
# # # # max(abs(test))
# # # #
# # # # sapply(intensity_tree2, FUN = function(i) {
# # # #   i$right_daughter
# # # # })
# # # # sapply(intensity_tree2, FUN = function(i) {
# # # #   i$left_daughter
# # # # })
# # # # sapply(intensity_tree2, FUN = function(i) {
# # # #   i$status
# # # # })
# # # # sapply(intensity_tree2, FUN = function(i) {
# # # #   i$already_split
# # # # })
# # # # sapply(B$tree, FUN = function(i) {
# # # #   i$intensity_pred
# # # # })
# # # # sapply(A$tree, FUN = function(i) {
# # # #   i$intensity_pred
# # # # })
# # # #
# # # # microbenchmark(ftree(), gtree(), times = 5)
# # # #
# # # # library(profvis)
# # # # profvis(ftree())
# # # # profvis(gtree())
# # # #
# # # #
# # # #
# # # # f <- function(){
# # # #   RforestPP(
# # # #     X = spatstat.data::bei,
# # # #     listcovariates = Rsandbox::beisoilres,
# # # #     Ntree = 10,
# # # #     minpts = 500,
# # # #     mtry = 1,
# # # #     cores_trees = 5
# # # #   )
# # # # }
# # # #
# # # # g <- function(){
# # # #   RforestPP2(
# # # #     X = spatstat.data::bei,
# # # #     listcovariates = Rsandbox::beisoilres,
# # # #     Ntree = 10,
# # # #     minpts = 500,
# # # #     mtry = 1,
# # # #     cores_trees = 5
# # # #   )
# # # # }
# # # #
# # # # library(microbenchmark)
# # # # microbenchmark(f(), g(), times = 5)
# # # #
