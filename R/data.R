#' @title Raw soil data for bei
#'
#' @description Contains the raw soil data of bei
#'
#' @format A list of 15 \code{im} with different resolutions, and given in the
#' same order as in Jeff's paper
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"beisoil"

#' @title Smooth soil data for bei
#'
#' @description Contains the soil data of bei put at the same resolution and smoothed 
#' by kriging.
#' It has been obtained with the following code. One need to download the bci 
#' data on \link{http://ctfs.si.edu/webatlas/datasets/bci/soilmaps/BCIsoil.html}.
#' library(gstat)
#' library(sp)
#' 
#' bci.blocks20 <- read_xls("bci.block20.data.xls", sheet = 2) |> 
#' rename("Nmin" = "N(min)")
#' 
#' boxcox <- function(x, lambda = 1) {
#'   if (lambda == 0) {
#'     return(log(x))
#'   } else {
#'     return((x^lambda - 1) / lambda)
#'   }
#' }
#' 
#' boxcoxinv <- function(x, lambda = 1) {
#'   if (lambda == 0) {
#'     return(exp(x))
#'   } else {
#'     return((lambda * x + 1)^(1 / lambda))
#'   }
#' }
#' 
#' smoothbci <- function(z,
#'                       lambda = 1,
#'                       degreg = 2,
#'                       model = "Exp",
#'                       nugget = 0.9552,
#'                       psill = 4.54,
#'                       effrange = 179.1,
#'                       block = c(20, 20)) {
#'   mod <- lm(boxcox(bci.blocks20[[z]]) ~ poly(x, y, degree = degreg),
#'     data = bci.blocks20
#'   )
#'   spatialres <- data.frame(
#'     x = bci.blocks20$x,
#'     y = bci.blocks20$y,
#'     z = residuals(mod)
#'   )
#'   g <- gstat(
#'     data = spatialres,
#'     formula = z ~ 1,
#'     locations = ~ x + y
#'   )
#'   coordinates(spatialres) <- ~ x + y
#' 
#'   res.fit <- fit.variogram(variogram(g),
#'     model = vgm(
#'       psill = psill,
#'       model = model,
#'       range = effrange,
#'       nugget = nugget
#'     )
#'   )
#' 
#'   ## Create new grid
#'   newseqx <- seq(0, 1000, by = 5)
#'   newseqy <- seq(0, 500, by = 5)
#'   newgrid <- expand.grid(x = newseqx, y = newseqy)
#'   gridded(newgrid) <- ~ x + y
#' 
#'   res.kriged <- krige(
#'     formula = z ~ 1,
#'     spatialres,
#'     newgrid,
#'     model = res.fit,
#'     block = block
#'   )
#' 
#'   newval <- res.kriged[["var1.pred"]] + predict(mod, as.data.frame(newgrid))
#' 
#'   mat <- matrix(boxcoxinv(newval),
#'     nrow = length(newseqx),
#'     ncol = length(newseqy),
#'     byrow = F
#'   )
#' 
#'   return(as.im(t(mat), W = as.owin(bei.extra[[1]])))
#' }
#' 
#' argu <- data.frame(
#'   z = c(
#'     "Al", "B", "Ca", "Cu", "Fe", "K",
#'     "Mg", "Mn", "P", "Zn", "N", "Nmin", "pH"
#'   ),
#'   lambda = c(1, 0.5, 0.5, 1, 0.5, 0, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 1),
#'   degreg = c(rep(2, 11), 1, 2),
#'   model = c(rep("Exp", 10), "Sph", "Sph", "Exp"),
#'   nugget = c(
#'     10801.7, 0.04053, 21.29, 0.9552, 2.95,
#'     0.0875, 5.73, 2.846, 0.4556, 0.00662,
#'     0.0754, 6.23, 0.0172
#'   ),
#'   psill = c(
#'     68415.3, 0.3649, 425.6, 4.54, 15.83,
#'     0.2193, 64.27, 133.5, 1.471, 2.135, 0.0986, 2.192, 0.1403
#'   ),
#'   effrange = c(
#'     155.0, 100.0, 56.5, 179.1, 146.6,
#'     78.5, 85.2, 212.8, 204.5, 33.3, 239.5, 141.1, 220.7
#'   )
#' )
#' 
#' tib <- purrr::pmap(.l = argu, .f = smoothbci)
#' 
#' bcismooth <- append(bei.extra, tib)
#' names(bcismooth) <- c("elev", "grad", argu$z)
#' @format A list of 15 \code{im} with same resolutions, and given in the
#' same order as in Jeff's paper
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"bcismooth"


#' @title soil data for bei
#'
#' @description Contains the soil data of bei put at the same resolution
#' It has been obtained with the following code
# beisoil <- list(
#   bei.extra$elev, bei.extra$grad,
#   Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, Nmin, pH
# )
#
# ## Put them in the same resolution ----
# beisoilres <- NULL
# beisoilres[[1]] <- beisoil[[1]]
# beisoilres[[2]] <- beisoil[[2]]
# for (i in 3:15) {
#   W <- commonGrid(beisoil[[1]], beisoil[[i]])
#   beisoilres[[i]] <- beisoil[[i]][W, drop = F]
#   rm(W)
# }
# names(beisoilres) <- c(
#   "elev", "grad", "Al", "B", "Ca",
#   "Cu", "Fe", "K", "Mg", "Mn",
#   "P", "Zn", "N", "Nmin", "pH"
# )
#
# ## center and scale
# beisoilnorm <- lapply(beisoilres, FUN = function(i) {
#   (i - mean(i)) / sd(i)
# })
#
# names(beisoilnorm) <- c(
#   "elev", "grad", "Al", "B", "Ca",
#   "Cu", "Fe", "K", "Mg", "Mn",
#   "P", "Zn", "N", "Nmin", "pH"
# )
#'
#' @format A list of 15 \code{im} with same resolutions, and given in the
#' same order as in Jeff's paper
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"beisoilres"

#' @title Downscaled soil data for bei
#'
#' @description Contains a downscaled version of the data for testing and vignette
#'
#' @format A list of 15 \code{im} with same resolutions, and given in the
#' same order as in Jeff's paper
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"beisoilsmall"

#' @title Extra Downscaled soil data for bei
#'
#' @description Contains an extra downscaled version of the data for testing and vignette
#'
#' @format A list of 15 \code{im} with same resolutions, and given in the
#' same order as in Jeff's paper
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"beisoilxsmall"

#' @title Normalised soil data for bei
#'
#' @description Contains the normalised soil data of bei
#'
#' @format A list of 15 \code{im} with same resolutions, in the
#' same order as in Jeff's paper, and center and normalised
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"beisoilnorm"

#' @title Data on position on trout
#'
#' @description Toy datasets
#'
#' @format TODO
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"trout"

#' @title Data for testing purposes
#'
#' @description Toy datasets
#'
#' @format TODO
#' \describe{
#'     \item{Al}{Al concentration ...}
#' }
#'
#'
"Xtest"


#' @title Simulated point pattern on a 3D mesh
#'
#' @description Contains a simulated inhomogeneous point pattern simulated on 
#' the mesh humface from the R package Rvcg.
#'library(Rvcg)
#'library(rgl)
#'library(spforest)
#'data("humface")
#'
#'features <- features_mesh(humface)
#'tricenter <- features$tricenter
#'nearest <- RANN::nn2(tricenter, rbind(humface$vb[, 7440][1:3]), k = 2000) # k = 1 => plus proche
#'meshspot <- humface
#'meshspot$it <- meshspot$it[, nearest$nn.idx]
#'
#'set.seed(28)
#'X <- dummypponmesh(mesh = meshspot, n = 300, weights = FALSE)
#'Y <- dummypponmesh(mesh = humface, n = 850, weights = FALSE)
#'Z <- Y
#'Z$pp <- rbind(X$pp, Y$pp)
#'
#'g <- function(pts) {
#'  (0.2 * pts[1] - 0.9 * pts[2] + pts[3] + 70) < 0
#'}
#'
#'idside <- sapply(1:dim(Z$pp)[1], FUN = function(i) {
#'  g(Z$pp[i, ])
#'})
#'
#'simppface <- list(mesh = Z$mesh, pp = Z$pp[!idside,])
#'
#'save(simppface, file = "simppface.rda", compress = "xz")
#' @format A list of two elements: "mesh" and "pp". 
#' 
#'
#'
"simppface"

