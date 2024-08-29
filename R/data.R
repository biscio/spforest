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


#' @title soil data for bei
#'
#' @description Contains the soil data of bei put at the same resolution
#'It has been obtained with the following code 
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
