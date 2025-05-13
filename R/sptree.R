#' sptree.object
#'
#' S3 class used to define random intensity trees.
#'
#' @name sptree.object
#' @rdname sptree.object
#' 
#' @details 
#' An object of this class represents a spatial intensity tree.
#' It contains 
#' \itemize{
#' \item \code{tree}, a list of nodes of the tree. The \eqn{i}-th 
#' element contains a list with entries \code{nodeID, nX, left_daughter, 
#' right_daughter, split_var, split_val, whystop}.
#' \item \code{X}, the point pattern used to estimate the tree.
#' \item{namecov}, the name of the covariates used in 
#' \code{listcovariates}.
#' \item{namelist}, TODO.
#' \item \code{listcov}, the list of covariates used.
#' \item \code{im}, a pixel image, as \code{\link[spatstat.geom]{im.object}},
#'  representing the estimated intensity.
#' }
NULL


#' Printing spatial intensity tree
#'
#' @param x A spatial intensity tree return by tesstree function
#' @param ... Additional arguments
#'
#' @return A description of the random intensity tree with
#' the number of points and covariates used.
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' arbre <- tesscovtree(
#'   X = spatstat.data::bei,
#'   vecval = vecval0,
#'   areapixel = beisoilres[[1]]$xstep * beisoilres[[1]]$ystep,
#'   dimcov = beisoilres[[1]]$dim,
#'   covrangex = beisoilres[[1]]$xrange,
#'   covrangey = beisoilres[[1]]$yrange,
#'   listcovariates = beisoilres,
#'   mtry = 1,
#'   minpts = 500
#' )
#' print(arbre)
print.sptree <- function(x, ...) {
  namecov <- x$namecov
  nb_termnode <- sum(sapply(x$tree, function(i) {
    i$status
  }) == 0)

  cat(paste(
    "Intensity tree estimate of point patterns with",
    x$X$n, "points.\n\n"
  ))

  cat(paste(length(namecov), "covariates used, with names: "))
  cat(namecov, "\n")

  cat(
    "Spatial intensity tree with", length(x$tree),
    "nodes and", nb_termnode, "terminal nodes.", "\n"
  )
  cat("The used covariates were in the object TODO", ".")
}

#' Summary of a spatial intensity tree
#'
#' @param object A spatial intensity tree return by tesstree function
#' @param fulltree Should we print all the column of the matrix?
#' @param ... Additional arguments
#'
#' @return A dataframe representing each node and spit of the intensity tree.
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' arbre <- tesscovtree(
#'   X = spatstat.data::bei,
#'   vecval = vecval0,
#'   areapixel = beisoilres[[1]]$xstep * beisoilres[[1]]$ystep,
#'   dimcov = beisoilres[[1]]$dim,
#'   covrangex = beisoilres[[1]]$xrange,
#'   covrangey = beisoilres[[1]]$yrange,
#'   listcovariates = beisoilres,
#'   mtry = 1,
#'   minpts = 500
#' )
#' summary(arbre)
summary.sptree <- function(object, fulltree = F, ...) {
  if (fulltree) {
    output <- Reduce(rbind, object$tree)
    rownames(output) <- NULL
  } else {
    # Just remove the imp and nodeID
    output <- Reduce(rbind, object$tree)[, c(2, 3, 4, 5, 6, 7, 8)]
    rownames(output) <- NULL
  }
  output
}



#' Plot spatial intensity tree
#'
#' @param x A spatial intensity tree return by tesstree function
#' @param ... additional arguments
#' @param main Title of the plot.
#'
#' @return A plot of the spatial intensity tree
#' and an \code{\link[spatstat.geom]{im}} object.
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' arbre <- tesscovtree(
#'   X = spatstat.data::bei,
#'   vecval = vecval0,
#'   areapixel = beisoilres[[1]]$xstep * beisoilres[[1]]$ystep,
#'   dimcov = beisoilres[[1]]$dim,
#'   covrangex = beisoilres[[1]]$xrange,
#'   covrangey = beisoilres[[1]]$yrange,
#'   listcovariates = beisoilres,
#'   mtry = 1,
#'   minpts = 500
#' )
#' plot(arbre)
plot.sptree <- function(x, ..., main = "Spatial Intensity Tree") {
  # Handling case if no main title is given for the plot
  # if (missing(main)) {
  #   main <- "Spatial Intensity Tree"
  # }

  spatstat.geom::plot.im(x$im, main = main, ...)

  return(invisible(x$im))
}


#' Tree prediction
#'
#' @param object A spatial intensity tree returned 
#' by tesscovtree of tesstree functions
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @return A number .....
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' arbre <- tesscovtree(
#'   X = spatstat.data::bei,
#'   vecval = vecval0,
#'   areapixel = beisoilres[[1]]$xstep * beisoilres[[1]]$ystep,
#'   dimcov = beisoilres[[1]]$dim,
#'   covrangex = beisoilres[[1]]$xrange,
#'   covrangey = beisoilres[[1]]$yrange,
#'   listcovariates = beisoilres,
#'   mtry = 1,
#'   minpts = 500
#' )
#' predict(object = arbre, newdata = c(100, 100))
predict.sptree <- function(object, newdata, ...) {
  # Handles the newdata to be in the correct form
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else if (!spatstat.geom::is.ppp(newdata)) {
    X <- spatstat.geom::ppp(
      x = newdata[1],
      y = newdata[2],
      window = object$X$window
    )
  } else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n == 0) {
      return(NULL)
    }
    X <- newdata
  }

  return(object$im[X])
}

