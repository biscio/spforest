#' Intensity tessellation tree
#'
#' @param X The observed data point pattern,
#' as a \code{\link[spatstat.geom]{ppp}} object .
#' @param lambda An integer. The number of points used for random tessellation.
#' @param dimyx A vector of two integers. The dimensions of the output image passed
#' to \code{\link[spatstat.geom]{as.im}}.
#' @param test.connected Logical. If \code{TRUE},
#' \code{\link[spatstat.geom]{connected}} is applied to the tessellation to split tiles
#' in different connected components. It is only useful if the windows of
#' observation od \code{X} is not convex.
#'
#' @return A pixel image, object of class \code{\link[spatstat.geom]{im.object}}.
#' @export
#'
#' @examples
#' Z <- tesstree2(
#'   X = bei,
#'   lambda = 100,
#'   dimyx = c(101, 201),
#'   test.connected = FALSE
#' )
#' plot(Z)
tesstree2 <- function(X,
                      lambda = 100,
                      dimyx = c(128, 128),
                      test.connected = FALSE) {
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

  del <- spatstat.geom::dirichlet(tess.points) # associated tessellation
  tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows
  if (test.connected) {
    tmp <- spatstat.geom::connected(tmp)
  }

  delarea <- spatstat.geom::tile.areas(tmp) # collect the areas

  # mX <- marks(cut(X, tmp)) ## the alternative is very slightly quicker
  mX <- tileindex(X$x, X$y, tmp)

  ptintess <- c()
  if (test.connected) {
    for (i in levels(tmp$image)) {
      ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
    }
  } else {
    for (i in names(tmp$tiles)) {
      ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
    }
  }
  return(as.im(tmp, values = ptintess / delarea, dimyx = dimyx))
}



#' Intensity tessellation forest
#'
#' @inheritParams tesstree2
#' @param Ntree A positive integer.
#' The number of trees in the random intensity forest.
#' @param cores A positive integer.
#' The number of cores used to computes the intensity trees.
#'
#' @return A pixel image, object of class \code{\link[spatstat.geom]{im.object}}.
#' @export
#'
#' @examples
#' Z <- tessforest2(
#'   X = bei,
#'   lambda = 100,
#'   dimyx = c(101, 201),
#'   test.connected = FALSE,
#'   Ntree = 5,
#'   cores = 1
#' )
#' plot(Z)
tessforest2 <- function(X,
                        Ntree = 1,
                        lambda = 100,
                        dimyx = c(50, 50),
                        test.connected = FALSE,
                        cores = 1) {
  if (is.null(lambda)) {
    lambda <- floor(mean(c(
      grDevices::nclass.FD(X$x),
      grDevices::nclass.FD(X$y)
    ))^2)
  }

  if (cores > 1) {
    listtree <- parallel::mclapply(1:Ntree, FUN = function(i) {
      tesstree2(
        X = X,
        lambda = lambda,
        dimyx = dimyx,
        test.connected = test.connected
      )
    }, mc.cores = cores)
  } else {
    listtree <- lapply(1:Ntree, FUN = function(i) {
      tesstree2(
        X = X,
        lambda = lambda,
        dimyx = dimyx,
        test.connected = test.connected
      )
    })
  }

  output <- list(
    imforest = Reduce("+", listtree) / length(listtree),
    trees = NULL,
    ntrees = length(listtree),
    pt_intree = rep(1, X$n),
    X = X,
    listcov = NULL,
    p = NULL,
    mtry = NULL
  )
  class(output) <- "spforest"
  
  return(output)
}
