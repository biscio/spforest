#' Spatial Intensity tree function without covariate
#'
#' @param X The observed data point pattern,
#' as a \code{\link[spatstat.geom]{ppp.object}}.
#' @param gamma An integer. The number of points used to generate random tessellation.
#' @param dimyx A vector of two integers representing the dimensions of the output image passed
#' to \code{\link[spatstat.geom]{as.im}}.
#' @param test.connected Logical. If \code{TRUE},
#' \code{\link[spatstat.geom]{connected}} is applied to the tessellation to split tiles
#' in different connected components. It is only useful if the windows of
#' observation of \code{X} is not convex.
#'
#' @details
#' The function generates a spatial intensity
#' tree from a point pattern \code{X} using an independently
#' generated Voronoi tessellation.
#' If the point pattern \code{X} is on a rectangular window,
#' the algorithm generates a Voronoi tessellation from
#' \code{gamma} points
#' uniformly sampled on the window of \code{X}, fattened by
#' \eqn{2 / \sqrt{\gamma}} to avoid border effect.
#' The intensity tree at coordinates \eqn{(x, y)}
#' is the number of points of \code{X} in the tile
#' of the tessellation containing \eqn{(x, y)}
#' divided by the area of the tile.
#'
#' If the point pattern \code{X} is on a non-rectangular window,
#' the algorithm generates a Voronoi tessellation
#' on the enclosing rectangle of the window of \code{X},
#' fattened by \eqn{2 / \sqrt{\gamma}}.
#' The tessellation is then intersected
#' with the window of \code{X} and,
#' if \code{test.connected = TRUE},
#' the function \code{\link[spatstat.geom]{connected}} is used
#' to split non-connected cells.
#'
#' @return A pixel image, object of class \code{\link[spatstat.geom]{im.object}}.
#' @export
#'
#' @examples
#' Z <- tesstree(
#'   X = bei,
#'   gamma = 100,
#'   dimyx = c(101, 201),
#'   test.connected = FALSE
#' )
#' plot(Z)
tesstree <- function(X,
                     gamma = 100,
                     dimyx = c(128, 128),
                     test.connected = FALSE) {
  # gamma is the nb of points (not the intensity)
  wind <- Window(X)
  enclose.rect <- spatstat.geom::owin(wind$xrange, wind$yrange)

  tol <- 2 / sqrt(gamma) # in order to add points outside the window
  # Alternative to avoid a polygonal windows.
  tess.points <- spatstat.geom::runifrect(
    gamma,
    owin(
      xrange = enclose.rect$xrange + c(-tol, tol),
      yrange = enclose.rect$yrange + c(-tol, tol)
    )
  )

  del <- spatstat.geom::dirichlet(tess.points) # associated tessellation
  tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows
  tmp0 <- tmp
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
  return(list(
    intensityim = as.im(tmp, values = ptintess / delarea, dimyx = dimyx),
    intensitytess = tmp0
  ))
}



#' Spatial Intensity tree function with covariates
#'
#' @param X The observed data point pattern.
#' @param vecval A list. The i-th element is the values of
#' a i-th covariate in \code{listcovariates} at the points of \code{X}.
#' @param areapixel The pixel area used for all covariate.
#' @param dimcov The element \code{dim} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#' @param covrangex The element \code{xrange} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#' @param covrangey The element \code{yrange} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#' @param listcovariates A list with all the covariates used for the tree. The
#' covariates must be given as im object from the package spatstat.
#' @param minpts A positive integer.
#' The minimum number of points allowed to try to split a cell.
#' @param mtry Probability of choosing a covariate.
#' @param randmtry Logical. If \code{TRUE}, \code{mtry} must be between 0 and 1 and
#' represents the probability to use each covariate at each split. If \code{FALSE}, \code{mtry}
#' covariates are randomly chosen at each split.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param threshold Minimum threshold to allow to split cell.
#' @param inforest Logical. Indicates if the function is run in a forest or not.
#'
#' @details
#' The function generates a spatial intensity
#' tree from a point pattern \code{X} using a tessellation
#' build from covariates. The covariates needs to be
#' input in a vectorised form in \code{vecval}.
#'
#' The tessellation is build iteratively, starting from
#' full window of \code{X} as the first cell, until there is less
#' than \code{minpts} points in a cell or
#' the area of the cell is less that \code{threshold}.
#'
#' At each step, we look to split the cell
#' according to the sub and sup level sets of
#' the vectorised covariate in \code{vecval}.
#' We choose among them the split that
#' gives the maximal score, given by \code{score}.
#' Moreover, at each step, only a randomly subset of
#' to covariates, each selected with probability \code{mtry}
#' are used.
#'
#' The arguments \code{dimcov}, \code{covrangex} and \code{covrangey},
#' are passed to \code{\link[spatstat.geom]{im}} to return
#' a pixel image of the spatial intensity tree.
#'
#' The function output an object of class  \code{sptree} which
#' contains the point pattern \code{X}.
#' This is useful for information purposes.
#' However, if the function is run
#' as part of an intensity tessellation forest,
#' it uses a lot of memory to create many copies of \code{X}.
#' To avoid this, set \code{inforest=T}.
#'
#' @return A pixel image as
#' \code{\link[spatstat.geom]{im.object}}.
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' mytree <- tesscovtree(
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
#' plot(mytree)
tesscovtree <- function(X,
                        vecval,
                        areapixel,
                        dimcov,
                        covrangex,
                        covrangey,
                        listcovariates,
                        minpts = 500,
                        mtry = 1,
                        randmtry = FALSE,
                        score = "lcv",
                        threshold = spatstat.geom::area(X) / 1e4,
                        inforest = F) {
  valpts <- lapply(listcovariates,
    FUN = function(i) {
      i[X]
    }
  )

  root <- list(
    nodeID = 1,
    nodeCov = vecval,
    nodeValpts = valpts,
    left_daughter = NA,
    right_daughter = NA,
    nX = spatstat.geom::npoints(X),
    split_var = NA,
    split_val = NA,
    status = 1,
    intensity_pred = spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
    already_split = FALSE,
    whystop = NULL,
    scrsplit = NA,
    scrdcr = NA
  )

  # Check if there is something to do on the initial point pattern.
  if (spatstat.geom::area.owin(X$window) <= threshold |
    spatstat.geom::npoints(X) <= minpts) {
    root$status <- 0
    output <- list(
      tree = list(root),
      X = X,
      namecov = names(listcovariates),
      im = spatstat.geom::as.im(spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
        W = X$window
      )
    )
    class(output) <- "sptree"
    return(output)
  }

  intensity_tree <- list(root)

  # Initialise the while loop
  k <- 0
  knew <- 1


  while (k != knew) {
    k <- length(intensity_tree)

    already_split_node <- sapply(intensity_tree,
      FUN = function(j) {
        j$already_split
      }
    )

    for (i in (1:k)[!already_split_node]) {
      ## Select randomly covariates
      usedcov <- rand_covar(
        listcovariates = listcovariates,
        mtry = mtry,
        randmtry = randmtry
      )

      if (intensity_tree[[i]]$nX <= minpts) {
        res.split <- "Not enough points to attempt to split"
      } else {
        # Split the cell, if the split is valid under the chosen parameters

        if (score == "palm") {
          res.split <- splitcellpalm(
            X = X,
            valpts = intensity_tree[[i]]$nodeValpts,
            vecval = intensity_tree[[i]]$nodeCov,
            listcovariates = listcovariates,
            usecovariates = usedcov,
            clustmodel = "LGCP",
            areapixel = areapixel,
            threshold = threshold
          )
        } else {
          res.split <- splitcell(
            X = X,
            valpts = intensity_tree[[i]]$nodeValpts,
            vecval = intensity_tree[[i]]$nodeCov,
            usecovariates = usedcov,
            areapixel = areapixel,
            score = score,
            threshold = threshold
          )
        }
      }

      if (is.character(res.split)) {
        intensity_tree[[i]]$status <- 0
        intensity_tree[[i]]$already_split <- TRUE
        intensity_tree[[i]]$whystop <- res.split
        intensity_tree[[i]]$scrsplit <- NA
      } else {
        # Update the parent
        intensity_tree[[i]]$left_daughter <- knew + 1
        intensity_tree[[i]]$right_daughter <- knew + 2
        intensity_tree[[i]]$split_var <- res.split$split_var
        intensity_tree[[i]]$split_val <- res.split$split_val
        intensity_tree[[i]]$already_split <- TRUE
        intensity_tree[[i]]$scrsplit <- res.split$scrsplit
        intensity_tree[[i]]$scrdcr <- res.split$scrdcr

        # Define the children
        areasub <- (areapixel * sum(!is.na(res.split$sublevels[[res.split$split_var]])))
        areasup <- (areapixel * sum(!is.na(res.split$suplevels[[res.split$split_var]])))

        childleft <- list(
          nodeID = knew + 1,
          nodeCov = res.split$sublevels,
          nodeValpts = res.split$valptssub,
          nX = res.split$nsub,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = res.split$nsub / areasub,
          already_split = FALSE,
          whystop = NULL,
          scrsplit = NA,
          scrdcr = NA
        )

        childright <- list(
          nodeID = knew + 2,
          nodeCov = res.split$suplevels,
          nodeValpts = res.split$valptssup,
          nX = res.split$nsup,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = res.split$nsup / areasup,
          already_split = FALSE,
          whystop = NULL,
          scrsplit = NA,
          scrdcr = NA
        )
        # append the children
        intensity_tree <- append(
          intensity_tree,
          list(childleft, childright)
        )
      }
      knew <- length(intensity_tree)
    }
  }

  idterm <- sapply(intensity_tree, FUN = function(i) {
    i$status
  })

  patchworks <- lapply(intensity_tree[idterm == 0],
    FUN = function(i) {
      patchs <- ifelse(!is.na(i$nodeCov[[1]]),
        i$intensity_pred,
        0
      )
      return(im(
        matrix(patchs,
          nrow = dimcov[1],
          ncol = dimcov[2], byrow = F
        ),
        xrange = covrangex,
        yrange = covrangey
      ))
    }
  )

  for (i in seq_along(intensity_tree)) {
    intensity_tree[[i]]$nodeCov <- NULL
    intensity_tree[[i]]$nodeValpts <- NULL
  }

  if (inforest) {
    output <- list(
      tree = intensity_tree,
      X = NULL,
      namecov = names(listcovariates),
      im = spatstat.geom::as.im(Reduce("+", patchworks), W = X$window)
    )
  } else {
    output <- list(
      tree = intensity_tree,
      X = X,
      namecov = names(listcovariates),
      listcov = listcovariates,
      im = spatstat.geom::as.im(Reduce("+", patchworks), W = X$window)
    )
  }

  class(output) <- "sptree" # For when I will define class

  return(output)
}
