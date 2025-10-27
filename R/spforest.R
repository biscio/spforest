#' Random intensity forest
#'
#' Estimate the intensity of a spatial point process by random intensity forest.
#'
#' @param X A spatial point process as a \code{\link[spatstat.geom]{ppp}} object.
#' @param listcovariates A list of covariates as \code{\link[spatstat.geom]{im}} objects.
#' @param Ntree Number of trees in the forest.
#' @param minpts A positive integer.
#' The minimum number of points after which we try to split a cell one last time.
#' @param mtry Probability of choosing a covariate.
#' @param randmtry Logical. If \code{TRUE}, \code{mtry} must be between 0 and 1 and
#' represents the probability to use each covariate at each split. If \code{FALSE}, \code{mtry}
#' covariates are randomly chosen at each split.
#' @param p A number in \eqn{[0,1)}
#' controlling the thinning process applied to \code{X}.
#' See \code{\link[spforest]{tesscovforest}} for details.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param threshold Minimum area to allow to split cell.
#' @param lambda An integer. The number of points used for random tessellation.
#' @param dimyx A vector of two integers. The dimensions of the output image passed
#' to \code{\link[spatstat.geom]{as.im}}.
#' @param test.connected Logical. If \code{TRUE},
#' \code{\link[spatstat.geom]{connected}} is applied to the tessellation to split tiles
#' in different connected components. It is only useful if the windows of
#' observation of \code{X} is not convex.
#' @param cores A positive integer. The number of cores to use.
#' If strictly larger than 1, parallel computing
#' is used to dispatch each intensity tree computation on different cores.
#'
#' @details
#' If \code{listcovariates} is not \code{NULL},
#' the function computes a random intensity forest using
#' \code{\link{tesscovforest}}. Otherwise, it computes a random intensity forest
#' using \code{\link{tessforest}}.
#' All arguments are passed to the corresponding function.
#'
#' @return An object of class \code{\link{spforest}}.
#' @export
#'
#' @examples
#' forestwithcov <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 5,
#'   minpts = 200,
#'   mtry = 2
#' )
#' plot(forestwithcov)
#' forestwithtoutcov <- tessforest(
#'   X = bei,
#'   lambda = 50,
#'   dimyx = c(101, 201),
#'   Ntree = 100
#' )
#' plot(forestwithtoutcov)
spforest <- function(X,
                     listcovariates = NULL,
                     Ntree = 10,
                     minpts = spatstat.geom::npoints(X) / 10,
                     mtry = 1 / 3,
                     randmtry = FALSE,
                     p = 0,
                     score = "lcv",
                     threshold = 2 * smallest_pixelarea(listcovariates),
                     lambda = NULL,
                     dimyx = c(50, 50),
                     test.connected = FALSE,
                     cores = 1) {
  if (!is.null(X$mesh)) {
    # lambda <- spatstat.geom::npoints(X)^(2 / 3)
    if (is.null(lambda)) {
      lambda0 <- nrow(X$pp)^(2 / 3)
    } else {
      lambda0 <- lambda
    }

    triangle_density <- manifold_forest(
      Ntrees = Ntree,
      intensity = lambda0,
      mesh = X$mesh,
      pointsech = X$pp
    )

    output <- list(
      tridensity = triangle_density,
      mesh = X$mesh,
      pp = X$pp
    )

    class(output) <- "spforestmesh"

    return(output)
  }

  if (is.null(listcovariates)) {
    output <- tessforest(X,
      Ntree = Ntree,
      lambda = lambda,
      dimyx = dimyx,
      test.connected = test.connected,
      cores = cores
    )
  } else {
    output <- tesscovforest(X,
      listcovariates = listcovariates,
      Ntree = Ntree,
      minpts = minpts,
      mtry = mtry,
      randmtry = randmtry,
      p = p,
      cores = cores,
      score = score,
      threshold = threshold
    )
  }

  return(output)
}


#' spforest.object
#'
#' S3 class used to define random forest intensity
#'
#' @name spforest.object
#' @rdname spforest.object
#'
#' @details
#' This class represents a spatial intensity forest and includes information
#' about the original point pattern and how the random intensity forest has been computed.
#' We recommend to always generate an \code{spforest.object} object with the function
#' \code{\link[spforest]{spforest}}.
#' If \code{RFI} is a spatial intensity forest, it contains
#' \itemize{
#' \item \code{imforest}, a pixel image as
#' \code{\link[spatstat.geom]{im.object}} object representing the value of
#' the estimated intensity on the window of \code{X}.
#' \item \code{trees}, a list of the spatial intensity trees as \code{\link{sptree.object}}.
#' \item \code{ntrees}, the number of trees in the forest.
#' \item \code{pt_intree}, a list where the \eqn{i}-th element
#' is a vector containing the index of the points of \code{X} used
#' to compute the \eqn{i}-th intensity tree.
#' \item \code{X}, the original point pattern as a \code{\link[spatstat.geom]{ppp}} object.
#' \item \code{listcov}, the list of the covariates passed in the argument \code{listcovariates}
#' of \code{\link[spforest]{spforest}} as \code{\link[spatstat.geom]{im}} objects.
#'  If \code{NULL}, independent random tessellations have been used instead of
#'  a given list of covariates.
#' \item \code{p}, a number in \eqn{[0,1)}
#' passed in the argument \code{p}
#' of \code{\link[spforest]{spforest}}.
#' controlling the thinning process applied to \code{X}
#' before computing a tree intensity estimator.
#' \item \code{mtry}, the argument \code{mtry}
#' of \code{\link[spforest]{spforest}}.
#' }
#' The class \code{spforest.object} has methods
#' \code{\link[spforest]{print.spforest}},  \code{\link[spforest]{plot.spforest}}
#' and  \code{\link[spforest]{predict.spforest}}.
NULL

#' Printing spatial intensity forest
#'
#' @param x A spatial intensity forest return by spforest function
#' @param ... Additional arguments
#'
#' @return A description of the random intensity forest with
#' the number of points, number of covariates, and number of trees used.
#' @export
#'
#' @examples
#' forest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' print(forest)
print.spforest <- function(x, ...) {
  cat(paste(
    "Spatial intensity forest with", x$ntrees,
    "trees of a point pattern with",
    x$X$n, "points.\n\n"
  ))


  if (is.null(x$listcov)) {
    cat("No covariate has been given: each tree has been generated with a random tesselation.\n")
  } else {
    ncov <- length(x$listcov)
    a <- paste0(names(x$listcov), collapse = "", sep = ", ")
    namecov <- paste0(substr(a, 1, nchar(a) - 2), ".")

    cat(paste(ncov, "covariates used: "))
    cat(namecov, "\n")
  }
}

#' Forest prediction
#'
#' Given a random forest intensity estimator obtained by
#' \code{\link[spforest]{spforest}}, evaluate the intensity at new locations.
#'
#' @param object A spatial intensity forest \code{\link[spforest]{spforest.object}}.
#' @param newdata A matrix, a \code{c(x,y)} vector or
#' \code{\link[spatstat.geom]{ppp.object}} to be taken
#' as the new locations where the estimated intensity is evaluated.
#' @param ... Ignored.
#'
#' @details \code{predict.spforest} return the values of the
#' estimated intensity, obtained by evaluating the pixels values of \code{as.im(object)}
#' at locations given by \code{newdata}.
#'
#' @return A vector of numeric values corresponding to the estimated intensity
#' at the locations given by \code{object}. Currently \code{newdata} can be:
#' \itemize{
#' \item a vector of length 2 \code{c(x,y)},
#' representing the \eqn{x} and \eqn{y} coordinates of a single location.
#' \item a matrix with two columns, representing the \eqn{x} and \eqn{y}
#' coordinates of the new locations.
#' \item a \code{\link[spatstat.geom]{ppp.object}} representing
#' the coordinates of the new locations.
#' }
#' If \code{newdata} is missing or \code{NULL},
#' the intensity is evaluated at the points of \code{object$X}.
#' @export
#'
#' @examples
#' RFI <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 10,
#'   minpts = 200,
#'   mtry = 1
#' )
#' predict(RFI, c(100, 100))
#' newloc <- spatstat.random::runifpoint(n = 10, win = spatstat.data::bei$w)
#' predict(RFI, newloc)
predict.spforest <- function(object, newdata, ...) {
  # Handles the newdata to be in the correct form
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else if (!spatstat.geom::is.ppp(newdata)) {
    if (is.vector(newdata)) {
      X <- spatstat.geom::ppp(
        x = newdata[1],
        y = newdata[2],
        window = object$X$window
      )
    } else if (is.matrix(newdata)) {
      X <- spatstat.geom::ppp(
        x = newdata[, 1],
        y = newdata[, 2],
        window = object$X$window
      )
    }
  } else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n == 0) {
      return(NULL)
    }
    X <- newdata
  }

  return(as.im(object)[X])
}


#' Plot spatial intensity forest
#'
#' Plot a \code{\link{spforest.object}}.
#'
#' @param x A spatial intensity tree return by spforest function
#' @param ... additional arguments
#' @param main A title for the plot.
#'
#' @details
#' This function first convert an \code{\link{spforest}}
#' as an \code{\link[spatstat.geom]{im}} object and then plot it with
#' \code{\link[spatstat.geom]{plot.im}}. All arguments in \code{...} are passed to
#' \code{\link[spatstat.geom]{plot.im}}.
#' @return Same as \code{\link[spatstat.geom]{plot.im}}.
#' @export
#'
#' @examples
#' forest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' plot(forest)
plot.spforest <- function(x, ..., main = "Spatial Intensity Forest") {
  # Handling case if no main title is given for the plot
  # if (missing(main)) {
  #   main <- "Spatial Intensity Forest"
  # }

  output <- spatstat.geom::plot.im(as.im.spforest(x), main = main, ...)

  return(invisible(output))
}

#' Plot spatial intensity forest on a manifold
#'
#' Plot a \code{\link{spforest.object}}.
#'
#' @param x A spforestmesh object
#' @param ... additional arguments
#' @param main A title for the plot.
#'
#' @details details
#' @export
#'
#' @examples
plot.spforestmesh <- function(x, points = FALSE, colorbar = FALSE, log=F, offset = exp(-8), ...) {
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("The package RANN must be installed.")
  }

  if (class(x) != "spforestmesh") {
    stop("The object to plot must be of the class spforestmesh")
  }
  
  if (isTRUE(log)) {
    y <- log(x+offset)
  } else {
    y <- x
  }
  
  if (points) {
    output <- plot_manifold_intensity(y,
      points = TRUE,
      colorbar = colorbar,
      zoom = 0.65,
      nticks = 10,
      lasttick = TRUE
    )
  } else {
    output <- plot_manifold_intensity(y,
      points = FALSE,
      colorbar = colorbar,
      zoom = 0.65,
      nticks = 10,
      lasttick = TRUE
    )
  }

  return(invisible(output))
}



#' Convert to Pixel Image
#'
#' Wrapper function to extract imforest from an \code{\link{spforest.object}}.
#'
#' @param X A spforest object
#' @param ...  ignored
#'
#' @return A pixel image \code{\link[spatstat.geom]{im.object}}.
#' @import spatstat.geom
#' @import spatstat.model
#' @import spatstat.random
#' @export
#'
#' @examples
#' forest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' as.im(forest)
as.im.spforest <- function(X, ...) {
  # list_im <- lapply(X$trees, FUN = function(i) {
  #   i$im
  # })
  #
  # if (X$p == 0) {
  #   output <- Reduce("+", list_im) / length(X$trees)
  # } else {
  #   output <- Reduce("+", list_im) / length(X$trees) / X$p
  # }

  return(X$imforest)
}


#' Boxplot forest
#'
#' Plot the boxplot of the importance of each covariate
#' in each tree of the random forest intensity
#'
#' @usage boxplot(x = forest)
#'
#' @param x A \code{\link{spforest.object}} with \code{listcov} not \code{NULL}.
#' @param ... Arguments passed to \code{\link[graphics]{boxplot}}.
#'
#' @details
#' First the function \code{\link[spforest]{importance}} is called to
#' compute the importance of each covariate in each tree \code{x}.
#' Then, a boxplot is drawn to visualize the distribution of the
#' importance of each covariate.
#'
#' @return A list returned by the function \code{\link[graphics]{boxplot}}.
#' @param viptype An integer in \{1,2,3,4\} passed to \code{\link[spforest]{importance}}.
#' @export
#'
#' @examples
#' forest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' boxplot(forest)
boxplot.spforest <- function(x, viptype = 4, ...) {
  vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
    importance(x, id_cov = i, viptype = viptype)
  })

  df <- data.frame(
    vipval = as.vector(vipval),
    var = rep(names(x$listcov), each = length(x$trees))
  )

  graphics::boxplot(vipval ~ var,
    data = df,
    xlab = "Variables",
    ylab = "Variable Importance",
    ...
  )
}


#' Vip barplot of forest
#'
#' @param x  A \code{\link{spforest.object}} with \code{listcov} not \code{NULL}.
#' @param sorted Logical. If TRUE, the variables are sorted by importance.
#' @param viptype An integer in \{1,2,3,4\} passed to \code{\link[spforest]{importance}}.
#' @param ... Arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @details
#' First the function \code{\link[spforest]{importance}} is called to
#' compute the importance of each covariate in each tree \code{x}.
#' Then, a barplot is drawn to visualize the distribution
#' of the importance of each covariate, averaged over all the trees.
#' @return Variable importance plot
#' @export
#'
#' @examples
#' forest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' vipplot(forest, sorted = TRUE)
vipplot <- function(x,
                    sorted = FALSE,
                    viptype = 4, ...) {
  output <- vip(x, viptype = viptype, treesdetails = FALSE)

  if (sorted) {
    output <- sort(output,
      decreasing = TRUE
    )
  }

  graphics::barplot(output,
    names.arg = names(output),
    ...
  )

  return(invisible(output))
}
# vipplot.spforest_old <- function(x, cores = 1, ...) {
#   vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
#     importance(x, id_cov = i, cores = cores)
#   })
#
#   df <- data.frame(
#     vipval = as.vector(vipval),
#     var = rep(names(x$listcov), each = length(x$trees))
#   )
#
#   input <- by(
#     data = df,
#     INDICES = df$var,
#     FUN = function(y) {
#       mean(y$vipval)
#     }
#   )
#
#   graphics::barplot(input,
#     names.arg = unique(df$var),
#     xlab = "Variables",
#     ylab = "Variable Importance"
#   )
# }


#' Merge two forests with same parameters
#'
#' Given two random forest intensity estimates,
#' merge them into one \code{\link{spforest}[spforest.object]}.
#'
#' @param x,y Two \code{\link{spforest}[spforest.object]} with identical entries
#' \code{p}, \code{mtry}, \code{listcov} and \code{X}.
#' @param ... Ignored
#'
#' @details
#' The function merges two \code{\link{spforest}[spforest.object]} by concatenating
#' their entires \code{trees}, \code{pt_intree} and
#' averaging their intensity estimates given by their entries \code{imforest}.
#' If any of the entries \code{p}, \code{mtry}, \code{listcov} or \code{X} is different,
#' the function returns an error.
#' @return An \code{\link[spforest]{spforest.object}}
#' combining the information of \code{x} and \code{y}.
#' @export
#'
#' @examples
#' forest1 <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' forest2 <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 5,
#'   minpts = 300,
#'   mtry = 1
#' )
#' merge(forest1, forest2)
merge.spforest <- function(x, y, ...) {
  # Check if x and y are valid forests
  # TODO: how to handle if they are not?
  # validate_spforest(x)
  # validate_spforest(y)

  # TODO: test the forests are mergeable
  if (x$p != y$p) {
    stop("The parameter p must be the same in both forests.")
  }
  if (x$mtry != y$mtry) {
    stop("The parameter mtry must be the same in both forests.")
  }
  if (x$X$n != y$X$n) {
    stop("The point patterns do not appear to be the same in both forests.
             They have different number of points.")
  }

  output <- list(
    imforest = (x$imforest + y$imforest) / 2,
    trees = c(x$trees, y$trees),
    ntrees = x$ntrees + y$ntrees,
    pt_intree = c(x$pt_intree, y$pt_intree),
    X = x$X,
    listcov = x$listcov,
    p = x$p,
    mtry = x$mtry
  )

  return(structure(output, class = "spforest"))
}


#' Extract a smaller forest
#'
#' Given a random forest intensity estimate, return a random forest intensity estimated
#' with a only subset of all the trees.
#'
#' @param forest A \code{\link{spforest}[spforest.object]} .
#' @param whichtrees A vector indicating which trees to use in the
#' new forest.
#'
#' @details The function extracts the trees indicated in \code{whichtrees}
#' from \code{forest} and returns a new \code{\link{spforest}[spforest.object]}.
#' The intensity estimate of the new forest is obtained by averaging the intensity
#' estimates of the selected trees. All the other arguments are updated accordingly.
#' @return A \code{\link{spforest}[spforest.object]} .
#' @export
#'
#' @examples
#' bigforest <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' smallforest <- extractforest(
#'   forest = bigforest,
#'   whichtrees = c(2, 5)
#' )
extractforest <- function(forest, whichtrees) {
  trees1 <- forest$trees[whichtrees]

  imsubtrees <- lapply(trees1, FUN = function(i) {
    i$im
  })

  imforest1 <- Reduce("+", imsubtrees) / length(whichtrees)

  output <- list(
    imforest = imforest1,
    trees = trees1,
    ntrees = length(whichtrees),
    pt_intree = forest$pt_intree[whichtrees],
    X = forest$X,
    listcov = forest$listcov,
    p = forest$p,
    mtry = forest$mtry
  )
  return(structure(output, class = "spforest"))
}
