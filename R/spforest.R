#' Random intensity forest
#'
#' @param X A spatial point process as a \code{\link[spatstat.geom]{ppp}} object.
#' @param listcovariates A list of covariates as \code{\link[spatstat.geom]{im}} objects.
#' @param Ntree Number of trees in the forest.
#' @param minpts A positive integer.
#' The minimum number of points allowed to try to split a cell.
#' @param mtry Probability of choosing a covariate.
#' @param p A number in \eqn{[0,1)}.
#' Control the thinning process applied to the original point pattern __X__ before
#' fitting a tree intensity estimator.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param threshold Minimum threshold to allow to split cell.
#' @param lambda An integer. The number of points used for random tessellation.
#' @param dimyx A vector of two integers. The dimensions of the output image passed
#' to \code{\link[spatstat.geom]{as.im}}.
#' @param test.connected Logical. If \code{TRUE},
#' \code{\link[spatstat.geom]{connected}} is applied to the tessellation to split tiles
#' in different connected components. It is only useful if the windows of
#' observation of \code{X} is not convex.
#' @param cores A positive integer.
#' The number of cores used to compute each intensity tree.
#'
#' @details 
#' This function is a wrapper around \code{\link{tesscovforest}} 
#' and \code{\link{tessforest}} to compute random intensity forest 
#' in the presence of covariates, 
#' or based on independent random tesselæation, respectively. 
#' If \code{listcovariates} is \code{NULL}, 
#' the function will call \code{\link{tessforest}}. Otherwise, it will call
#' \code{\link{tesscovforest}}.
#'
#' @return An object of class \code{\link{spforest}}.
#' @export
#'
#' @examples
#' forestwithcov <- spforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 200,
#'   mtry = 1 / 3,
#'   p = 0,
#'   cores = 1
#' )
#' plot(forestwithcov)
#' forestwithtoutcov <- tessforest(
#'   X = bei,
#'   lambda = 50,
#'   dimyx = c(101, 201),
#'   test.connected = FALSE,
#'   Ntree = 5,
#'   cores = 1
#' )
#' plot(forestwithtoutcov)
spforest <- function(X,
                     listcovariates = NULL,
                     Ntree = 10,
                     minpts = spatstat.geom::npoints(X) / 10,
                     mtry = 1 / 3,
                     p = 0,
                     score = "lcv",
                     threshold = smallest_pixelarea(listcovariates),
                     lambda = 100,
                     dimyx = c(50, 50),
                     test.connected = FALSE,
                     cores = 1) {
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
#' An object of this class represents a spatial intensity forest.
#' It contains 
#' \itemize{
#' \item \code{imforest}, an pixel image as 
#' \code{\link[spatstat.geom]{im.object}} object representing the value of 
#' the estimated intensity on the window of \code{X}.
#' \item \code{trees}, a list of the spatial intensity trees as \code{\link{sptree.object}}.
#' \item \code{ntrees}, the number of trees in the forest.
#' \item \code{pt_intree}, a list where the \eqn{i}-th element 
#' is a vector containing the index of the points of \code{X} used 
#' to compute the \eqn{i}-th intensity tree.
#' \item \code{X}, the original point pattern as a \code{\link[spatstat.geom]{ppp}} object.
#' \item \code{listcov}, a list of the covariates used to compute the forest
#' as \code{\link[spatstat.geom]{im}} objects. If \code{NULL}, 
#' independently random tessellation have been used.
#' \item \code{p}, a number in \eqn{[0,1)} 
#' controlling the thinning process applied to \code{X} 
#' before computing a tree intensity estimator.
#' \item \code{mtry}, probability that a covariate is used at each a split.
#' }
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
    "trees, of a point pattern with",
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
#' @param object A spatial intensity forest return by spforest function
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @return A vector of predicted intensity
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
#' predict(forest, c(100, 100))
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


#' Convert to Pixel Image
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
#' @param x the forest
#' @param cores To compute faster
#' @param ... ignored.
#'
#' @return Boxplot of the variable importances.
#' @importFrom graphics boxplot
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
boxplot.spforest <- function(x, cores = 1, ...) {
  vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
    importance(x, id_cov = i, cores = cores)
  })

  df <- data.frame(
    vipval = as.vector(vipval),
    var = rep(names(x$listcov), each = length(x$trees))
  )

  graphics::boxplot(vipval ~ var,
    data = df,
    xlab = "Variables",
    ylab = "Variable Importance"
  )
}


#' Vip barplot of forest
#'
#' @param x  An spforest object
#' @param sorted Logical. If TRUE, the variables are sorted by importance.
#' @param cores To compute faster
#' @param ... ignoted
#'
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
vipplot <- function(x, sorted = FALSE, cores = 1, ...) {
  vipval <- lapply(X = seq_along(x$listcov), FUN = function(i) {
    importance(x, id_cov = i, cores = cores)
  })

  avvip <- unlist(lapply(vipval, mean))

  avvipsort <- sort(avvip,
    decreasing = TRUE,
    index.return = TRUE
  )

  if (sorted) {
    graphics::barplot(avvipsort$x,
      names.arg = names(x$listcov)[avvipsort$ix],
      xlab = "Variables",
      ylab = "Variable Importance"
    )
    return(invisible(data.frame(
      Variable = names(x$listcov)[avvipsort$ix],
      Importance = avvipsort$x
    )))
  } else {
    graphics::barplot(avvip,
      names.arg = names(x$listcov),
      xlab = "Variables",
      ylab = "Variable Importance"
    )
    return(invisible(data.frame(
      Variable = names(x$listcov),
      Importance = avvip
    )))
  }
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
#' @param x First forest
#' @param y Second forest
#' @param ... ignored
#'
#' @return An \code{\link{spforest}} object
#' the intensity trees from both forest.
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
#'   Ntree = 3,
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
    trees = c(x$trees, y$trees),
    pt_intree = c(x$pt_intree, y$pt_intree),
    X = x$X,
    listcov = x$listcov,
    p = x$p,
    mtry = x$mtry
  )

  return(structure(output, class = "spforest"))
}



# Function below to remove ?? ----
#'
#' #' Constructor for spforest
#' new_spforest <- function(trees = list(),
#'                          pt_intree = list(),
#'                          X = ppp(),
#'                          listcov = list(),
#'                          p = double(),
#'                          mtry = double()) {
#'   output <- list(
#'     trees = trees,
#'     pt_intree = pt_intree,
#'     X = X,
#'     listcov = listcov,
#'     p = p,
#'     mtry = mtry
#'   )
#'
#'   return(structure(output, class = "spforest"))
#' }
#'
#'
#'
#'
#' #' Validator for spforest
#' #'
#' #' @param x
#' #'
#' #' @return The input forest if it is valid
#' #' @export
#' #'
#' #' @examples
#' #' forest <- spforest(
#' #'   X = spatstat.data::bei,
#' #'   listcovariates = spatstat.data::bei.extra,
#' #'   Ntree = 3,
#' #'   minpts = 300,
#' #'   mtry = 1
#' #' )
#' #' validate_spforest(forest)
#' validate_spforest <- function(x) {
#'   # TODO: add conformity checks on the entries of the forest
#'   values <- unclass(x)
#'
#'   # Check for correct entries' name in the forest
#'   if (sum(is.element(
#'     c(
#'       "trees", "pt_intree", "X",
#'       "listcov", "p", "mtry"
#'     ),
#'     names(x)
#'   )) != 6) {
#'     stop(
#'       "A spforest object should have entries
#'             trees, pt_intree, X, listcov, p, mtry",
#'       call. = FALSE
#'     )
#'   }
#'   # Check for correct type in the forest's entries
#'   if (!is.list(values$trees) | !is.list(values$listcov)) {
#'     stop("A spforest object must have entries trees and listcov as list.",
#'       call. = FALSE
#'     )
#'   }
#'   if (!is.ppp(values$X) | !is.list(values$listcov)) {
#'     stop("A spforest object must have entry X as a ppp object.",
#'       call. = FALSE
#'     )
#'   }
#'   if (!is.numeric(values$p) | !is.numeric(values$mtry)) {
#'     stop("A spforest object must have entries p and mtry as numeric.",
#'       call. = FALSE
#'     )
#'   }
#'
#'   # check if entries have possible values
#'   if (values$p < 0 | values$p > 1 | values$mtry > 1 | values$mtry < 0) {
#'     stop("p and mtry must be between 0 and 1.",
#'       call. = FALSE
#'     )
#'   }
#'
#'   alllength_ptintree <- sapply(values$pt_intree, length)
#'   if (!all(alllength_ptintree == npoints(values$X))) {
#'     stop("pt_intree should be a vector with
#'              same length as the number of points in X",
#'       call. = FALSE
#'     )
#'   }
#'
#'   x
#' }
