#' Printing spatial intensity forest
#'
#' @param x A spatial intensity forst return by RforestPP function
#' @param ... Additional arguments
#'
#' @return A description of the random intensity forest with
#' the number of points, number of covariates, and number of trees used.
#' @export
#'
#' @examples
#' forest <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' print(forest)
print.spforest <- function(x, ...) {
  ncov <- length(x$listcov)
  a <- paste0(names(x$listcov), collapse = "", sep = ", ")
  namecov <- paste0(substr(a, 1, nchar(a) - 2), ".")

  cat(paste(
    "Intensity forest estimate of point patterns with",
    x$X$n, "points.\n\n"
  ))

  cat(paste(ncov, "covariables used, with names: "))
  cat(namecov, "\n")

  cat(
    "Spatial intensity forest with", length(x$trees), "trees."
  )
}

#' Forest prediction
#'
#' @param object A spatial intensity forest return by RforestPP function
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @return A vector of predicted intensity
#' @export
#'
#' @examples
#' forest <- RforestPP(
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

  return(as.im(forest)[X])
}


#' Boxplot forest
#'
#' @param x the forest
#' @param cores To compute faster
#' @param ... ignored.
#'
#' @return Boxplot of the variable importances.
#' @export
#'
#' @examples
#' forest <- RforestPP(
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
#' @param cores To compute faster
#' @param ... ignoted
#'
#' @return Variable importance plot
#' @export
#'
#' @examples
#' forest <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' vipplot(forest, sorted = T)
vipplot <- function(x, sorted = F, cores = 1, ...) {
  vipval <- lapply(X = seq_along(x$listcov), FUN = function(i) {
    importance(x, id_cov = i, cores = cores)
  })

  avvip <- unlist(lapply(vipval, mean))

  avvipsort <- sort(avvip,
    decreasing = T,
    index.return = T
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
#' forest1 <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' forest2 <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' merge.spforest(forest1, forest2)
merge.spforest <- function(x, y, ...) {
  # Check if x and y are valid forests
  # TODO: how to handle if they are not?
  validate_spforest(x)
  validate_spforest(y)

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

  output <- new_spforest(c(x$trees, y$trees),
    pt_intree = c(x$pt_intree, y$pt_intree),
    X = x$X,
    listcov = x$listcov,
    p = x$p,
    mtry = x$mtry
  )

  return(validate_spforest(output))
}



# Function below to remove ?? ----

#' Constructor for spforest
new_spforest <- function(trees = list(),
                         pt_intree = list(),
                         X = ppp(),
                         listcov = list(),
                         p = double(),
                         mtry = double()) {
  output <- list(
    trees = trees,
    pt_intree = pt_intree,
    X = X,
    listcov = listcov,
    p = p,
    mtry = mtry
  )
  
  return(structure(output, class = "spforest"))
}




#' Validator for spforest
#'
#' @param x
#'
#' @return The input forest if it is valid
#' @export
#'
#' @examples
#' forest <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 300,
#'   mtry = 1
#' )
#' validate_spforest(forest)
validate_spforest <- function(x) {
  # TODO: add conformity checks on the entries of the forest
  values <- unclass(x)

  # Check for correct entries' name in the forest
  if (sum(is.element(
    c(
      "trees", "pt_intree", "X",
      "listcov", "p", "mtry"
    ),
    names(x)
  )) != 6) {
    stop(
      "A spforest object should have entries
            trees, pt_intree, X, listcov, p, mtry",
      call. = FALSE
    )
  }
  # Check for correct type in the forest's entries
  if (!is.list(values$trees) | !is.list(values$listcov)) {
    stop("A spforest object must have entries trees and listcov as list.",
      call. = FALSE
    )
  }
  if (!is.ppp(values$X) | !is.list(values$listcov)) {
    stop("A spforest object must have entry X as a ppp object.",
      call. = FALSE
    )
  }
  if (!is.numeric(values$p) | !is.numeric(values$mtry)) {
    stop("A spforest object must have entries p and mtry as numeric.",
      call. = FALSE
    )
  }

  # check if entries have possible values
  if (values$p < 0 | values$p > 1 | values$mtry > 1 | values$mtry < 0) {
    stop("p and mtry must be between 0 and 1.",
      call. = FALSE
    )
  }

  alllength_ptintree <- sapply(values$pt_intree, length)
  if (!all(alllength_ptintree == npoints(values$X))) {
    stop("pt_intree should be a vector with
             same length as the number of points in X",
      call. = FALSE
    )
  }

  x
}
