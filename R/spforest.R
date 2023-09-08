#' Printing spatial intensity forest
#'
#' @param x A spatial intensity forst return by RforestPP function
#' @param ... Additional arguments
#'
#' @return
#' @export
#'
#' @examples
print.spforest <- function(x, ...) {
  ncov <- length(x$listcov)

  cat(paste(
    "Intensity forest estimate of point patterns with",
    x$X$n, "points.\n\n"
  ))

  cat(paste(ncov, "covariables used, with names: "))
  cat(names(x$listcov), "\n")

  cat(
    "Spatial intensity forest with", length(x$trees), "trees."
  )
}


#' Compute image of a spatial intensity forest
#'
#' @param x A spatial intensity tree return by RforestPP function
#' @param ... additional arguments
#'
#' @return
#' @export
#'
#' @examples
imspforest <- function(x, ...) {
  list_im <- lapply(x$trees, FUN = function(i) {
    i$im
  })

  output <- Reduce("+", list_im) / length(x$trees) / x$p

  return(output)
}


#' Plot spatial intensity forest
#'
#' @param x A spatial intensity tree return by RforestPP function
#' @param ... additional arguments
#'
#' @return
#' @export
#'
#' @examples
plot.spforest <- function(x, ..., main) {
  # Handling case if no main title is given for the plot
  if (missing(main)) {
    main <- "Spatial Intensity Forest"
  }

  list_im <- lapply(x$trees, FUN = function(i) {
    i$im
  })

  if (x$p == 0) {
    output <- Reduce("+", list_im) / length(x$trees)
  } else {
    output <- Reduce("+", list_im) / length(x$trees) / x$p
  }

  spatstat.geom::plot.im(output, main = main, ...)

  return(invisible(output))
}




#' Importance of one covariable
#'
#' @param forest A spforest object
#' @param id_cov The id in forest$listcov of the covariable looked at.
#' @param cores how many cores to use to speed up computation
#'
#' @return
#' @export
#'
#' @examples
importance.spforest <- function(forest, id_cov, cores = 1) {
  # listZ <- forest$listcovsp
  X <- forest$X # this is alway the root
  Z <- forest$listcov[[id_cov]] # list of cov

  vip_tree <- NULL

  Zfun <- lapply(forest$listcov, spatstat.geom::as.function.im) # to remove?

  vip_tree <- parallel::mclapply(1:length(forest$trees), FUN = function(i) {
    # Shuffle the chosen covariate value
    dimmat <- Z$dim
    Z$v <- matrix(Z$v[sample.int(length(Z$v))],
      ncol = dimmat[2]
    )

    # OOB sample
    Xout <- X[forest$pt_intree[[i]] != 1]

    listZ_shuf <- forest$listcov
    listZ_shuf[[id_cov]] <- Z

    # Shuffle the tree
    treepert <- forest$trees[[i]]
    treepert$listcov <- listZ_shuf

    ### OOB prediction for shuffled covariable
    pts_pred_OOB_pert <- predict.sptree(
      object = treepert,
      newdata = Xout
    )

    ### OOB prediction
    # tree when in a forest do not have listcov. I add it
    forest$trees[[i]]$listcov <- forest$listcov
    pts_pred_OOB <- predict.sptree(
      object = forest$trees[[i]],
      newdata = Xout
    )

    # OOB mean square error of the tree
    return(sqrt(mean((pts_pred_OOB_pert - pts_pred_OOB)^2,
      na.rm = TRUE
    )))
  }, mc.cores = cores)


  # Return the error of all the trees
  return(unlist(vip_tree))
}




#' Boxplot forest
#'
#' @param x the forest
#' @param cores To compute faster
#' @param ... why note
#'
#' @return
#' @export
#'
#' @examples
boxplot_spforest <- function(x, cores = 1, ...) {
  vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
    importance.spforest(x, id_cov = i, cores = cores)
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
#' @param ... why not
#'
#' @return
#' @export
#'
#' @examples
vipplot.spforest <- function(x, cores = 1, ...) {
  vipval <- lapply(X = seq_along(x$listcov), FUN = function(i) {
    importance.spforest(x, id_cov = i, cores = cores)
  })

  avvip <- unlist(lapply(vipval, mean))

  graphics::barplot(avvip,
    names.arg = names(x$listcov),
    xlab = "Variables",
    ylab = "Variable Importance"
  )
}
# vipplot.spforest_old <- function(x, cores = 1, ...) {
#   vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
#     importance.spforest(x, id_cov = i, cores = cores)
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




#' Constructor for spforest
#'
#' @param trees
#' @param pt_intree
#' @param X
#' @param listcov
#' @param p
#' @param mtry
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
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


#' Merge two forests with same parameters
#'
#' @param x
#' @param y
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
