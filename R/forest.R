#' Random Intensity Forest
#'
#' @param X A spatial point process
#' \code{\link[spatstat.geom]{ppp}} object from spatstat.
#' @param listcovariates A list of covariates as
#' \code{\link[spatstat.geom]{im}} objects from spatstat.
#' @param Ntree A positive integer. 
#' The number of trees in the random intensity forest.
#' @param minpts Numeric. The minimum number of points in a region to allow a split.
#' @param mtry A number in \eqn{[0,1)}.
#' Probability that a covariate is used at each a split.
#' @param p A number in \eqn{[0,1)}.
#' Control the thinning process applied to the original point pattern X before
#' fitting a tree intensity estimator.
#' @param cores_trees A positive integer.
#' The number of cores used to computes the intensity trees.
#' @param score String specifying the score used to choose among splits, see details.
#' @param threshold A positive number.
#' The minimum area of a region for which we allow at most one split.
#' @param tol unused
#' @param minsplitq unused
#' @param maxsplitq unused
#'
#' @details
#' This function is the core of the package and compute random intensity
#' trees on a point pattern in presence of covariates.
#'
#' @return An object of class \code{spforest}.
#' @export
#'
#' @examples
#' forest <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = list(
#'     grad = spatstat.data::bei.extra$grad,
#'     elev = spatstat.data::bei.extra$elev
#'   ),
#'   Ntree = 3,
#'   minpts = 100,
#'   mtry = 1 / 3,
#'   p = 0,
#'   cores_trees = 1
#' )
RforestPP <- function(X,
                      listcovariates = NULL,
                      Ntree = 10,
                      minpts = spatstat.geom::npoints(X) / 10,
                      mtry = 1 / 3,
                      p = 0,
                      cores_trees = 1,
                      score = "lcv",
                      threshold = smallest_pixelarea(listcovariates),
                      tol = Inf,
                      minsplitq = 0.5,
                      maxsplitq = 0.5) {
  nbcov <- length(listcovariates)
  namescov <- names(listcovariates)

  if (!do.call(spatstat.geom::compatible.im, unname(listcovariates))) {
    listcovariates <- do.call(
      spatstat.geom::harmonise.im,
      listcovariates
    )
    warning("The im objects in listcovariates have been
    harmonised with the function harmonise.im.")
  }

  # Compute the forest's trees
  treeinforest <- parallel::mclapply(1:Ntree, FUN = function(i) {
    # Determine points in and out

    #########################
    # TODO: ALLOW FOR NO POINTS in the random selection OR NOT !!!!!! ??????
    ##################################

    if (p == 0) { # bootstrap case, with replacement
      ptintree <- sample.int(n = X$n, size = X$n, replace = T)
      Xintree <- X[ptintree]
    } else {
      ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)
      safety <- 1
      while (sum(ptintree) == 0 & safety <= 1e3) {
        ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)
        safety <- safety + 1
      }
      if (safety > 1e3) {
        stop("Check your parameters, there is no point sampled in the trees")
      }
      Xintree <- X[ptintree == 1]
      # ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)

      # Xintree <- X[ptintree == 1]
      # Xout <- X[ptintree != 1]
    }

    tree <- treerec(
      X = Xintree,
      threshold = threshold,
      score = score,
      listcovariates = listcovariates,
      mtry = mtry,
      tol = tol,
      minpts = minpts,
      minsplitq = minsplitq,
      maxsplitq = maxsplitq,
      inforest = T
    )

    # TODO: change printing methods to take tree from a forest into account
    # To save a lot of space in memory
    # emptylist <- vector("list", nbcov)
    # names(emptylist) <- namescov
    # tree$listcov <- emptylist

    return(list(
      sptree = tree,
      pt_intree = ptintree
    ))
  }, mc.cores = cores_trees)


  output <- list(
    trees = lapply(treeinforest,
      FUN = function(i) {
        i$sptree
      }
    ),
    pt_intree = lapply(treeinforest,
      FUN = function(i) {
        i$pt_intree
      }
    ),
    X = X,
    listcov = listcovariates,
    p = p,
    mtry = mtry
  )

  class(output) <- "spforest"

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

  output <- spatstat.geom::plot.im(as.im(x), main = main, ...)

  return(invisible(output))
}


#' Convert to Pixel Image
#'
#' @param X A spforest object
#' @param ...  ignored
#'
#' @return
#' @export
#'
#' @examples
as.im.spforest <- function(X, ...) {
  list_im <- lapply(X$trees, FUN = function(i) {
    i$im
  })

  if (X$p == 0) {
    output <- Reduce("+", list_im) / length(X$trees)
  } else {
    output <- Reduce("+", list_im) / length(X$trees) / x$p
  }

  return(output)
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

  if (x$p == 0) {
    output <- Reduce("+", list_im) / length(x$trees)
  } else {
    output <- Reduce("+", list_im) / length(x$trees) / x$p
  }

  return(output)
}

# `[.spforest` <- function(x, i, ...) {
#    as.im.spforest(x)[i]
#   return(output)
# }


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
importance <- function(forest, id_cov, cores = 1) {
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




#' OOB forest
#'
#' @param forest spforest
#' @param cores to speed up
#'
#' @return
#' @export
#'
#' @examples
OOBscr <- function(forest, cores = 1) {
  X <- forest$X # this is always the root

  # Put listcov back in the sptree object, required in predict.sptree
  for (i in 1:length(forest$trees)) {
    forest$trees[[i]]$listcov <- forest$listcov
  }

  OOBscr <- parallel::mclapply(1:length(forest$trees), FUN = function(i) {
    OOBval <- rep(NA, X$n)

    if (forest$p == 0) {
      torm <- unique(forest$pt_intree[[i]])
      if (length(torm) == X$n) { # If all points are drawn in the bootstrap,
        # then nothing do do.
        return(OOBval)
      }
      OOBpts <- 1:X$n
      OOBpts <- OOBpts[!(OOBpts %in% torm)]
    } else {
      # vector of same length as number of pts in X
      OOBpts <- (forest$pt_intree[[i]] != 1)
      if (all(!OOBpts)) { # If no points in OOBpts, then nothing do do.
        return(OOBval)
      }
    }

    # OOB sample
    Xout <- X[OOBpts]

    ### OOB prediction
    pts_pred_OOB <- predict.sptree(
      object = forest$trees[[i]],
      newdata = Xout
    )

    OOBval[OOBpts] <- pts_pred_OOB


    # OOB score
    return(OOBval)
  }, mc.cores = cores)

  logterm <- log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE))
  if (all(is.na(logterm))) {
    output <- NA
  } else {
    output <- sum(logterm, na.rm = TRUE)
  }
  # output <- sum(log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE)), na.rm = TRUE)

  # Return the average error of all the trees
  return(output)
}


#' OOB ppm
#'
#' @param forest  spforest
#' @param covariates list of covariates considered
#' @param cores  to speed up
#'
#' @return
#' @export
#'
#' @examples
OOBppmscr <- function(forest, covariates, cores = 1) {
  X <- forest$X # this is always the root

  OOBppm <- parallel::mclapply(1:length(forest$trees), FUN = function(i) {
    # vector of same length as number of pts in X
    OOBpts <- (forest$pt_intree[[i]] != 1)
    OOBval <- rep(NA, X$n)

    # OOB sample
    Xout <- X[OOBpts]
    Xin <- X[!OOBpts]

    modfit <- spatstat.model::ppm(Xin, ~., data = covariates)

    OOBval[OOBpts] <- spatstat.model::predict.ppm(modfit, locations = Xout)

    # OOB score
    return(OOBval)
  }, mc.cores = cores)

  # Return the average error of all the trees
  output <- sum(log(rowMeans(do.call(cbind, OOBppm), na.rm = TRUE)), na.rm = TRUE)
  return(output)
}
