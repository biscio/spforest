#' Find the small pixel area in list of images
#'
#' @param x an im object
#'
#' @return A number.
#' @export
#'
#' @examples
#' smallest_pixelarea(spatstat.data::bei.extra)
smallest_pixelarea <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }

  allarea <- sapply(x, FUN = function(i) {
    y <- unclass(i)[c("xstep", "ystep")]
    pixelarea <- y$xstep * y$ystep
  })

  return(min(allarea))
}


#' Randomise input for tree
#'
#' @param listcovariates list of covariates used
#' @param mtry Probability of choosing a covariate.
#'
#' @return a vector of 0 and 1
#' @export
#'
#' @examples
#' rand_covar(
#'   listcovariates = list(1, 1, 1, 1, 1),
#'   mtry = 1 / 2
#' )
rand_covar <- function(listcovariates, mtry = 1) {
  nbcov <- 0

  while (nbcov == 0) {
    usedcov <- sample(c(0, 1),
      size = length(listcovariates),
      replace = T, prob = c(1 - mtry, mtry)
    )
    nbcov <- sum(usedcov)
  }

  return(usedcov)
}



#' Importance of one covariate
#'
#' @param forest A spforest object
#' @param id_cov The position in the list of covariates of \code{forest}.
#' @param cores How many cores to use to speed up computation
#'
#' @return A vector of the importance of the covariate \code{id_cov}
#' in the list of covariates of \code{forest}, for each tree.
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
#' importance(forest, id_cov = 1)
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

    ### OOB prediction for shuffled covariate
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
#' @return A number.
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
#' OOBscr(forest, cores = 1)
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
#' @return A number
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
# OOBppmscr(forest, covariates = spatstat.data::bei.extra)
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




# Function that to a nodeID find the parent ID. (Deprecated ??)
# Require the vectors of left_daughter and right_daughter
# from the tree
findparent <- function(ID, idleft, idright) {
  sapply(ID, FUN = function(i) {
    if (i %% 2 == 0) {
      parent <- which(i == idleft)
    } else {
      parent <- which(i == idright)
    }

    if (length(parent) != 1) {
      return(NULL)
    } else {
      return(parent)
    }
  })
}
