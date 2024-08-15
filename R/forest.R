#' Title
#'
#' @param X XXX
#' @param listcovariates XXX
#' @param score XXX
#' @param p XXX
#' @param Ntree XXX
#' @param threshold XXX
#' @param cores_trees XXX
#' @param mtry XXX
#' @param tol XXX
#' @param minpts XXX
#' @param minsplitq XXX
#' @param maxsplitq XXX
#'
#' @return
#' @export
#'
#' @examples
#' forest <- RforestPP(
#'   X = spatstat.data::bei,
#'   listcovariates = list(
#'     grad = spatstat.data::bei.extra$grad,
#'     elev = spatstat.data::bei.extra$elev
#'   ),
#'   score = "lcv2",
#'   p = 1,
#'   Ntree = 3,
#'   threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
#'   cores_trees = 1,
#'   mtry = 1 / 3,
#'   tol = Inf,
#'   minpts = 50,
#'   minsplitq = 0.5,
#'   maxsplitq = 0.5
#' )
RforestPP <- function(X,
                      listcovariates = NULL,
                      score = "lcv",
                      p = 0.8,
                      Ntree = 10,
                      threshold = spatstat.geom::area(X) / 1e4,
                      cores_trees = 1,
                      mtry = 1 / 3,
                      tol = Inf,
                      minpts = spatstat.geom::npoints(X) / 10,
                      minsplitq = 0.5,
                      maxsplitq = 0.5) {
  nbcov <- length(listcovariates)
  namescov <- names(listcovariates)

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





#' OOB forest
#'
#' @param forest spforest
#' @param cores to speed up
#'
#' @return
#' @export
#'
#' @examples
OOBscr.spforest <- function(forest, cores = 1) {
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
