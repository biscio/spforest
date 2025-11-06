#' Intensity tessellation forest without covariate
#'
#' @inheritParams tesstree
#' @param Ntree A positive integer.
#' The number of trees in the random intensity forest.
#' @param verbose If TRUE, display progress bar.
#' @param parallel Check if future::plan has been called.
#'
#' @details
#' This function compute a random intensity forest
#' using \code{Ntree} independent spatial Intensity tree
#' computed with the function \code{\link{tesstree}}.
#'
#' The arguments \code{X}, \code{gamma}, \code{dimyx},
#' and \code{test.connected} are passed to \code{\link{tesstree}}.
#'
#' @return A pixel image, object of class \code{\link[spatstat.geom]{im.object}}.
#' @export
#'
#' @examples
#' Z <- tessforest(
#'   X = bei,
#'   gamma = 100,
#'   dimyx = c(101, 201),
#'   test.connected = FALSE,
#'   Ntree = 5
#' )
#' plot(Z)
tessforest <- function(X,
                       Ntree = 1,
                       gamma = 100,
                       dimyx = c(50, 50),
                       test.connected = FALSE,
                       verbose = FALSE,
                       parallel = TRUE) {
  if (is.null(gamma)) {
    gamma <- gamma_choice(X)
  }

  
  
  if (parallel & is(future::plan(), "sequential")) {
    stop('"parallel = T" but no parallel backend has been set with future::plan. To enable parallelism, we rely on the package future which allows for a fine control of the parallism backend: https://future.futureverse.org/index.html. Adding "future::plan("multisession", workers = N)" before the function, where N is the desired number of cores to use, will work in most cases.')
  }

  if (!parallel & !is(future::plan(), "sequential")) {
    message('"parallel=F" but a parallel backend with future::plan has been detected. No parallelism will be applied. Consider "parallel = T” to enable parallelism.')
    if (verbose) {
      progressr::handlers(global = TRUE)
      forestpgr <- function(x) {
        p <- progressr::progressor(along = x)
        lapply(x, FUN = function(i) {
          output <- tesstree(
            X = X,
            gamma = gamma,
            dimyx = dimyx,
            test.connected = test.connected
          )
          p(sprintf("x=%g", x))
          output
        })
      }
      listtree <- forestpgr(x = 1:Ntree)
      progressr::handlers(global = FALSE)
    } else {
      listtree <- lapply(1:Ntree, FUN = function(i) {
        tesstree(
          X = X,
          gamma = gamma,
          dimyx = dimyx,
          test.connected = test.connected
        )
      })
    }
  }

  if (parallel | is(future::plan(), "sequential")) {
    seedparallel <- TRUE
    if (is(future::plan(), "sequential")) {
      seedparallel <- NULL
    }

    if (verbose) {
      progressr::handlers(global = TRUE)
      forestpgr <- function(x) {
        p <- progressr::progressor(along = x)
        future.apply::future_lapply(x, FUN = function(i) {
          output <- tesstree(
            X = X,
            gamma = gamma,
            dimyx = dimyx,
            test.connected = test.connected
          )
          p(sprintf("x=%g", x))
          output
        }, future.seed = seedparallel)
      }
      listtree <- forestpgr(x = 1:Ntree)
      progressr::handlers(global = FALSE)
    } else {
      listtree <- future.apply::future_lapply(1:Ntree, FUN = function(i) {
        tesstree(
          X = X,
          gamma = gamma,
          dimyx = dimyx,
          test.connected = test.connected
        )
      }, future.seed = seedparallel)
    }
  }
  listim <- lapply(listtree, FUN = function(i) i$intensityim)
  listtess <- lapply(listtree, FUN = function(i) i$intensitytess)


  output <- list(
    imforest = Reduce("+", listim) / length(listim),
    trees = NULL,
    ntrees = length(listim),
    pt_intree = rep(1, X$n),
    X = X,
    listcov = NULL,
    p = NULL,
    mtry = NULL,
    tesselations = listtess
  )
  class(output) <- "spforest"

  return(output)
}


#' Intensity tessellation forest with covariates
#'
#' Compute a random intensity forest from an observed point pattern
#' and a list of covariates.
#' @param X A spatial point process as a
#' \code{\link[spatstat.geom]{ppp}} object from spatstat.
#' @param listcovariates A list of covariates as
#' \code{\link[spatstat.geom]{im}} objects from spatstat.
#' @param Ntree A positive integer.
#' The number of trees in the random intensity forest.
#' @param minpts A positive integer.
#' The minimum number of points in a region to allow a split.
#' @param mtry A number in \eqn{[0,1)}.
#' Probability that a covariate is used at each a split.
#' @param randmtry Logical. If \code{TRUE}, \code{mtry} must be between 0 and 1 and
#' represents the probability to use each covariate at each split. If \code{FALSE}, \code{mtry}
#' covariates are randomly chosen at each split.
#' @param p A number in \eqn{[0,1)}.
#' Control the thinning process applied to the original point pattern __X__ before
#' fitting a tree intensity estimator.
#' @param score String specifying the score used to choose among splits, see details.
#' @param threshold A positive number.
#' The minimum area of a region for which we allow at most one split.
#' @param verbose If TRUE, display progress bar.
#' @param parallel Check if future::plan has been called.
#'
#' @details
#' This function compute a random intensity forest using the covariates given
#' in \code{listcovariates}.
#' First the points of \code{X} are thinned with probability
#' \code{p} if \eqn{p>0}, or bootstrapped if \eqn{p=0}.
#' Then the function call \code{\link{tesscovtree}} to compute
#' an intensity tree estimate. This is repeated \code{Ntree} times.
#' The result returns in an \code{\link{spforest}} object.
#'
#' When computing a random intensity, one need a stopping criterion after
#' which we do not try to split the domain of observation any more.
#' Two stopping criterion are implemented: based on the minimal area of the
#' window and on the minimal number of points.
#' When there is less than \code{minpts} in a cell,
#' we do not try to split the cell any more.
#' When a cell has an area below \code{threshold} we do not attempt to split.
#' In principle, this criterion can be used jointly but we have empirically
#' found that it was more convenient and similar in term of performances
#' to consider only \code{minpts} and by default, \code{threshold}
#' is the smallest area of a pixel in the covariates in \code{listcovariates}.
#' Note that it does not make sense
#' to have \code{threshold} smaller than the smallest pixel area in the covariates.
#'
#' To increase the diversity between trees in random forest,
#' a common practice is to consider different covariates at each split.
#' This is controlled by \code{mtry}
#' which is the probability that a covariate, independently of the others,
#' is used at each split.
#'
#' We let \eqn{n_1} and \eqn{n_2} be the number of points in the regions \eqn{W1},
#' \eqn{W2}, respectively. We let |W| be the area of a region W.
#' Currently, several scores implemented to compute the performance of a split:
#' \itemize{
#' \item "lcv", \eqn{\frac{n_1 \log{n_1-1}}{|W_1|}\mathbf{1}_{n_1>1} +
#' \frac{n_2\log{n_2-1}}{|W_2|} \mathbf{1}_{n_2>1}}
#' \item "lcv2", \eqn{\frac{n_1\log{n_1-1}}{|W_1|} +
#' \frac{n_2\log{n_2-1}}{|W_2|} + \infty \mathbf{1}_{n_1=1} + \infty \mathbf{1}_{n_2=1}}
#' \item "ent", \eqn{\frac{n_1\log{n_1}}{|W_1|}\mathbf{1}_{n_1>1} +
#' \frac{n_2\log{n_2}}{|W_2|} \mathbf{1}_{n_2>1}}
#' \item "star", \eqn{|\frac{n_1}{n_1+n_2} - \frac{|W_1|}{(|W_1| + |W_2|)}|}
#' \item "ise", \eqn{-\frac{n_1^2}{|W_1|} - \frac{n_2^2}{|W_2|}}
#' \item "isecv", \eqn{-\frac{n_1^2-2n_1}{|W_1|} - \frac{n_2^2-2n_2}{|W_2|}}
#' }
#' Empirically, there has been no score outputting better
#' integrated mean squares error than others.
#'
#'
#' @return An object of class \code{\link{spforest}}.
#' @export
#'
#' @references \url{Article arxiv version}
#'
#' @examples
#' forest <- tesscovforest(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   Ntree = 3,
#'   minpts = 100,
#'   mtry = 1 / 3,
#'   p = 0
#' )
tesscovforest <- function(X,
                          listcovariates = NULL,
                          Ntree = 10,
                          minpts = spatstat.geom::npoints(X) / 10,
                          mtry = 1 / 3,
                          randmtry = FALSE,
                          p = 0,
                          score = "lcv",
                          threshold = smallest_pixelarea(listcovariates),
                          verbose = FALSE,
                          parallel = TRUE) {
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

  # Used in all trees
  areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
  vecval <- lapply(listcovariates, FUN = function(i) {
    if (!spatstat.geom::is.im(i)) {
      stop("Elements of listcovar must be an im object")
    }
    c(as.matrix.im(i))
  })
  dimcov <- listcovariates[[1]]$dim
  covrangex <- listcovariates[[1]]$xrange
  covrangey <- listcovariates[[1]]$yrange

  plantingtree <- function(i) {
    # Determine points in and out

    #########################
    # TODO: ALLOW FOR NO POINTS in the random selection OR NOT !!!!!! ??????
    ##################################

    if (p == 0) { # bootstrap case, with replacement
      ptintree <- sample.int(n = X$n, size = X$n, replace = T)
      Xintree <- X[ptintree] # Should I use unique(X[ptintree]) ??
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

    tree <- tesscovtree(
      X = Xintree,
      vecval = vecval,
      areapixel = areapixel,
      dimcov = dimcov,
      covrangex = covrangex,
      covrangey = covrangey,
      listcovariates = listcovariates,
      minpts = minpts,
      mtry = mtry,
      randmtry = randmtry,
      score = score,
      threshold = threshold,
      inforest = T
    )

    # summary(tree)[,"split_var"] |> unlist() |> table()

    # TODO: change printing methods to take tree from a forest into account
    # To save a lot of space in memory
    # emptylist <- vector("list", nbcov)
    # names(emptylist) <- namescov
    # tree$listcov <- emptylist

    return(list(
      sptree = tree,
      pt_intree = ptintree
    ))
  }

  if (parallel & is(future::plan(), "sequential")) {
    stop('"parallel = T" but no parallel backend has been set with future::plan. To enable parallelism, we rely on the package future which allows for a fine control of the parallism backend: https://future.futureverse.org/index.html. Adding "future::plan("multisession", workers = N)" before the function, where N is the desired number of cores to use, will work in most cases.')
  }

  if (parallel == FALSE & !is(future::plan(), "sequential")) {
    message('"parallel=F" but a parallel backend with future::plan has been detected. No parallelism will be applied. Consider "parallel = T” to enable parallelism.')
    if (verbose) {
      progressr::handlers(global = TRUE)
      forestpgr <- function(x) {
        p <- progressr::progressor(along = x)
        lapply(x, FUN = function(i) {
          output <- plantingtree(i)
          p(sprintf("x=%g", x))
          output
        })
      }
      treeinforest <- forestpgr(x = 1:Ntree)
      progressr::handlers(global = FALSE)
    } else {
      treeinforest <- lapply(1:Ntree,
        FUN = plantingtree
      )
    }
  }

  if (parallel | is(future::plan(), "sequential")) {
    # Compute the forest's trees - check if need to parallel
    seedparallel <- TRUE
    if (is(future::plan(), "sequential")) {
      seedparallel <- NULL
    }

    if (verbose) {
      progressr::handlers(global = TRUE)
      forestpgr <- function(x) {
        p <- progressr::progressor(along = x)
        future.apply::future_lapply(x, FUN = function(i) {
          output <- plantingtree(i)
          p(sprintf("x=%g", x))
          output
        }, future.seed = seedparallel)
      }
      treeinforest <- forestpgr(x = 1:Ntree)
      progressr::handlers(global = FALSE)
    } else {
      treeinforest <- future.apply::future_lapply(1:Ntree,
        FUN = plantingtree,
        future.seed = seedparallel
      )
    }
  }

  # Computation of the image
  list_im <- lapply(treeinforest, FUN = function(i) {
    i$sptree$im
  })
  if (p == 0) {
    imforest <- Reduce("+", list_im) / length(treeinforest)
  } else {
    imforest <- Reduce("+", list_im) / length(treeinforest) / p
  }

  output <- list(
    imforest = imforest,
    trees = lapply(treeinforest,
      FUN = function(i) {
        i$sptree
      }
    ),
    ntrees = length(list_im),
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
