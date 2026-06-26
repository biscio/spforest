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



#' Tesselation intensityy by rule of thumb
#'
#' @param X The observed data point pattern,
#' as a \code{\link[spatstat.geom]{ppp.object}}.
#'
#' @return A number.
#' @export
#'
#' @examples
#' gamma_choice(spatstat.random::rpoispp(100))
gamma_choice <- function(X) {
  gamma0 <- floor(mean(c(
    grDevices::nclass.FD(X$x),
    grDevices::nclass.FD(X$y)
  ))^2)

  return(gamma0)
}


#' Randomise input for tree
#'
#' @param listcovariates list of covariates used
#' @param mtry Probability of choosing a covariate.
#' @param randmtry Logical. If \code{TRUE}, \code{mtry} must be between 0 and 1 and
#' represents the probability to use each covariate at each split. If \code{FALSE}, \code{mtry}
#' covariates are randomly chosen at each split.
#'
#' @return a vector of 0 and 1
#' @export
#'
#' @examples
#' rand_covar(
#'   listcovariates = list(1, 1, 1, 1, 1),
#'   mtry = 3
#' )
rand_covar <- function(listcovariates, mtry = 1, randmtry = FALSE) {
  nbcov <- 0

  if (randmtry) {
    if (mtry > 1) {
      stop("mtry is strictly larger than one.
           Decrease it or set randmtry to FALSE")
    }
    while (nbcov == 0) {
      usedcov <- stats::rbinom(
        n = length(listcovariates),
        size = 1,
        prob = mtry
      )
      # usedcov <- sample(c(0, 1),
      #   size = length(listcovariates),
      #   replace = T, prob = c(1 - mtry, mtry)
      # )
      nbcov <- sum(usedcov)
    }
  } else {
    if (mtry > length(listcovariates)) {
      stop("mtry is larger than the number of covariates.")
    }
    if (mtry < 1) {
      stop("mtry is strictly smaller than one. Increase it or set randmtry to TRUE")
    }
    usedcov <- sample(
      c(
        rep(1, mtry),
        rep(0, length(listcovariates) - mtry)
      ),
      replace = FALSE
    )
  }

  return(usedcov)
}


#' Tree prediction for importance
#'
#' @param object A spatial intensity tree returned
#' by tesscovtree of tesstree functions
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @details
#' To compute the importance, we need to keep the structure
#' of a tree and not just the pixel image of the intensity that
#' it outputs.
#'
#' @return A number .....
#' @export
#'
#' @examples
#' vecval0 <- lapply(beisoilres, FUN = function(i) {
#'   c(as.matrix.im(i))
#' })
#' arbre <- tesscovtree(
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
#' predicttree(object = arbre, newdata = c(100, 100))
predicttree <- function(object, newdata, ...) {
  # Test if the covariates are im object
  whichcovim <- unlist(lapply(
    object$listcov,
    spatstat.geom::is.im
  ))
  if (!all(whichcovim)) {
    stop("It appears that in predict.sptree, the covariates of the
             tree are not spatstat im objects")
  }

  # Handles the newdata to be in the correct form
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else if (!spatstat.geom::is.ppp(newdata)) {
    X <- spatstat.geom::ppp(
      x = newdata[1],
      y = newdata[2],
      window = object$X$window
    )
  } else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n == 0) {
      return(NULL)
    }
    X <- newdata
  }

  valsplits <- lapply(object$listcov, FUN = function(j) {
    j[X]
  })

  trees <- lapply(object$tree, FUN = function(i) {
    c(
      i$status,
      i$split_var,
      i$split_val,
      i$intensity_pred,
      i$left_daughter,
      i$right_daughter
    )
  })

  output <- NULL
  for (i in 1:spatstat.geom::npoints(X)) {
    node <- trees[[1]]

    while (node[1] == 1) {
      if (valsplits[[node[2]]][i] < node[3]) {
        child <- node[5]
      } else {
        child <- node[6]
      }
      node <- trees[[child]]
    }
    output[i] <- node[4]
  }

  return(output)
}

#' OOB forest
#'
#' @param forest spforest
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
#' OOBscr(forest)
OOBscr <- function(forest) {
  X <- forest$X # this is always the root

  OOBscr <- lapply(1:length(forest$trees),
    FUN = function(i) {
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
        if (all(!OOBpts)) {
          # If no points in OOBpts, then nothing do do.
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

      return(OOBval)
    }
  )

  logterm <- log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE))
  if (all(is.na(logterm))) {
    output <- NA
  } else {
    output <- sum(logterm, na.rm = TRUE)
  }
  # output <- sum(log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE)), na.rm = TRUE)

  return(output)
}



#' Find best OOB parameter combination
#'
#' @param X The point pattern as a \code{\link[spatstat.geom]{ppp}} .
#' @param listcovariates list of covariates
#' @param params A list containing entries mtry, minpts and, optionally, Ntree.
#' @param ... Other arguments passed to spforest.
#'
#' @details
#' A spatial random forest is generated from the combinations of all the
#' arguments given in params. If Ntree is not given, it is set to 50 by defaults.
#' The arguments ... are passed to the function \code{\link{spforest}}.
#'
#' @return A dataframe with 3 columns: mtry, minpts and OOB
#' @export
#'
#' @examples
#' X <- spatstat.data::bei
#' listcovariates <- spatstat.data::bei.extra
#' params <- list(mtry = c(1, 2), minpts = c(50, 100, 200))
#' OOBoptim(X = X, listcovariates = listcovariates, Ntree = 50, params = params)
OOBoptim <- function(X, listcovariates, params, ...) {
  if (!"mtry" %in% names(params)) {
    stop("The arguments 'params' must have an entry named 'mtry'.")
  }
  if (!"minpts" %in% names(params)) {
    stop("The arguments 'params' must have an entry named 'minpts'.")
  }
  if (!"Ntree" %in% names(params)) {
    nbtree <- 50
  } else {
    nbtree <- params$Ntree
  }

  argu <- expand.grid(params)

  allforest <- mapply(spforest,
    mtry = argu$mtry,
    minpts = argu$minpts,
    Ntree = nbtree,
    MoreArgs = list(
      X = X,
      listcovariates = listcovariates,
      ...
    )
  )

  allOOB <- apply(allforest, 2, OOBscr)

  return(cbind(argu, OOB = allOOB))
}


#' OOB ppm
#'
#' @param forest  spforest
#' @param covariates list of covariates considered
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

  OOBppm <- future.apply::future_lapply(1:length(forest$trees), FUN = function(i) {
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
  })

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