#' Splitting cell
#'
#' @param X The observed data point pattern.
#' @param valpts A list. Values of each covariates at the points of \code{X}.
#' @param vecval A list. Values of the covariates at each pixel.
#' @param usecovariates  Vector of 0 and 1 values indicating which covariates
#' to choose
#' @param areapixel The pixel area used in each covariates.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param threshold Minimum threshold to allow to split cell.
#' @param dimcov The element \code{dim} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#' @param covrangex The element \code{xrange} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#' @param covrangey The element \code{yrange} of the covariates which are
#' \code{\link[spatstat.geom]{im}}.
#'
#' @return A list with  split_var = split_var,
#' split_val = split_val,
#' sublevels = sublevels,
#' nsub = nsub,
#' suplevels = suplevels,
#' nsup = nsup,
#' whystop = NULL
#' @export
#'
#' @examples
#' listcovariates <- Rsandbox::beisoilres
#' X <- bei
#' valpts <- lapply(listcovariates,
#'   FUN = function(i) {
#'     i[X]
#'   }
#' )
#' vecval <- lapply(Rsandbox::beisoilres, FUN = function(i) {
#'   if (!spatstat.geom::is.im(i)) {
#'     stop("Elements of listcovar must be an im object")
#'   }
#'   c(as.matrix.im(i))
#' })
#' usecovariates <- rep(1, 15)
#' areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
#' score <- "lcv"
#' dimcov <- listcovariates[[1]]$dim
#' covrangex <- listcovariates[[1]]$xrange
#' covrangey <- listcovariates[[1]]$yrange
#' A<-splitcell2(
#'   X=X,
#'   valpts = valpts,
#'   vecval = vecval,
#'   usecovariates = usecovariates,
#'   dimcov = dimcov,
#'   covrangex = covrangex,
#'   covrangey = covrangey,
#'   areapixel = areapixel,
#'   threshold = 100
#' )
splitcell2 <- function(X,
                       valpts,
                       vecval,
                       usecovariates,
                       areapixel,
                       score = "lcv",
                       threshold = spatstat.geom::area(X) / 1e4,
                       dimcov,
                       covrangex,
                       covrangey) {
  whynot <- NULL
  vecvalused <- vecval[usecovariates == 1]
  # listcovar <- listcovariates[usecovariates == 1]

  mediancov <- lapply(vecvalused, FUN = function(j) {
    stats::median(j, na.rm = TRUE)
  })

  sublvl <- lapply(1:length(mediancov),
    FUN = function(j) {
      vecvalused[[j]] < mediancov[[j]]
    }
  )

  # Determination of level sets for all covariables
  scr_cov <- NULL
  for (i in 1:sum(usecovariates)) {
    Wsub <- areapixel * sum(sublvl[[i]], na.rm = TRUE)
    Wsup <- areapixel * sum(!sublvl[[i]], na.rm = TRUE)

    # Test if they are too small and computation of the score if not
    if (Wsub <= threshold | Wsup <= threshold) {
      scr_cov[i] <- -Inf
    } else {
      n1 <- sum(valpts[[i]] < mediancov[[i]], na.rm = T)
      n2 <- sum(valpts[[i]] >= mediancov[[i]], na.rm = T)
      # tempsublvl <- im(
      #   matrix(sublvl[[i]],
      #     nrow = dimcov[1],
      #     ncol = dimcov[2], byrow = F
      #   ),
      #   xrange = covrangex, yrange = covrangey
      # )
      # tempsuplvl <- im(
      #   matrix(!sublvl[[i]],
      #          nrow = dimcov[1],
      #          ncol = dimcov[2], byrow = F
      #   ),
      #   xrange = covrangex, yrange = covrangey
      # )
      # n1 <- npoints(X[tempsublvl])
      # n2 <- npoints(X[tempsuplvl])
      scr_cov[i] <- score.split(
        n1 = n1,
        n2 = n2,
        W1area = Wsub,
        W2area = Wsup,
        score = score
      )
    }
  }

  ## Go out if all the score are -Inf
  if (all(is.infinite(scr_cov))) {
    whynot <- c("All split scores are -Inf, cell too small")
  }

  id_best_scr <- sort(scr_cov,
    index.return = T,
    decreasing = T
  )$ix[1]

  if (!is.null(whynot)) {
    return(whynot)
  } else {
    split_var <- which(usecovariates == 1)[id_best_scr]
    split_val <- mediancov[[id_best_scr]]
    splitsub <- (vecval[[split_var]] < split_val)

    imsublvl <- im(
      matrix(splitsub,
        nrow = dimcov[1],
        ncol = dimcov[2], byrow = F
      ),
      xrange = covrangex, yrange = covrangey
    )

    subvalpts <- (valpts[[split_var]] < split_val)
    nsub <- sum(subvalpts, na.rm=T)
    nsup <- sum(!subvalpts, na.rm=T)
    
    valptssub <- lapply(valpts, FUN = function(j){
      ifelse(subvalpts, 
             j, NA)
    })
    valptssup <- lapply(valpts, FUN = function(j){
      ifelse(!subvalpts, 
             j, NA)
    })
    
    # nsub <- spatstat.geom::npoints(X[imsublvl])
    # nsup <- spatstat.geom::npoints(X[!imsublvl])

    splitsup <- !splitsub
    splitsub[!splitsub] <- NA
    splitsup[!splitsup] <- NA
    sublevels <- lapply(vecval, FUN = function(j) {
      return(j * splitsub)
    })
    suplevels <- lapply(vecval, FUN = function(j) {
      return(j * splitsup)
    })

    nodeChilds <- list(
      split_var = split_var,
      split_val = split_val,
      valptssub = valptssub, 
      valptssup = valptssup,
      sublevels = sublevels,
      nsub = nsub,
      suplevels = suplevels,
      nsup = nsup,
      whystop = NULL
    )
  }
}




#' Spatial Intensity tree function (improved)
#'
#' @param X The observed data point pattern.
#' @param listcovariates A list with all the covariates used for the tree. The
#' covariates must be given as im object from the package spatstat.
#' @param minpts A positive integer. 
#' The minimum number of points allowed to try to split a cell.
#' @param mtry Probability of choosing a covariable.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param threshold Minimum threshold to allow to split cell.
#' @param inforest Logical. Indicates if the function is run in a forest or not. 
#'
#' @return An object of class \code{sptree}.
#' @export
#'
#' @examples
#' mytree <- intensitytree(
#'   X = spatstat.data::bei,
#'   listcovariates = beisoilres,
#'   mtry = 1,
#'   minpts = 500
#' )
#' plot(mytree)
intensitytree <- function(X,
                          listcovariates,
                          minpts = 500,
                          mtry = 1,
                          score = "lcv",
                          threshold = spatstat.geom::area(X) / 1e4,
                          inforest = F) {
  valpts <- lapply(listcovariates,
    FUN = function(i) {
      i[X]
    }
  )

  areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep

  vecval <- lapply(listcovariates, FUN = function(i) {
    if (!spatstat.geom::is.im(i)) {
      stop("Elements of listcovar must be an im object")
    }
    c(as.matrix.im(i))
  })

  
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
    whystop = NULL
  )

  dimcov <- listcovariates[[1]]$dim
  covrangex <- listcovariates[[1]]$xrange
  covrangey <- listcovariates[[1]]$yrange

  if (spatstat.geom::area.owin(X$window) <= threshold |
    spatstat.geom::npoints(X) <= minpts) {
    root$status <- 0
    output <- list(
      tree = list(root),
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      im = spatstat.geom::as.im(spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
        W = X$window
      )
    )
    class(output) <- "sptree"
    return(output)
  }

  # Check if there is something to do on the initial point pattern.
  if (spatstat.geom::area.owin(X$window) <= threshold |
    spatstat.geom::npoints(X) <= minpts) {
    root$status <- 0
    output <- list(
      tree = list(root),
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
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
        mtry = mtry
      )

      if (intensity_tree[[i]]$nX <= minpts) {
        res.split <- "Not enough points to attempt to split"
      } else {
        # Split the cell, if the split is valid under the chosen parameters
        res.split <- splitcell2(
          X = X,
          valpts = intensity_tree[[i]]$nodeValpts,
          vecval = intensity_tree[[i]]$nodeCov,
          usecovariates = usedcov,
          areapixel = areapixel,
          score = score,
          threshold = threshold,
          dimcov = dimcov,
          covrangex = covrangex,
          covrangey = covrangey
        )
      }

      if (is.character(res.split)) {
        intensity_tree[[i]]$status <- 0
        intensity_tree[[i]]$already_split <- TRUE
        intensity_tree[[i]]$whystop <- res.split
      } else {
        # Update the parent
        intensity_tree[[i]]$left_daughter <- knew + 1
        intensity_tree[[i]]$right_daughter <- knew + 2
        intensity_tree[[i]]$split_var <- res.split$split_var
        intensity_tree[[i]]$split_val <- res.split$split_val
        intensity_tree[[i]]$already_split <- TRUE

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
          whystop = NULL
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
          whystop = NULL
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

  if (inforest) {
    output <- list(
      tree = intensity_tree,
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      im = spatstat.geom::as.im(Reduce("+", patchworks), W = X$window)
    )
  } else {
    output <- list(
      tree = intensity_tree,
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      listcov = listcovariates,
      im = spatstat.geom::as.im(Reduce("+", patchworks), W = X$window)
    )
  }

  class(output) <- "sptree" # For when I will define class

  return(output)

  # return(spatstat.geom::as.im(Reduce("+", patchworks), W = X$window))
}


# plot(spatstat.geom::as.im(Reduce("+", A), W = X$window))
#
# f <- function() {
#   treerec(
#     X = spatstat.data::bei,
#     threshold = 1000,
#     listcovariates = beisoilres,
#     mtry = 1,
#     minpts = 100
#   )
# }
#
# g <- function() {
#   intensitytree(
#     X = spatstat.data::bei,
#     listcovariates = beisoilres,
#     mtry = 1,
#     minpts = 100,
#     threshold = 1000
#   )
# }
#
# A <- f()
# B <- g()
# plot(A)
# plot(B)

#
# library(microbenchmark)
# microbenchmark(f(), g())
#
# library(profvis)
# profvis(f())
# profvis(g())

#
# plot(im(matrix(intensity_tree[idterm == 0][[1]]$nodeCov[[1]],
#                nrow = dimcov[1],
#                ncol = dimcov[2], byrow = F),
#         xrange = covrangex, yrange = covrangey) )
#
