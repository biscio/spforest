#' Randomise input for tree
#'
#' @param listcovariates list of covariatesd used
#' @param mtry selecting proba
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


#' Spatial Intensity tree function
#'
#' @inheritParams splitcell
#' @param mtry Probability of choosing a covariable
#' @param threshold minimum threshold to allow to split cell
#' @param inforest Logical. Indicates if the function is run inside RforestPP.
#'
#' @return An object of class \code{sptree}.
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   threshold = 1000,
#'   score = "lcv",
#'   listcovariates = list(
#'     grad = spatstat.data::bei.extra$grad,
#'     elev = spatstat.data::bei.extra$elev
#'   ),
#'   mtry = 1,
#'   tol = Inf,
#'   minpts = 50,
#'   minsplitq = 0.5,
#'   maxsplitq = 0.5
#' )
treerec <- function(X,
                    score = "lcv",
                    threshold = spatstat.geom::area(X) / 1e4,
                    listcovariates = NULL,
                    mtry = 1,
                    tol = Inf,
                    minpts = spatstat.geom::npoints(X) / 10,
                    minsplitq = 0.5,
                    maxsplitq = 0.5,
                    inforest = F) {
  # Sanity checks
  if (threshold <= 0 & minpts <= 0) {
    stop("Either threshold or minpts must be strictly greater than 0.")
  }

  # The code works quicker when the windows is on a mask
  spatstat.geom::Window(X) <- spatstat.geom::as.mask(spatstat.geom::Window(X),
    xy = listcovariates[[1]]
  )

  # Initiating the root
  root <- list(
    nodeID = 1,
    nodePP = X,
    left_daughter = NA,
    right_daughter = NA,
    split_var = NA,
    split_val = NA,
    status = 1,
    intensity_pred = spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
    scr_parent = NA,
    improvement = NULL,
    already_split = FALSE,
    whystop = NULL
  )

  # Check if there is something to do on the initial point pattern.
  # if (spatstat.geom::area.owin(X$window) <= threshold |
  #   spatstat.geom::npoints(X) <= minpts) {
  #   root$status <- 0
  #   return(list(root))
  # }
  if (spatstat.geom::area.owin(X$window) <= threshold |
    spatstat.geom::npoints(X) <= minpts) {
    root$status <- 0
    output <- list(
      tree = list(root),
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      im = as.im(spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
        W = X$window
      )
    )
    class(output) <- "sptree" # For when I will define class
    return(output)
  }


  intensity_tree <- list(root)

  # Initialise the while loop
  k <- 0
  knew <- 1
  # Core computation of the tree
  while (k != knew) {
    k <- length(intensity_tree)

    already_split_node <- sapply(intensity_tree,
      FUN = function(j) {
        j$already_split
      }
    )

    # Never used visibly
    # if (all(already_split_node) & k > 1) {
    #   return(intensity_tree)
    # }

    for (i in (1:k)[!already_split_node]) {
      ## Select randomly covariates
      usedcov <- rand_covar(
        listcovariates = listcovariates,
        mtry = mtry
      )

      # Split the cell, if the split is valid under the chosen parameters
      res.split <- splitcell(
        X = intensity_tree[[i]]$nodePP,
        score = score,
        listcovariates = listcovariates,
        usecovariates = usedcov,
        thres.cell = threshold,
        minpts = minpts,
        tol = tol,
        imp = intensity_tree[[i]]$improvement,
        minsplitq = minsplitq,
        maxsplitq = maxsplitq
      )

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
        childleft <- list(
          nodeID = knew + 1,
          nodePP = res.split$PPleft,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = spatstat.geom::intensity.ppp(res.split$PPleft),
          scr_parent = res.split$scr_parent,
          improvement = res.split$improvement,
          already_split = FALSE,
          whystop = NULL
        )
        childright <- list(
          nodeID = knew + 2,
          nodePP = res.split$PPright,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = spatstat.geom::intensity.ppp(res.split$PPright),
          scr_parent = res.split$scr_parent,
          improvement = res.split$improvement,
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

  # Compute the image
  idterm <- sapply(intensity_tree, FUN = function(i) {
    i$status
  })

  listmask <- lapply(intensity_tree[idterm == 0], FUN = function(i) {
    i$nodePP$window$m # It is a mask usually but if no split is done then it is not.
  })

  intensity_prediction <- sapply(intensity_tree[idterm == 0], FUN = function(i) {
    i$intensity_pred
  })

  matim <- mapply("*", listmask, intensity_prediction, SIMPLIFY = FALSE)

  imoutput <- spatstat.geom::as.im(Reduce("+", matim),
    W = intensity_tree[[1]]$nodePP
  )

  # Remove all the intermediary PP in nodePP
  for (i in seq_along(intensity_tree)) {
    intensity_tree[[i]]$nodePP <- NULL
    intensity_tree[[i]]$already_split <- NULL
    # x$tree[[i]]$improvement <- NULL
  }

  if (inforest) {
    output <- list(
      tree = intensity_tree,
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      im = imoutput
    )
  } else {
    output <- list(
      tree = intensity_tree,
      X = X,
      namecov = names(listcovariates),
      namelist = as.character(match.call()[4]),
      listcov = listcovariates,
      im = imoutput
    )
  }

  class(output) <- "sptree" # For when I will define class

  return(output)
}





# Function that to a nodeID find the parent ID.
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
