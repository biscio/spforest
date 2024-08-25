# Original point pattern
X <- bei
# Covariates
listcovariates <- Rsandbox::beisoilres
# Choose a score
score <- "lcv"
# Used a subset each time
usecovariates <- rep(c(1, 0, 1), each = 5)
# limit small cell
threshold <- 100



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

# In while loop in tree function

splitcell3 <- function(valpts,
                       vecval, 
                       usecovariates, 
                       areapixel,
                       threshold) {
  whynot = NULL
  vecvalused <- vecval[usecovariates == 1]
  listcovar <- listcovariates[usecovariates == 1]

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
    Wleft <- areapixel * sum(sublvl[[i]], na.rm = TRUE)
    Wright <- areapixel * sum(!sublvl[[i]], na.rm = TRUE)

    # Test if they are too small and computation of the score if not
    if (Wleft <= threshold | Wright <= threshold) {
      scr_cov[i] <- -Inf
    } else {
      n1 <- sum(valpts[usecovariates == 1][[i]] < mediancov[[i]], na.rm = TRUE)
      scr_cov[i] <- score.split(
        n1 = n1,
        n2 = spatstat.geom::npoints(X) - n1,
        W1area = Wleft,
        W2area = Wright
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
      sublevels = sublevels,
      suplevels = suplevels,
      whystop = NULL
    )
  }
}

root <- list(
  nodeID = 1,
  nodeCov = listcovariates,
  left_daughter = NA,
  right_daughter = NA,
  split_var = NA,
  split_val = NA,
  status = 1,
  intensity_pred = spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
  already_split = FALSE,
  whystop = NULL
)

# Check if there is something to do on the initial point pattern.
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
  class(output) <- "sptree"
  return(output)
}

intensity_tree2 <- list(root)

# Initialise the while loop
k <- 0
knew <- 1

while (k != knew) {TODO stop("TODO")}


#### treerec2 code----
treerec2 <- function(X,
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

  areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep

  # Initiating the root

  # vecval <- lapply(listcovariates, FUN = function(i) {
  #   if (!spatstat.geom::is.im(i)) {
  #     stop("Elements of listcovar must be an im object")
  #   }
  #   c(i$v)
  # })

  root <- list(
    nodeID = 1,
    nodePP = rep(T, spatstat.geom::npoints(X)),
    nodeCov = listcovariates,
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
    class(output) <- "sptree"
    return(output)
  }

  intensity_tree2 <- list(root)

  # Initialise the while loop
  k <- 0
  knew <- 1
  dimcov <- listcovariates[[1]]$dim
  covrangex <- listcovariates[[1]]$xrange
  covrangey <- listcovariates[[1]]$yrange
  # Core computation of the tree
  ### TODO
  ### CHECK THAT WE CAN STOP !!!!
  while (k != knew) {
    k <- length(intensity_tree2)

    already_split_node <- sapply(intensity_tree2,
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

      # Split the cell, if the split is valid under the chosen parameters
      res.split <- splitcell2(
        X = X[intensity_tree2[[i]]$nodePP],
        score = score,
        listcovariates = intensity_tree2[[i]]$nodeCov,
        usecovariates = usedcov,
        thres.cell = threshold,
        minpts = minpts,
        tol = tol,
        imp = intensity_tree2[[i]]$improvement,
        minsplitq = minsplitq,
        maxsplitq = maxsplitq
      )

      if (is.character(res.split)) {
        intensity_tree2[[i]]$status <- 0
        intensity_tree2[[i]]$already_split <- TRUE
        intensity_tree2[[i]]$whystop <- res.split
      } else {
        # Update the parent
        intensity_tree2[[i]]$left_daughter <- knew + 1
        intensity_tree2[[i]]$right_daughter <- knew + 2
        intensity_tree2[[i]]$split_var <- res.split$split_var
        intensity_tree2[[i]]$split_val <- res.split$split_val
        intensity_tree2[[i]]$already_split <- TRUE

        newsublvl <- lapply(res.split$sublevels, FUN = function(jj) {
          A <- matrix(jj, nrow = dimcov[1], ncol = dimcov[2], byrow = F)
          im(A, xrange = covrangex, yrange = covrangey)
        })
        newsuplvl <- lapply(res.split$suplevels, FUN = function(jj) {
          A <- matrix(jj, nrow = dimcov[1], ncol = dimcov[2], byrow = F)
          im(A, xrange = covrangex, yrange = covrangey)
        })

        # Define the children
        childleft <- list(
          nodeID = knew + 1,
          nodePP = res.split$PPleft,
          nodeCov = newsublvl,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = sum(res.split$PPleft) / (areapixel * sum(res.split$Wleft)),
          scr_parent = res.split$scr_parent,
          improvement = res.split$improvement,
          already_split = FALSE,
          whystop = NULL
        )
        childright <- list(
          nodeID = knew + 2,
          nodePP = !res.split$PPleft,
          nodeCov = newsuplvl,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = sum(!res.split$PPleft) / (areapixel * sum(!res.split$Wleft)),
          scr_parent = res.split$scr_parent,
          improvement = res.split$improvement,
          already_split = FALSE,
          whystop = NULL
        )
        # append the children
        intensity_tree2 <- append(
          intensity_tree2,
          list(childleft, childright)
        )
      }
      knew <- length(intensity_tree2)
    }
    # intensity_tree2
    sapply(intensity_tree2, FUN = function(i) {
      i$right_daughter
    })
    sapply(intensity_tree2, FUN = function(i) {
      i$left_daughter
    })
    sapply(intensity_tree2, FUN = function(i) {
      sum(i$nodePP)
    })
    sapply(intensity_tree2, FUN = function(i) {
      i$status
    })
    sapply(intensity_tree2, FUN = function(i) {
      i$already_split
    })
  }

  # # Compute the image
  # idterm <- sapply(intensity_tree2, FUN = function(i) {
  #   i$status
  # })

  ### TODO
  ### HOW TO COMPUTE THE IMAGE NOW?
  #
  # listmask <- lapply(intensity_tree2[idterm == 0], FUN = function(i) {
  #   i$nodePP$window$m # It is a mask usually but if no split is done then it is not.
  # })
  #
  # intensity_prediction <- sapply(intensity_tree2[idterm == 0], FUN = function(i) {
  #   i$intensity_pred
  # })
  #
  # matim <- mapply("*", listmask, intensity_prediction, SIMPLIFY = FALSE)
  #
  # imoutput <- spatstat.geom::as.im(Reduce("+", matim),
  #                                  W = intensity_tree2[[1]]$nodePP
  # )
  #
  # # Remove all the intermediary PP in nodePP
  # for (i in seq_along(intensity_tree2)) {
  #   intensity_tree2[[i]]$nodePP <- NULL
  #   intensity_tree2[[i]]$already_split <- NULL
  #   # x$tree[[i]]$improvement <- NULL
  # }
  #
  # if (inforest) {
  #   output <- list(
  #     tree = intensity_tree2,
  #     X = X,
  #     namecov = names(listcovariates),
  #     namelist = as.character(match.call()[4]),
  #     im = imoutput
  #   )
  # } else {
  #   output <- list(
  #     tree = intensity_tree2,
  #     X = X,
  #     namecov = names(listcovariates),
  #     namelist = as.character(match.call()[4]),
  #     listcov = listcovariates,
  #     im = imoutput
  #   )
  # }
  #
  # class(output) <- "sptree" # For when I will define class
  #
  # return(output)
}
#
#
# treerec2(
#     X = spatstat.data::bei,
#     minpts = 50,
#     threshold = 1000,
#     listcovariates = bei.extra,
#     mtry = 1,
#   )
