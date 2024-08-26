# Original point pattern
X <- bei
# Covariates
listcovariates <- Rsandbox::beisoilres
# Choose a score
score <- "lcv"
# Used a subset each time
# usecovariates <- rep(c(1, 0, 1), each = 5)
usecovariates <- rep(1, 15)
# limit small cell
threshold <- 100
# Limit number points
minpts <- 500
# randomisation covariates
mtry <- 1

arbre <- treerec(
  X = spatstat.data::bei,
  threshold = 100,
  listcovariates = Rsandbox::beisoilres,
  mtry = 1,
  minpts = 500
)


g <- function() {

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
dimcov <- listcovariates[[1]]$dim
covrangex <- listcovariates[[1]]$xrange
covrangey <- listcovariates[[1]]$yrange

#### TODO, add X as arguments, for consistency
splitcell3 <- function(X,
                       valpts,
                       vecval,
                       usecovariates,
                       areapixel,
                       threshold,
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
      tempsublvl <- im(
        matrix(sublvl[[i]],
          nrow = dimcov[1],
          ncol = dimcov[2], byrow = F
        ),
        xrange = covrangex, yrange = covrangey
      )
      n1 <- npoints(X[tempsublvl])
      n2 <- npoints(X[!tempsublvl])
      scr_cov[i] <- score.split(
        n1 = n1,
        n2 = n2,
        W1area = Wsub,
        W2area = Wsup
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

    nsub <- spatstat.geom::npoints(X[imsublvl])
    nsup <- spatstat.geom::npoints(X[!imsublvl])

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
      nsub = nsub,
      suplevels = suplevels,
      nsup = nsup,
      whystop = NULL
    )
  }
}

# A<-splitcell3(
#   valpts = valpts,
#   vecval = vecval,
#   usecovariates = usecovariates,
#   areapixel = areapixel,
#   threshold = threshold
# )

root <- list(
  nodeID = 1,
  nodeCov = vecval,
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

    if (intensity_tree2[[i]]$nX <= minpts) {
      res.split <- "Not enough points to attempt to split"
    } else {
      # Split the cell, if the split is valid under the chosen parameters
      res.split <- splitcell3(
        X = X, 
        valpts = valpts,
        vecval = intensity_tree2[[i]]$nodeCov,
        usecovariates = usedcov,
        areapixel = areapixel,
        threshold = threshold,
        dimcov = dimcov,
        covrangex = covrangex,
        covrangey = covrangey
      )
    }

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

      # Define the children
      areasub <- (areapixel * sum(!is.na(res.split$sublevels[[res.split$split_var]])))
      areasup <- (areapixel * sum(!is.na(res.split$suplevels[[res.split$split_var]])))

      childleft <- list(
        nodeID = knew + 1,
        nodeCov = res.split$sublevels,
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
      intensity_tree2 <- append(
        intensity_tree2,
        list(childleft, childright)
      )
    }
    knew <- length(intensity_tree2)
  }

  # sapply(intensity_tree2, FUN = function(i) {
  #   i$right_daughter
  # })
  # sapply(intensity_tree2, FUN = function(i) {
  #   i$left_daughter
  # })
  # sapply(intensity_tree2, FUN = function(i) {
  #   i$status
  # })
  # sapply(intensity_tree2, FUN = function(i) {
  #   i$already_split
  # })
  # sapply(intensity_tree2, FUN = function(i) {
  #   i$nX
  # })
}

idterm <- sapply(intensity_tree2, FUN = function(i) {
  i$status
})
A <- lapply(intensity_tree2[idterm == 0],
  FUN = function(i) {
    A <- ifelse(!is.na(i$nodeCov[[1]]),
      i$intensity_pred,
      0
    )
    return(im(
      matrix(A,
        nrow = dimcov[1],
        ncol = dimcov[2], byrow = F
      ),
      xrange = covrangex, 
      yrange = covrangey
    ))
  }
)
return(spatstat.geom::as.im(Reduce("+", A), W = X$window))
}
plot(spatstat.geom::as.im(Reduce("+", A), W = X$window))

f <- function(){treerec(
  X = spatstat.data::bei,
  threshold = 1000,
  listcovariates = beisoilres,
  mtry = 1,
  minpts = 100
)
}

library(microbenchmark)
microbenchmark(f(),g())

#### treerec2 code old ----

#
# plot(im(matrix(intensity_tree2[idterm == 0][[1]]$nodeCov[[1]],
#                nrow = dimcov[1],
#                ncol = dimcov[2], byrow = F),
#         xrange = covrangex, yrange = covrangey) )
#

# # Compute the image
# idterm <- sapply(intensity_tree2, FUN = function(i) {
#   i$status
# })

### TODO
### HOW TO COMPUTE THE IMAGE NOW?
#
# listmask <- lapply(intensity_tree2[idterm == 0],
# FUN = function(i) {
#   i$nodePP$window$m # It is a mask usually but if no split is done then it is not.
# })
#
# intensity_prediction <- sapply(intensity_tree2[idterm == 0],
# FUN = function(i) {
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
# }
#
#
# treerec2(
#     X = spatstat.data::bei,
#     minpts = 50,
#     threshold = 1000,
#     listcovariates = bei.extra,
#     mtry = 1,
#   )
