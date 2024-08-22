splitcell2 <- function(X,
                       score = "lcv",
                       listcovariates = NULL,
                       usecovariates = rep(1, length(listcovariates)),
                       thres.cell = 25,
                       minpts = 50,
                       tol = Inf,
                       imp = NULL,
                       minsplitq = 0.5,
                       maxsplitq = 0.5) {
  whynot <- NULL
  # Check for incorrect input
  stopifnot(spatstat.geom::is.ppp(X))
  stopifnot(spatstat.geom::is.im(listcovariates[[1]]))

  # Cas pathologique si pas de covariable - no split
  if (sum(usecovariates) == 0 | spatstat.geom::npoints(X) <= minpts) { # normally I always have one
    whynot <- c("no covariate or less than minpts")
    # return(NULL)
  }

  listcovar <- listcovariates[usecovariates == 1]
  vecval <- lapply(listcovar, FUN = function(i) {
    c(i$v)
  })

  valpts <- lapply(listcovar, FUN = function(i) {
    i[X]
  })

  scr_cov <- NULL
  id_best_scr <- NULL
  scr_parent <- score.pp(X = X, score = score)

  # Calculate the splitting scores of all covariates
  split_scrs <- lapply(vecval, FUN = function(j) {
    if (minsplitq != maxsplitq) {
      rand_q <- stats::runif(n = 1, min = minsplitq, max = maxsplitq)
      return(stats::quantile(j, probs = rand_q))
    } else {
      return(stats::median(j, na.rm = TRUE))
    }
  })

  # Remove names of columns due to quantile fct
  names(split_scrs) <- NULL

  leftlvl0 <- lapply(1:length(split_scrs),
    FUN = function(j) {
      vecval[[j]] < split_scrs[[j]]
    }
  )

  areapixel <- listcovar[[1]]$xstep * listcovar[[1]]$ystep

  # Determination of level sets for all covariables
  for (i in 1:length(listcovar)) {
    if (!spatstat.geom::is.im(listcovar[[i]])) {
      stop("Elements of listcovar must be an im object")
    }
    maxsplitq = 0.5
    Wleft <- areapixel * sum(leftlvl0[[i]])
    Wright <- areapixel * sum(!leftlvl0[[i]])

    # Test if they are too small and computation of the score if not
    if (Wleft <= thres.cell |
      Wright <= thres.cell) {
      scr_cov[i] <- -Inf
    } else {
      n1 <- sum(valpts[[i]]<split_scrs[[i]])
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
    # return(NULL)
  }

  ## If there in the code
  ## not all the score are -Inf so there is one best score
  id_best_scr <- sort(scr_cov,
    index.return = T,
    decreasing = T
  )$ix[1]

  # Look at the parents
  # Look if there was any improvements on the score in the last splits
  newimp <- append(
    imp,
    scr_cov[id_best_scr] > scr_parent
  )
  L_imp <- length(newimp)

  if (L_imp >= tol & all(!utils::tail(newimp, n = tol))) {
    whynot <- c("no improvement since tol splits")
    # return(NULL)
  }

  # For split_var,
  # I determine the index of the split var among all covariables

  if (!is.null(whynot)) {
    return(whynot)
  } else {
    nodeChilds <- list(
      PPleft = valpts[[id_best_scr]] < split_scrs[[id_best_scr]],
      Wleft =  (vecval[[id_best_scr]] < split_scrs[[id_best_scr]]),
      split_var = which(usecovariates == 1)[id_best_scr],
      split_val = split_scrs[id_best_scr],
      improvement = newimp,
      scr_parent = scr_parent,
      whystop = NULL
    )
  }
  return(nodeChilds)
}

# 
# f = function(){ splitcell2(
#   X = spatstat.data::bei,
#   score = "lcv",
#   listcovariates = beisoilres,
#   usecovariates = rep(1,15),
#   thres.cell = 100,
#   minpts = 10,
#   tol = Inf,
#   imp = NULL,
#   minsplitq = 0.5,
#   maxsplitq = 0.5
# )}
# g <- function() {splitcell(
#   X = spatstat.data::bei,
#   score = "lcv",
#   listcovariates = beisoilres,
#   usecovariates = rep(1,15),
#   thres.cell = 100,
#   minpts = 10,
#   tol = Inf,
#   imp = NULL,
#   minsplitq = 0.5,
#   maxsplitq = 0.5
# )}
# 
# library(microbenchmark)
# 
# microbenchmark(f(),g())


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
  
  areapixel <-  listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
  
  # The code works quicker when the windows is on a mask
  spatstat.geom::Window(X) <- spatstat.geom::as.mask(spatstat.geom::Window(X),
                                                     xy = listcovariates[[1]]
  )
  
  
  #### TODO 
  #### 
  ####  in place of a ppp in nodePP, 
  ####  put a logical vector with T if the point is present.
  #### 
  #### 
  # Initiating the root
  root <- list(
    nodeID = 1,
    nodePP = rep(T, spatstat.geom::npoints(X)),
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
  
  intensity_tree <- list(root)
  
  # Initialise the while loop
  k <- 0
  knew <- 1
  # Core computation of the tree
  ### TODO 
  ### CHECK THAT WE CAN STOP !!!!
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
      res.split <- splitcell2(
        X = X[intensity_tree[[i]]$nodePP],
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
          intensity_pred = sum(res.split$PPleft)/(areapixel*sum(res.split$Wleft)),
          scr_parent = res.split$scr_parent,
          improvement = res.split$improvement,
          already_split = FALSE,
          whystop = NULL
        )
        childright <- list(
          nodeID = knew + 2,
          nodePP = !res.split$PPleft,
          left_daughter = NA,
          right_daughter = NA,
          split_var = NA,
          split_val = NA,
          status = 1,
          intensity_pred = sum(!res.split$PPleft)/(areapixel*sum(!res.split$Wleft)),
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
  
  ### TODO 
  ### HOW TO COMPUTE THE IMAGE NOW?
  
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

