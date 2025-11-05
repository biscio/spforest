#' Importance of one covariate
#'
#' @param forest A spforest object
#' @param id_cov The index in the list of covariates used in the argument
#' \code{listcov} passed in \code{\link[spforest]{spforest}}.
#' @param viptype An integer in \{1,2,3,4\} representing the type of vip used.
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
importance <- function(forest, id_cov, viptype = 4) {
  # listZ <- forest$listcovsp
  X <- forest$X # this is always the root
  Z <- forest$listcov[[id_cov]] # list of cov

  vip_tree <- NULL

  for (i in 1:length(forest$trees)) {
    if (viptype == 3) {
      # Impurity like the randomForest R package - page 6 of the help

      splitvar <- lapply(forest$trees[[i]]$tree, FUN = function(j) {
        j$split_var
      }) |> unlist()

      splitscr <- lapply(forest$trees[[i]]$tree, FUN = function(j) {
        j$scrsplit
      }) |> unlist()

      vip_tree[[i]] <- sum(splitscr[splitvar == id_cov], na.rm = TRUE)
    }

    if (viptype == 4) {
      # Impurity like the randomForest R package - page 6 of the help
      # Look at the decrease of the scores between parent and children

      splitvar <- lapply(forest$trees[[i]]$tree, FUN = function(j) {
        j$split_var
      }) |> unlist()

      splitdcr <- lapply(forest$trees[[i]]$tree, FUN = function(j) {
        j$scrdcr
      }) |> unlist()

      vip_tree[[i]] <- sum(splitdcr[splitvar == id_cov], na.rm = TRUE)
    }


    # Shuffle the chosen covariate value
    dimmat <- Z$dim
    Z$v <- matrix(Z$v[sample.int(length(Z$v))],
      ncol = dimmat[2]
    )

    # OOB sample
    if (forest$p == 0) {
      torm <- unique(forest$pt_intree[[i]])
      if (length(torm) == X$n) { # If all points are drawn in the bootstrap,
        # then nothing do do.
        warning("Some out of bag sample were empty")
        return(0)
      }
      OOBpts <- 1:X$n
      OOBpts <- OOBpts[!(OOBpts %in% torm)]
    } else {
      # vector of same length as number of pts in X
      OOBpts <- (forest$pt_intree[[i]] != 1)
      if (all(!OOBpts)) {
        # If no points in OOBpts, then nothing do do.
        warning("Some out of bag sample were empty")
        return(0)
      }
    }
    Xout <- X[OOBpts]

    listZ_shuf <- forest$listcov
    listZ_shuf[[id_cov]] <- Z

    # Shuffle the tree
    treepert <- forest$trees[[i]]
    treepert$listcov <- listZ_shuf

    ### OOB prediction for shuffled covariate
    pts_pred_OOB_pert <- predicttree(
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

    differr <- pts_pred_OOB_pert - pts_pred_OOB

    if (viptype == 1) {
      # OOB mean square error of the tree
      vip_tree[[i]] <- sqrt(mean(differr^2, na.rm = TRUE))
    }

    if (viptype == 2) {
      # like the randomForest R package - page 6 of the help
      errsd <- sd(differr,
        na.rm = TRUE
      )
      vip_tree[[i]] <- ifelse(errsd == 0,
        mean(abs(differr), na.rm = TRUE),
        mean(abs(differr), na.rm = TRUE) / errsd
      )
    }
  }

  # Return the error of all the trees
  return(unlist(vip_tree))
}


#' Importance of all the covariates
#'
#' @param forest A spforest object
#' @param viptype Always set to 4. Argument passed to
#' \code{\link[spforest]{importance}}.
#' @param treesdetails Boolean. If TRUE, returns a matrix where each column
#' contains the importance of a variable for each trees. Otherwise, return
#' the mean importance of each variable over all trees.
#'
#' @return A vector of the importance of all the covariates of
#' \code{forest} as returned by \code{\link[spforest]{importance}}.
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
#' vip(forest)
vip <- function(forest, viptype = 4, treesdetails = FALSE) {
  if (is.null(forest$listcov)) {
    stop("There is no covariates on which computing the importance.")
  }

  if (viptype != 4) {
    warning("viptype has been set to 4. Other values are not implemented yet.")
  }

  vipvals <- sapply(
    X = seq_along(forest$listcov),
    FUN = function(i) {
      importance(forest,
        id_cov = i,
        viptype = 4
      )
    }
  )

  if (!is.null(names(forest$listcov))) {
    if (is.matrix(vipvals)) {
      colnames(vipvals) <- names(forest$listcov)
    } else {
      names(vipvals) <- names(forest$listcov)
    }
  }

  if (treesdetails) {
    return(vipvals)
  } else {
    output <- colMeans(vipvals)
    return(output)
  }
}


#' Vip barplot of forest
#'
#' @param x  A \code{\link{spforest.object}} with \code{listcov} not \code{NULL}.
#' @param sorted Logical. If TRUE, the variables are sorted by importance.
#' @param viptype An integer in \{1,2,3,4\} passed to \code{\link[spforest]{importance}}.
#' @param ... Arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @details
#' First the function \code{\link[spforest]{importance}} is called to
#' compute the importance of each covariate in each tree \code{x}.
#' Then, a barplot is drawn to visualize the distribution
#' of the importance of each covariate, averaged over all the trees.
#' @return Variable importance plot
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
#' vipplot(forest, sorted = TRUE)
vipplot <- function(x,
                    sorted = FALSE,
                    viptype = 4, ...) {
  
  if (is(x, "spforest")) {
    output <- vip(x, viptype = viptype, treesdetails = FALSE)
  } 
  
  if (is.vector(x)) {
    output <- x
    if (is.null(names(output))) {
      names(output) <- paste0("var", seq_along(x))
    }
  }
  
  if (is.matrix(x)) {
    graphics::boxplot(x, ...)
    return(invisible(x))
  }

  if (sorted) {
    output <- sort(output,
      decreasing = TRUE
    )
  }

  graphics::barplot(output,
    names.arg = names(output),
    ...
  )

  return(invisible(output))
}


#' Boxplot forest
#'
#' Plot the boxplot of the importance of each covariate
#' in each tree of the random forest intensity
#'
#' @usage boxplot(x = forest)
#'
#' @param x A \code{\link{spforest.object}} with \code{listcov} not \code{NULL}.
#' @param ... Arguments passed to \code{\link[graphics]{boxplot}}.
#'
#' @details
#' First the function \code{\link[spforest]{importance}} is called to
#' compute the importance of each covariate in each tree \code{x}.
#' Then, a boxplot is drawn to visualize the distribution of the
#' importance of each covariate.
#'
#' @return A list returned by the function \code{\link[graphics]{boxplot}}.
#' @param viptype An integer in \{1,2,3,4\} passed to \code{\link[spforest]{importance}}.
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
#' boxplot(forest)
boxplot.spforest <- function(x, viptype = 4, ...) {
  vipval <- sapply(X = seq_along(x$listcov), FUN = function(i) {
    importance(x, id_cov = i, viptype = viptype)
  })

  df <- data.frame(
    vipval = as.vector(vipval),
    var = rep(names(x$listcov), each = length(x$trees))
  )

  graphics::boxplot(vipval ~ var,
    data = df,
    xlab = "Variables",
    ylab = "Variable Importance",
    ...
  )
}
