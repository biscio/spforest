#' Printing spatial intensity tree
#'
#' @param x A spatial intensity tree return by treerec function
#' @param ... Additional arguments
#'
#' @return A description of the random intensity tree with
#' the number of points and covariates used.
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   minpts = 300,
#'   mtry = 1
#' )
#' print(arbre)
print.sptree <- function(x, ...) {
  namecov <- x$namecov
  nb_termnode <- sum(sapply(x$tree, function(i) {
    i$status
  }) == 0)

  cat(paste(
    "Intensity tree estimate of point patterns with",
    x$X$n, "points.\n\n"
  ))

  cat(paste(length(namecov), "covariables used, with names: "))
  cat(namecov, "\n")

  cat(
    "Spatial intensity tree with", length(x$tree),
    "nodes and", nb_termnode, "terminal nodes.", "\n"
  )
  cat("The used covariates were in the object", x$namelist, ".")
}

#' Summary of a spatial intensity tree
#'
#' @param object A spatial intensity tree return by treerec function
#' @param fulltree Should we print all the column of the matrix?
#' @param ... Additional arguments
#'
#' @return A dataframe representing each node and spit of the intensity tree.
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   minpts = 300,
#'   mtry = 1
#' )
#' summary(arbre)
summary.sptree <- function(object, fulltree = F, ...) {
  if (fulltree) {
    output <- Reduce(rbind, object$tree)
    rownames(output) <- NULL
  } else {
    # Just remove the imp and nodeID
    output <- Reduce(rbind, object$tree)[, c(2, 3, 4, 5, 6, 7, 8)]
    rownames(output) <- NULL
  }
  output
}



#' Plot spatial intensity tree
#'
#' @param x A spatial intensity tree return by treerec function
#' @param ... additional arguments
#'
#' @return A plot of the spatial intensity tree 
#' and an \code{\link[spatstat]{im}} object.
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   minpts = 300,
#'   mtry = 1
#' )
#' plot(arbre)
plot.sptree <- function(x, ..., main) {
  # Handling case if no main title is given for the plot
  if (missing(main)) {
    main <- "Spatial Intensity Tree"
  }

  spatstat.geom::plot.im(x$im, main = main, ...)

  return(invisible(x$im))
}


#' Tree prediction
#'
#' @param object A spatial intensity tree return by treerec function
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @return A number .....
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   mtry = 1,
#'   minpts = 300
#' )
#' arbre$listcov <- list(
#'   grad = spatstat.data::bei.extra$grad,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' predict(object = arbre, newdata = c(100, 100))
predict.sptree <- function(object, newdata, ...) {
  # Test if the covariates are im object
  whichcovim <- unlist(lapply(
    object$listcov,
    spatstat.geom::is.im
  ))
  if (!all(whichcovim)) {
    stop("It appears that in predict.sptree, the covariables of the
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

  Zfun <- lapply(object$listcov, spatstat.geom::as.function.im)
  ptxy <- cbind(X$x, X$y)
  valsplits <- lapply(Zfun, FUN = function(j) {
    j(X)
  }) # FIXME What I am doing there ????

  temp <- lapply(object$tree, FUN = function(i) {
    c(
      i$status,
      i$split_var,
      i$split_val,
      i$intensity_pred,
      i$left_daughter,
      i$right_daughter
    )
  })

  treemat <- do.call(rbind, temp)

  output <- lapply(1:nrow(ptxy),
    FUN = function(i, ...) {
      node <- treemat[1, ]

      while (node[1] == 1) {
        if (valsplits[[node[2]]][i] < node[3]) {
          child <- node[5]
        } else {
          child <- node[6]
        }
        node <- treemat[child, ]
      }

      return(node[4])
    }
  )

  return(unlist(output))
}




#' Tree prediction (Save when attempting another version)
#'
#' @param object A spatial intensity tree return by treerec function
#' @param newdata a xy vector or a ppp object
#' @param ... Additional argument
#'
#' @return A number .....
#' @export
#'
#' @examples
#' arbre <- treerec(
#'   X = spatstat.data::bei,
#'   listcovariates = spatstat.data::bei.extra,
#'   mtry = 1,
#'   minpts = 300
#' )
#' arbre$listcov <- list(
#'   grad = spatstat.data::bei.extra$grad,
#'   elev = spatstat.data::bei.extra$elev
#' )
#' predict(object = arbre, newdata = c(100, 100))
predict.sptree_save <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else if (!spatstat.geom::is.ppp(newdata)) {
    X <- spatstat.geom::ppp(x = newdata[1], y = newdata[2], window = object$X$window)
  } else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n == 0) {
      return(NULL)
    }
    X <- newdata
  }


  Zfun <- lapply(object$listcov, spatstat.geom::as.function.im)
  ptxy <- cbind(X$x, X$y)
  valsplits <- lapply(Zfun, FUN = function(j) {
    j(X)
  })

  output <- sapply(1:nrow(ptxy),
    FUN = function(i, ...) {
      node <- object$tree[[1]]

      while (node$status == 1) {
        # child <- ifelse(Zfun[[node$split_var]](ptxy[i, 1], ptxy[i, 2]) < node$split_val,
        #                 node$left_daughter,
        #                 node$right_daughter
        # )

        child <- ifelse(valsplits[[node$split_var]][i] < node$split_val,
          node$left_daughter,
          node$right_daughter
        )

        node <- object$tree[[child]]
      }

      return(node$intensity_pred)
    }
  )

  return(output)
}
