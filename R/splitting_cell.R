#' Splitting cell
#'
#' @param X The observed data point pattern.
#' @param score A score to choose among "lcv", "lcv2", "ent", "star", "ise", "isecv".
#' @param listcovariates A list with all the covariates used for the tree. The
#' covariates must be given as im object from the package spatstat.
#' @param usecovariates  Vector of 0,1 with covariate to choose
#' @param thres.cell Minimum size of a cell
#' @param minpts A positive integer. The minimum number of points allowed to try to split a cell.
#' @param tol A positive integer.
#' Look at the tol last gen to see if it was improved
#' @param imp A positive integer.
#' History of the improvement up to now. Used in random forest algorithm.
#' @param minsplitq A number in $(0,1)$. Determine the minimal allowed splitting value
#' of a covariable.
#' @param maxsplitq A number in $(0,1)$. Determine the maximal allowed splitting value
#' of a covariable.
#'
#' @return A list with PPleft = split.left[[id_best_scr]],
#' PPright = X[!leftlvl[[id_best_scr]]],
#' split_var = which(usecovariates == 1)[id_best_scr],
#' split_val = split_scrs[id_best_scr],
#' improvement = newimp,
#' scr_parent = scr_parent
#' @export
#'
#' @examples
#' splitcell(
#'   X = spatstat.data::bei,
#'   score = "lcv",
#'   listcovariates = list(
#'     grad = spatstat.data::bei.extra$grad,
#'     elev = spatstat.data::bei.extra$elev
#'   ),
#'   usecovariates = c(1, 1),
#'   thres.cell = 100,
#'   minpts = 10,
#'   tol = Inf,
#'   imp = NULL,
#'   minsplitq = 0.5,
#'   maxsplitq = 0.5
#' )
splitcell <- function(X,
                      score = "lcv",
                      listcovariates = NULL,
                      usecovariates = rep(1, length(listcovariates)),
                      thres.cell = 100,
                      minpts = 0,
                      tol = Inf,
                      imp = NULL,
                      minsplitq = 0.2,
                      maxsplitq = 0.8) {

  whynot=NULL
  # Check for incorrect input
  stopifnot(spatstat.geom::is.ppp(X))
  stopifnot(spatstat.geom::is.im(listcovariates[[1]]))

  # Cas pathologique si pas de covariable - no split
  if (sum(usecovariates) == 0 | spatstat.geom::npoints(X) <= minpts) { # normally I always have one
      whynot=c("no covariate or less than minpts")
      # return(NULL)
  }

  split.left <- NULL
  leftlvl <- NULL
  # split.right <- NULL
  scr_cov <- NULL
  id_best_scr <- NULL
  scr_parent <- score.pp(X = X, score = score)


  # Restriction of the used covariates to the windows of X

  ## Selection of the used covariates
  listcovar0 <- listcovariates[which(usecovariates == 1)]
  listcovar <- list()

  ## Restriction on the windows with optimised version for mask
  for (i in 1:length(listcovar0)) {
    if (spatstat.geom::is.mask(X$window)) {
      coval <- listcovar0[[i]]$v * (1 * X$window$m)
      coval[coval == 0] <- NA
      listcovar[[i]] <- listcovar0[[i]]
      listcovar[[i]]$v <- coval # todo : check same dimension
    } else {
      listcovar[[i]] <- listcovar0[[i]][X$window, drop = F]
    }
  }

  # Calculate the splitting scores of all covariates
  split_scrs <- sapply(listcovar, FUN = function(j) {
    if (minsplitq != maxsplitq) {
      rand_q <- stats::runif(n = 1, min = minsplitq, max = maxsplitq)
      return(stats::quantile(j, probs = rand_q))
    } else {
      return(stats::median(j))
    }
    # median(j)
    # rand_q <- stats::runif(n = 1, min = minsplitq, max = maxsplitq)
    # return(stats::quantile(j, probs = rand_q))
  })
  # Remove names of columns due to quantile fct
  names(split_scrs) <- NULL

  # Determination of level sets for all covariables
  for (i in 1:length(listcovar)) {
    if (!spatstat.geom::is.im(listcovar[[i]])) {
      stop("Elements of listcovar must be an im object")
    }

    leftlvl0 <- (listcovar[[i]] < split_scrs[i])

    areapixel <- leftlvl0$xstep * leftlvl0$ystep

    Wleft <- areapixel * sum(leftlvl0$v, na.rm = TRUE)
    Wright <- areapixel * sum(!leftlvl0$v, na.rm = TRUE)

    # Test if they are too small and computation of the score if not
    if (Wleft <= thres.cell |
      Wright <= thres.cell) {
      scr_cov[i] <- -Inf
      split.left[[i]] <- NULL
      # split.right[[i]] <- NULL
      leftlvl[[i]] <- NULL
    } else {
      split.left[[i]] <- X[leftlvl0]
      leftlvl[[i]] <- leftlvl0
      # split.right[[i]] <- X[!leftlvl]
      scr_cov[i] <- score.split(
        n1 = spatstat.geom::npoints(split.left[[i]]),
        n2 = spatstat.geom::npoints(X) - spatstat.geom::npoints(split.left[[i]]),
        W1area = Wleft,
        W2area = Wright
      )
    }
  }

  ## Go out if all the score are -Inf
  if (all(is.infinite(scr_cov))) {
      whynot=c("All split scores are -Inf, cell too small")
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

  # I can look tol steps in the past only if there has been enough split
  # faster if large vector (maybe)
  # if (L_imp >= tol) {
  #   if (all(!imp[seq(L_imp - tol + 1, L_imp)])) {
  #     return(output_nosplit)
  #   }
  # }
  # more readable
  if (L_imp >= tol & all(!utils::tail(newimp, n = tol))) {
      whynot=c("no improvement since tol splits")
    # return(NULL)
  }

  # For split_var, I determine the index of the split var among all covariables

  if (!is.null(whynot)) {
      return(whynot)
  } else {
  nodeChilds <- list(
    PPleft = split.left[[id_best_scr]],
    PPright = X[!leftlvl[[id_best_scr]]],
    # PPright = split.right[[id_best_scr]],
    split_var = which(usecovariates == 1)[id_best_scr],
    split_val = split_scrs[id_best_scr],
    improvement = newimp,
    scr_parent = scr_parent,
    whystop = NULL
  )
}
  return(nodeChilds)
}
