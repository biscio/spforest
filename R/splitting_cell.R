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
#' listcovariates <- spforest::beisoilres
#' X <- bei
#' valpts <- lapply(listcovariates,
#'   FUN = function(i) {
#'     i[X]
#'   }
#' )
#' vecval <- lapply(spforest::beisoilres, FUN = function(i) {
#'   if (!spatstat.geom::is.im(i)) {
#'     stop("Elements of listcovar must be an im object")
#'   }
#'   c(as.matrix.im(i))
#' })
#' usecovariates <- c(rep(1, 5), rep(0, 5), rep(1, 5))
#' areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
#' score <- "lcv"
#' A <- splitcell(
#'   X = X,
#'   valpts = valpts,
#'   vecval = vecval,
#'   usecovariates = usecovariates,
#'   areapixel = areapixel,
#'   threshold = 100
#' )
splitcell <- function(X,
                      valpts,
                      vecval,
                      usecovariates,
                      areapixel,
                      score = "lcv",
                      threshold = spatstat.geom::area(X) / 1e4) {
  whynot <- NULL
  vecvalused <- vecval[usecovariates == 1]
  valptsused <- valpts[usecovariates == 1]
  # listcovar <- listcovariates[usecovariates == 1]

  mediancov <- lapply(vecvalused, FUN = function(j) {
    stats::median(j, na.rm = TRUE)
  })

  sublvl <- lapply(1:length(mediancov),
    FUN = function(j) {
      vecvalused[[j]] < mediancov[[j]]
    }
  )

  # Determine level sets for all covariates
  scr_cov <- NULL
  scr_sub <- NULL
  scr_sup <- NULL
  for (i in 1:sum(usecovariates)) {
    Wsub <- areapixel * sum(sublvl[[i]], na.rm = TRUE)
    Wsup <- areapixel * sum(!sublvl[[i]], na.rm = TRUE)

    # Test if they are too small and computation of the score if not
    if (Wsub <= threshold | Wsup <= threshold) {
      scr_cov[i] <- -Inf
    } else {
      
      n1 <- sum(valptsused[[i]] < mediancov[[i]], na.rm = T)
      n2 <- sum(valptsused[[i]] >= mediancov[[i]], na.rm = T)
      
      if (i==1) {
        scr_parent <- score.pp(n0 = n1+n2,
                               W0area = Wsub + Wsup,
                               score = score)
      }
      
      scr_cov[i] <- score.split(
        n1 = n1,
        n2 = n2,
        W1area = Wsub,
        W2area = Wsup,
        score = score
      )
      
      scr_sub[i] <- score.pp(n0 = n1,
                             W0area = Wsub,
                             score = score)
      scr_sup[i] <- score.pp(n0 = n2,
                             W0area = Wsup,
                             score = score)
    }
  }
  
  ## Go out if all the score are -Inf
  if (all(is.infinite(scr_cov))) {
    whynot <- c("All split scores are -Inf, cell too small")
  }

  allscr <- rep(NA, length(vecval))
  allscr[usecovariates == 1] <- scr_cov
  
  id_best_scr <- sort(scr_cov,
    index.return = T,
    decreasing = T
  )$ix[1]


  if (!is.null(whynot)) {
    return(whynot)
  } else {
    # split_var <- which(usecovariates == 1)[id_best_scr]
    split_var <- which.max(allscr)
    split_val <- mediancov[[id_best_scr]]
    splitsub <- (vecval[[split_var]] < split_val)

    subvalpts <- (valpts[[split_var]] < split_val)
    nsub <- sum(subvalpts, na.rm = T)
    nsup <- sum(!subvalpts, na.rm = T)

    valptssub <- lapply(valpts, FUN = function(j) {
      ifelse(subvalpts,
        j, NA
      )
    })
    valptssup <- lapply(valpts, FUN = function(j) {
      ifelse(!subvalpts,
        j, NA
      )
    })

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
      whystop = NULL,
      scrsplit = max(allscr, na.rm = TRUE),
      scrdcr =  scr_parent - scr_sub[id_best_scr] - scr_sup[id_best_scr]
    )

    return(nodeChilds)
  }
}
