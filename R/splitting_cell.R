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
#' A <- splitcell2(
#'   X = X,
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
