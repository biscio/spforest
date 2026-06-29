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

  scr_parent <- score.pp(
    n0 = length(stats::na.omit(valptsused[[1]])),
    W0area = areapixel * sum(!is.na(vecvalused[[1]])),
    score = score
  )

  ### Pre-compute areas for all used covariates
  nsub_pix <- vapply(sublvl, function(s) sum(s, na.rm = TRUE), numeric(1L))
  nsup_pix <- vapply(sublvl, function(s) sum(!s, na.rm = TRUE), numeric(1L))
  Wsub_all <- areapixel * nsub_pix
  Wsup_all <- areapixel * nsup_pix

  too_small <- (Wsub_all <= threshold) | (Wsup_all <= threshold)

  ### Point counts per covariate
  n1_all <- vapply(
    seq_along(valptsused),
    function(i) sum(valptsused[[i]] < mediancov[[i]], na.rm = TRUE), numeric(1L)
  )
  n2_all <- vapply(
    seq_along(valptsused),
    function(i) sum(valptsused[[i]] >= mediancov[[i]], na.rm = TRUE), numeric(1L)
  )

  ### Score all covariates at once (score.split / score.pp are already vectorised)
  scr_cov <- score.split(n1_all, n2_all, Wsub_all, Wsup_all, score)
  scr_sub <- score.pp(n1_all, Wsub_all, score)
  scr_sup <- score.pp(n2_all, Wsup_all, score)

  ### Apply threshold mask
  scr_cov[too_small] <- -Inf
  scr_sub[too_small] <- NA
  scr_sup[too_small] <- NA

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

    idx_sub <- which(subvalpts)
    valptssub <- lapply(valpts, function(j) {
      j[-idx_sub] <- NA
      j
    })
    valptssup <- lapply(valpts, function(j) {
      j[idx_sub] <- NA
      j
    })
    # valptssub <- lapply(valpts, FUN = function(j) {
    #   A <- rep(NA, length(j))
    #   A[which(subvalpts)] <- j[which(subvalpts)]
    #   return(A)
    # })
    # valptssup <- lapply(valpts, FUN = function(j) {
    #   B <- j
    #   B[which(subvalpts)] <- NA
    #   return(B)
    # })

    # Precompute integer index vectors once — reused across all covariates
    idx_not_sub <- which(!splitsub | is.na(splitsub))
    idx_not_sup <- which(splitsub  | is.na(!splitsub))
    
    sublevels <- lapply(vecval, function(j) { j[idx_not_sub] <- NA_real_; j })
    suplevels <- lapply(vecval, function(j) { j[idx_not_sup] <- NA_real_; j })
    
    ## Older code start
    # splitsup <- !splitsub
    # splitsub[!splitsub] <- NA
    # splitsup[!splitsup] <- NA
    # sublevels <- lapply(vecval, FUN = function(j) {
    #   return(j * splitsub)
    # })
    # suplevels <- lapply(vecval, FUN = function(j) {
    #   return(j * splitsup)
    # })
    ## Older code end 
    
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
      scrdcr = scr_sub[id_best_scr] + scr_sup[id_best_scr] - scr_parent
    )

    return(nodeChilds)
  }
}





#' Test palm split
#'
#' @param X 
#' @param valpts 
#' @param vecval 
#' @param listcovariates 
#' @param usecovariates 
#' @param areapixel 
#' @param clustmodel 
#' @param threshold 
#'
#' @returns 
#' @export
#'
#' @examples
splitcellpalm <- function(X,
                          valpts,
                          vecval,
                          listcovariates,
                          usecovariates,
                          areapixel,
                          clustmodel = "LGCP",
                          threshold = spatstat.geom::area(X) / 1e4) {
  whynot <- NULL
  vecvalused <- vecval[usecovariates == 1]
  valptsused <- valpts[usecovariates == 1]
  listused <- listcovariates[usecovariates == 1]
  
  mediancov <- lapply(vecvalused, FUN = function(j) {
    stats::median(j, na.rm = TRUE)
  })
  
  sublvl <- lapply(1:length(mediancov),
                   FUN = function(j) {
                     solutionset(listused[[j]] < mediancov[[j]])
                   }
  )
  suplvl <- lapply(1:length(mediancov),
                   FUN = function(j) {
                     solutionset(listused[[j]] >= mediancov[[j]])
                   }
  )
  
  # Compute palm score
  scr_parent <- as.numeric(logLik.kppm(
    kppm.ppp(X, ~., clusters = clustmodel, method = "palm", data = listused)
  ))
  
  
  scr_cov <- sapply(1:length(listused), FUN = function(j) {
    as.numeric(logLik.kppm(
      kppm.ppp(X, as.formula(paste0("~", names(listused)[j])), clusters = clustmodel, method = "palm", data = listused)
    ))
  })
  
  scr_sup <- sapply(1:length(listused), FUN = function(j) {
    as.numeric(logLik.kppm(
      kppm.ppp(X[suplvl[[j]]], as.formula(paste0("~", names(listused)[j])), clusters = clustmodel, method = "palm", data = listused)
    ))
  })
  scr_sub <- sapply(1:length(listused), FUN = function(j) {
    as.numeric(logLik.kppm(
      kppm.ppp(X[sublvl[[j]]], as.formula(paste0("~", names(listused)[j])), clusters = clustmodel, method = "palm", data = listused)
    ))
  })
  
  ### Pre-compute areas for all used covariates
  nsub_pix <- vapply(sublvl, function(s) area.owin(s), numeric(1L))
  nsup_pix <- vapply(suplvl, function(s) area.owin(s), numeric(1L))
  Wsub_all <- areapixel * nsub_pix
  Wsup_all <- areapixel * nsup_pix
  
  too_small <- (Wsub_all <= threshold) | (Wsup_all <= threshold)
  
  ### Apply threshold mask
  scr_cov[too_small] <- -Inf
  scr_sub[too_small] <- NA
  scr_sup[too_small] <- NA
  
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
    
    idx_sub <- which(subvalpts)
    valptssub <- lapply(valpts, function(j) {
      j[-idx_sub] <- NA
      j
    })
    valptssup <- lapply(valpts, function(j) {
      j[idx_sub] <- NA
      j
    })
    
    # Precompute integer index vectors once — reused across all covariates
    idx_not_sub <- which(!splitsub | is.na(splitsub))
    idx_not_sup <- which(splitsub | is.na(!splitsub))
    
    sublevels <- lapply(vecval, function(j) {
      j[idx_not_sub] <- NA_real_
      j
    })
    suplevels <- lapply(vecval, function(j) {
      j[idx_not_sup] <- NA_real_
      j
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
      scrdcr = scr_sub[id_best_scr] + scr_sup[id_best_scr] - scr_parent
    )
    
    return(nodeChilds)
  }
}



