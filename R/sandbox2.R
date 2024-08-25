
################### 333

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
  stopifnot(spatstat.geom::is.ppp(X))
  stopifnot(spatstat.geom::is.im(listcovariates[[1]]))
  
  # Cas pathologique si pas de covariable - no split
  if (sum(usecovariates) == 0 | spatstat.geom::npoints(X) <= minpts) {
    whynot <- c("no covariate or less than minpts")
  }
  
  vecval <- lapply(listcovariates, FUN = function(i) {
    if (!spatstat.geom::is.im(i)) {
      stop("Elements of listcovar must be an im object")
    }
    c(as.matrix.im(i))
  })
  
  vecvalused <- vecval[usecovariates == 1]
  leftlvl0 <- lapply(1:length(split_scrs),
                     FUN = function(j) {
                       vecvalused[[j]] < split_scrs[[j]]
                     }
  )
  listcovar <- listcovariates[usecovariates == 1]
  
  valpts <- lapply(listcovariates[usecovariates == 1],
                   FUN = function(i) {
                     i[X]
                   }
  )
  
  scr_cov <- NULL
  id_best_scr <- NULL
  scr_parent <- score.pp(X = X, score = score)
  
  # Calculate the splitting scores of all covariates
  split_scrs <- lapply(vecvalused, FUN = function(j) {
    if (minsplitq != maxsplitq) {
      rand_q <- stats::runif(n = 1, min = minsplitq, max = maxsplitq)
      return(stats::quantile(j, probs = rand_q, na.rm = TRUE))
    } else {
      return(stats::median(j, na.rm = TRUE))
    }
  })
  
  # Remove names of columns due to quantile fct
  names(split_scrs) <- NULL
  
  leftlvl0 <- lapply(1:length(split_scrs),
                     FUN = function(j) {
                       vecvalused[[j]] < split_scrs[[j]]
                     }
  )
  
  areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
  
  # Determination of level sets for all covariables
  for (i in 1:sum(usecovariates)) {
    maxsplitq <- 0.5
    Wleft <- areapixel * sum(leftlvl0[[i]], na.rm = TRUE)
    Wright <- areapixel * sum(!leftlvl0[[i]], na.rm = TRUE)
    
    # Test if they are too small and computation of the score if not
    if (Wleft <= thres.cell |
        Wright <= thres.cell) {
      scr_cov[i] <- -Inf
    } else {
      n1 <- sum(valpts[[i]] < split_scrs[[i]], na.rm = TRUE)
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
    split_var <- which(usecovariates == 1)[id_best_scr]
    split_val <- split_scrs[[id_best_scr]]
    
    splitsub <- (vecval[[split_var]] < split_scrs[[id_best_scr]])
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
      PPleft = (valpts[[id_best_scr]] < split_val),
      Wleft = (vecval[[id_best_scr]] < split_val),
      split_var = split_var,
      split_val = split_val,
      sublevels = sublevels,
      suplevels = suplevels,
      improvement = newimp,
      scr_parent = scr_parent,
      whystop = NULL
    )
  }
  return(nodeChilds)
}

#### Test for splitcell2 ----
# f <- function() {
#   splitcell2(
#     X = spatstat.data::bei,
#     score = "lcv",
#     listcovariates = beisoilres,
#     usecovariates = rep(c(1, 0, 1), each = 5),
#     thres.cell = 100,
#     minpts = 10,
#     tol = Inf,
#     imp = NULL,
#     minsplitq = 0.5,
#     maxsplitq = 0.5
#   )
# }
# g <- function() {splitcell(
#   X = spatstat.data::bei,
#   score = "lcv",
#   listcovariates = beisoilres,
#   usecovariates = rep(c(1,0,1),each=5),
#   thres.cell = 100,
#   minpts = 10,
#   tol = Inf,
#   imp = NULL,
#   minsplitq = 0.5,
#   maxsplitq = 0.5
# )}
# A <- f()
# B <- g()
#
# library(microbenchmark)
#
# microbenchmark(f(),g())









# Problem with Y, it split differently between splitcell and splitcell2. 
# It also triggers an infinite loop in treerec2
# library(spatstat)
# 
# 
# X <- bei
# 
# XX <- bei.extra[[1]]
# pixelarea <- XX$xstep*XX$ystep
# 
# 
# valpts <- XX[X]
# vecval <- c(XX$v)
# 
# thres <- median(vecval)
# 
# ptsbelow <- (valpts<=thres)
# 
# n1 <- sum(ptsbelow)
# n2 <- npoints(X) - n1
# 
# npix1 <- sum(vecval <=thres)
# npix2 <- length(vecval) - npix1
# 
# 
# scr <- ifelse(n1 > 1, n1 * log((n1 - 1) / (npix1*pixelarea)), 0) +
#   ifelse(n2 > 1, n2 * log((n2 - 1) / (npix1*pixelarea)), 0)
# 
# Wleft <- (vecval <=thres)
# Wright <- (vecval > thres)
# 
# f=function(){
# XX[vecval <=thres]
# A<-matrix(vecval <=thres, nrow=XX$dim[1], ncol=XX$dim[2], byrow = F)
# im(A*XX$v, xrange=XX$xrange, yrange=XX$yrange)
# }
# 
# g <- function(){XX[XX<=thres,drop=F]}
# microbenchmark(f())
# microbenchmark(g())
# plot(XX[XX<=thres,drop=F])
# 
# 
# plot(XX[vecval <=thres, drop=F])
# 
# B<-splitcell(
#   X = X,
#   listcovariates = bei.extra,
#   usecovariates = c(1,1),
#   thres.cell = 25,
#   minpts = 10
# )
# 
# plot(B$PPleft)
# 
# vec <- rep(0,12)
# mat <- matrix(vec, nrow=3, ncol=4)
# mat[1,1] <- 1
# # plot(im(mat))
# 
# mat[1,2] <- 2
# mat[2,1] <- 3
# m <- im(mat)
# # plot(m)
# m$v
# c(m$v)
# library(microbenchmark)
# 
#  
# m$v
# 
