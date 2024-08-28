ftree <- function(){
  treerec(
    X = spatstat.data::bei,
    listcovariates = Rsandbox::beisoilres,
    mtry = 1,
    minpts = 500
  )
}

gtree <- function(){
  intensitytree(
    X = spatstat.data::bei,
    listcovariates = beisoilres,
    mtry = 1,
    minpts = 500
  )
}

A<-ftree()
B <- gtree()
plot(A)
plot(B)
test<-A$im-B$im ### Small difference why ??
max(abs(test))

microbenchmark(ftree(), gtree(), times = 5)

library(profvis)
profvis(ftree())
profvis(gtree())



f <- function(){
  RforestPP(
    X = spatstat.data::bei,
    listcovariates = Rsandbox::beisoilres,
    Ntree = 10,
    minpts = 500,
    mtry = 1,
    cores_trees = 1
  )
}

g <- function(){
  RforestPP2(
    X = spatstat.data::bei,
    listcovariates = Rsandbox::beisoilres,
    Ntree = 10,
    minpts = 500,
    mtry = 1, 
    cores_trees = 1
  )
}

library(microbenchmark)
microbenchmark(f(), g(), times = 2)

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

# 
# 
# sum(valpts[[res.split$split_var]]<=res.split$split_val)
# 
# sum(valpts[[res.split$split_var]]
#     <=res.split$split_val) 
# 
# (areapixel * sum(res.split$Wsub))
# 
# sum(!is.na(res.split$sublevels[[res.split$split_var]]))
# 
# A<-vecval[[res.split$split_var]] < res.split$split_val
# dimcov <- listcovariates[[1]]$dim
# covrangex <- listcovariates[[1]]$xrange
# covrangey <- listcovariates[[1]]$yrange
# temp <- im(matrix(A, nrow = dimcov[1], ncol = dimcov[2], byrow = F), 
#      xrange = covrangex, yrange = covrangey)
# 
# 
# 
# 
# 
# # splitcell2 - Code ---- 
# 
# splitcell2 <- function(X,
#                        score = "lcv",
#                        listcovariates = NULL,
#                        usecovariates = rep(1, length(listcovariates)),
#                        thres.cell = 25,
#                        minpts = 50,
#                        tol = Inf,
#                        imp = NULL,
#                        minsplitq = 0.5,
#                        maxsplitq = 0.5) {
#   whynot <- NULL
#   stopifnot(spatstat.geom::is.ppp(X))
#   stopifnot(spatstat.geom::is.im(listcovariates[[1]]))
#   
#   # Cas pathologique si pas de covariable - no split
#   if (sum(usecovariates) == 0 | spatstat.geom::npoints(X) <= minpts) {
#     whynot <- c("no covariate or less than minpts")
#   }
#   
#   vecval <- lapply(listcovariates, FUN = function(i) {
#     if (!spatstat.geom::is.im(i)) {
#       stop("Elements of listcovar must be an im object")
#     }
#     c(as.matrix.im(i))
#   })
#   
#   vecvalused <- vecval[usecovariates == 1]
#   leftlvl0 <- lapply(1:length(split_scrs),
#                      FUN = function(j) {
#                        vecvalused[[j]] < split_scrs[[j]]
#                      }
#   )
#   listcovar <- listcovariates[usecovariates == 1]
#   
#   valpts <- lapply(listcovariates[usecovariates == 1],
#                    FUN = function(i) {
#                      i[X]
#                    }
#   )
#   
#   scr_cov <- NULL
#   id_best_scr <- NULL
#   scr_parent <- score.pp(X = X, score = score)
#   
#   # Calculate the splitting scores of all covariates
#   split_scrs <- lapply(vecvalused, FUN = function(j) {
#     if (minsplitq != maxsplitq) {
#       rand_q <- stats::runif(n = 1, min = minsplitq, max = maxsplitq)
#       return(stats::quantile(j, probs = rand_q, na.rm = TRUE))
#     } else {
#       return(stats::median(j, na.rm = TRUE))
#     }
#   })
#   
#   # Remove names of columns due to quantile fct
#   names(split_scrs) <- NULL
#   
#   leftlvl0 <- lapply(1:length(split_scrs),
#                      FUN = function(j) {
#                        vecvalused[[j]] < split_scrs[[j]]
#                      }
#   )
#   
#   areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
#   
#   # Determination of level sets for all covariables
#   for (i in 1:sum(usecovariates)) {
#     maxsplitq <- 0.5
#     Wleft <- areapixel * sum(leftlvl0[[i]], na.rm = TRUE)
#     Wright <- areapixel * sum(!leftlvl0[[i]], na.rm = TRUE)
#     
#     # Test if they are too small and computation of the score if not
#     if (Wleft <= thres.cell |
#         Wright <= thres.cell) {
#       scr_cov[i] <- -Inf
#     } else {
#       n1 <- sum(valpts[[i]] < split_scrs[[i]], na.rm = TRUE)
#       scr_cov[i] <- score.split(
#         n1 = n1,
#         n2 = spatstat.geom::npoints(X) - n1,
#         W1area = Wleft,
#         W2area = Wright
#       )
#     }
#   }
#   
#   ## Go out if all the score are -Inf
#   if (all(is.infinite(scr_cov))) {
#     whynot <- c("All split scores are -Inf, cell too small")
#     # return(NULL)
#   }
#   
#   ## If there in the code
#   ## not all the score are -Inf so there is one best score
#   id_best_scr <- sort(scr_cov,
#                       index.return = T,
#                       decreasing = T
#   )$ix[1]
#   
#   # Look at the parents
#   # Look if there was any improvements on the score in the last splits
#   newimp <- append(
#     imp,
#     scr_cov[id_best_scr] > scr_parent
#   )
#   L_imp <- length(newimp)
#   
#   if (L_imp >= tol & all(!utils::tail(newimp, n = tol))) {
#     whynot <- c("no improvement since tol splits")
#     # return(NULL)
#   }
#   
#   # For split_var,
#   # I determine the index of the split var among all covariables
#   
#   if (!is.null(whynot)) {
#     return(whynot)
#   } else {
#     split_var <- which(usecovariates == 1)[id_best_scr]
#     split_val <- split_scrs[[id_best_scr]]
#     
#     splitsub <- (vecval[[split_var]] < split_scrs[[id_best_scr]])
#     splitsup <- !splitsub
#     
#     splitsub[!splitsub] <- NA
#     splitsup[!splitsup] <- NA
#     sublevels <- lapply(vecval, FUN = function(j) {
#       return(j * splitsub)
#     })
#     suplevels <- lapply(vecval, FUN = function(j) {
#       return(j * splitsup)
#     })
#     
#     nodeChilds <- list(
#       PPleft = (valpts[[id_best_scr]] < split_val),
#       Wleft = (vecval[[id_best_scr]] < split_val),
#       split_var = split_var,
#       split_val = split_val,
#       sublevels = sublevels,
#       suplevels = suplevels,
#       improvement = newimp,
#       scr_parent = scr_parent,
#       whystop = NULL
#     )
#   }
#   return(nodeChilds)
# }
# 
# #Test for splitcell2 ----
# # f <- function() {
# #   splitcell2(
# #     X = spatstat.data::bei,
# #     score = "lcv",
# #     listcovariates = beisoilres,
# #     usecovariates = rep(c(1, 0, 1), each = 5),
# #     thres.cell = 100,
# #     minpts = 10,
# #     tol = Inf,
# #     imp = NULL,
# #     minsplitq = 0.5,
# #     maxsplitq = 0.5
# #   )
# # }
# # g <- function() {splitcell(
# #   X = spatstat.data::bei,
# #   score = "lcv",
# #   listcovariates = beisoilres,
# #   usecovariates = rep(c(1,0,1),each=5),
# #   thres.cell = 100,
# #   minpts = 10,
# #   tol = Inf,
# #   imp = NULL,
# #   minsplitq = 0.5,
# #   maxsplitq = 0.5
# # )}
# # A <- f()
# # B <- g()
# #
# # library(microbenchmark)
# #
# # microbenchmark(f(),g())
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Problem with Y, it split differently between splitcell and splitcell2. 
# # It also triggers an infinite loop in treerec2
# # library(spatstat)
# # 
# # 
# # X <- bei
# # 
# # XX <- bei.extra[[1]]
# # pixelarea <- XX$xstep*XX$ystep
# # 
# # 
# # valpts <- XX[X]
# # vecval <- c(XX$v)
# # 
# # thres <- median(vecval)
# # 
# # ptsbelow <- (valpts<=thres)
# # 
# # n1 <- sum(ptsbelow)
# # n2 <- npoints(X) - n1
# # 
# # npix1 <- sum(vecval <=thres)
# # npix2 <- length(vecval) - npix1
# # 
# # 
# # scr <- ifelse(n1 > 1, n1 * log((n1 - 1) / (npix1*pixelarea)), 0) +
# #   ifelse(n2 > 1, n2 * log((n2 - 1) / (npix1*pixelarea)), 0)
# # 
# # Wleft <- (vecval <=thres)
# # Wright <- (vecval > thres)
# # 
# # f=function(){
# # XX[vecval <=thres]
# # A<-matrix(vecval <=thres, nrow=XX$dim[1], ncol=XX$dim[2], byrow = F)
# # im(A*XX$v, xrange=XX$xrange, yrange=XX$yrange)
# # }
# # 
# # g <- function(){XX[XX<=thres,drop=F]}
# # microbenchmark(f())
# # microbenchmark(g())
# # plot(XX[XX<=thres,drop=F])
# # 
# # 
# # plot(XX[vecval <=thres, drop=F])
# # 
# # B<-splitcell(
# #   X = X,
# #   listcovariates = bei.extra,
# #   usecovariates = c(1,1),
# #   thres.cell = 25,
# #   minpts = 10
# # )
# # 
# # plot(B$PPleft)
# # 
# # vec <- rep(0,12)
# # mat <- matrix(vec, nrow=3, ncol=4)
# # mat[1,1] <- 1
# # # plot(im(mat))
# # 
# # mat[1,2] <- 2
# # mat[2,1] <- 3
# # m <- im(mat)
# # # plot(m)
# # m$v
# # c(m$v)
# # library(microbenchmark)
# # 
# #  
# # m$v
# # 
# 
# 
# #### treerec2 code----
# treerec2 <- function(X,
#                      score = "lcv",
#                      threshold = spatstat.geom::area(X) / 1e4,
#                      listcovariates = NULL,
#                      mtry = 1,
#                      tol = Inf,
#                      minpts = spatstat.geom::npoints(X) / 10,
#                      minsplitq = 0.5,
#                      maxsplitq = 0.5,
#                      inforest = F) {
#   # Sanity checks
#   if (threshold <= 0 & minpts <= 0) {
#     stop("Either threshold or minpts must be strictly greater than 0.")
#   }
#   
#   areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
#   
#   # Initiating the root
#   
#   # vecval <- lapply(listcovariates, FUN = function(i) {
#   #   if (!spatstat.geom::is.im(i)) {
#   #     stop("Elements of listcovar must be an im object")
#   #   }
#   #   c(i$v)
#   # })
#   
#   root <- list(
#     nodeID = 1,
#     nodePP = rep(T, spatstat.geom::npoints(X)),
#     nodeCov = listcovariates,
#     left_daughter = NA,
#     right_daughter = NA,
#     split_var = NA,
#     split_val = NA,
#     status = 1,
#     intensity_pred = spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
#     scr_parent = NA,
#     improvement = NULL,
#     already_split = FALSE,
#     whystop = NULL
#   )
#   # Check if there is something to do on the initial point pattern.
#   # if (spatstat.geom::area.owin(X$window) <= threshold |
#   #   spatstat.geom::npoints(X) <= minpts) {
#   #   root$status <- 0
#   #   return(list(root))
#   # }
#   if (spatstat.geom::area.owin(X$window) <= threshold |
#       spatstat.geom::npoints(X) <= minpts) {
#     root$status <- 0
#     output <- list(
#       tree = list(root),
#       X = X,
#       namecov = names(listcovariates),
#       namelist = as.character(match.call()[4]),
#       im = as.im(spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
#                  W = X$window
#       )
#     )
#     class(output) <- "sptree"
#     return(output)
#   }
#   
#   intensity_tree2 <- list(root)
#   
#   # Initialise the while loop
#   k <- 0
#   knew <- 1
#   dimcov <- listcovariates[[1]]$dim
#   covrangex <- listcovariates[[1]]$xrange
#   covrangey <- listcovariates[[1]]$yrange
#   # Core computation of the tree
#   ### TODO
#   ### CHECK THAT WE CAN STOP !!!!
#   while (k != knew) {
#     k <- length(intensity_tree2)
#     
#     already_split_node <- sapply(intensity_tree2,
#                                  FUN = function(j) {
#                                    j$already_split
#                                  }
#     )
#     
#     for (i in (1:k)[!already_split_node]) {
#       ## Select randomly covariates
#       usedcov <- rand_covar(
#         listcovariates = listcovariates,
#         mtry = mtry
#       )
#       
#       # Split the cell, if the split is valid under the chosen parameters
#       res.split <- splitcell2(
#         X = X[intensity_tree2[[i]]$nodePP],
#         score = score,
#         listcovariates = intensity_tree2[[i]]$nodeCov,
#         usecovariates = usedcov,
#         thres.cell = threshold,
#         minpts = minpts,
#         tol = tol,
#         imp = intensity_tree2[[i]]$improvement,
#         minsplitq = minsplitq,
#         maxsplitq = maxsplitq
#       )
#       
#       if (is.character(res.split)) {
#         intensity_tree2[[i]]$status <- 0
#         intensity_tree2[[i]]$already_split <- TRUE
#         intensity_tree2[[i]]$whystop <- res.split
#       } else {
#         # Update the parent
#         intensity_tree2[[i]]$left_daughter <- knew + 1
#         intensity_tree2[[i]]$right_daughter <- knew + 2
#         intensity_tree2[[i]]$split_var <- res.split$split_var
#         intensity_tree2[[i]]$split_val <- res.split$split_val
#         intensity_tree2[[i]]$already_split <- TRUE
#         
#         newsublvl <- lapply(res.split$sublevels, FUN = function(jj) {
#           A <- matrix(jj, nrow = dimcov[1], ncol = dimcov[2], byrow = F)
#           im(A, xrange = covrangex, yrange = covrangey)
#         })
#         newsuplvl <- lapply(res.split$suplevels, FUN = function(jj) {
#           A <- matrix(jj, nrow = dimcov[1], ncol = dimcov[2], byrow = F)
#           im(A, xrange = covrangex, yrange = covrangey)
#         })
#         
#         # Define the children
#         childleft <- list(
#           nodeID = knew + 1,
#           nodePP = res.split$PPleft,
#           nodeCov = newsublvl,
#           left_daughter = NA,
#           right_daughter = NA,
#           split_var = NA,
#           split_val = NA,
#           status = 1,
#           intensity_pred = sum(res.split$PPleft) / (areapixel * sum(res.split$Wsub)),
#           scr_parent = res.split$scr_parent,
#           improvement = res.split$improvement,
#           already_split = FALSE,
#           whystop = NULL
#         )
#         childright <- list(
#           nodeID = knew + 2,
#           nodePP = !res.split$PPleft,
#           nodeCov = newsuplvl,
#           left_daughter = NA,
#           right_daughter = NA,
#           split_var = NA,
#           split_val = NA,
#           status = 1,
#           intensity_pred = sum(!res.split$PPleft) / (areapixel * sum(!res.split$Wsub)),
#           scr_parent = res.split$scr_parent,
#           improvement = res.split$improvement,
#           already_split = FALSE,
#           whystop = NULL
#         )
#         # append the children
#         intensity_tree2 <- append(
#           intensity_tree2,
#           list(childleft, childright)
#         )
#       }
#       knew <- length(intensity_tree2)
#     }
#     # intensity_tree2
#     sapply(intensity_tree2, FUN = function(i) {
#       i$right_daughter
#     })
#     sapply(intensity_tree2, FUN = function(i) {
#       i$left_daughter
#     })
#     sapply(intensity_tree2, FUN = function(i) {
#       sum(i$nodePP)
#     })
#     sapply(intensity_tree2, FUN = function(i) {
#       i$status
#     })
#     sapply(intensity_tree2, FUN = function(i) {
#       i$already_split
#     })
#   }
#   
#   # # Compute the image
#   # idterm <- sapply(intensity_tree2, FUN = function(i) {
#   #   i$status
#   # })
#   
#   ### TODO
#   ### HOW TO COMPUTE THE IMAGE NOW?
#   #
#   # listmask <- lapply(intensity_tree2[idterm == 0], FUN = function(i) {
#   #   i$nodePP$window$m # It is a mask usually but if no split is done then it is not.
#   # })
#   #
#   # intensity_prediction <- sapply(intensity_tree2[idterm == 0], FUN = function(i) {
#   #   i$intensity_pred
#   # })
#   #
#   # matim <- mapply("*", listmask, intensity_prediction, SIMPLIFY = FALSE)
#   #
#   # imoutput <- spatstat.geom::as.im(Reduce("+", matim),
#   #                                  W = intensity_tree2[[1]]$nodePP
#   # )
#   #
#   # # Remove all the intermediary PP in nodePP
#   # for (i in seq_along(intensity_tree2)) {
#   #   intensity_tree2[[i]]$nodePP <- NULL
#   #   intensity_tree2[[i]]$already_split <- NULL
#   #   # x$tree[[i]]$improvement <- NULL
#   # }
#   #
#   # if (inforest) {
#   #   output <- list(
#   #     tree = intensity_tree2,
#   #     X = X,
#   #     namecov = names(listcovariates),
#   #     namelist = as.character(match.call()[4]),
#   #     im = imoutput
#   #   )
#   # } else {
#   #   output <- list(
#   #     tree = intensity_tree2,
#   #     X = X,
#   #     namecov = names(listcovariates),
#   #     namelist = as.character(match.call()[4]),
#   #     listcov = listcovariates,
#   #     im = imoutput
#   #   )
#   # }
#   #
#   # class(output) <- "sptree" # For when I will define class
#   #
#   # return(output)
# }
# #
# #
# # treerec2(
# #     X = spatstat.data::bei,
# #     minpts = 50,
# #     threshold = 1000,
# #     listcovariates = bei.extra,
# #     mtry = 1,
# #   )
