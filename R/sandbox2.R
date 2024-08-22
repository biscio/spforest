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
