# f <- function() {
#   spforest(
#     X = spatstat.data::bei,
#     listcovariates = spforest::beisoilres,
#     mtry = 1,
#     minpts = 200,
#     Ntree = 1,
#     parallel = F
#   )
#   future::plan("sequential")
# }
# f()
# A<-microbenchmark::microbenchmark(f(), times=1000)
# print(A)
# 
# # Unit: milliseconds
# # expr     min       lq     mean   median       uq      max neval
# # f() 150.689 369.3376 383.0576 377.7269 403.8627 471.3426   100
# 
# # Unit: milliseconds
# # expr      min       lq     mean   median       uq      max neval
# # f() 149.1744 369.7257 363.8065 374.9942 381.9154 538.5577  1000
# 
# #After modif which
# # Unit: milliseconds
# # expr     min       lq    mean   median      uq      max neval
# # f() 123.014 354.7156 342.918 369.4774 386.779 774.6513  1000
# 
# #After modif Vectorisation
# # Unit: milliseconds
# # expr      min       lq     mean   median       uq      max neval
# # f() 112.0398 128.4505 194.6082 136.1301 170.0241 490.8887  1000
# 
# # After removing splitsup <- !splitsub; splitsub[!splitsub] <- NA;  splitsup[!splitsup] <- NA
# # Unit: milliseconds
# # expr      min       lq     mean  median       uq    max neval
# # f() 117.1908 128.0696 175.8012 129.679 145.5564 776.25  1000
# 
# 
# # Profvis
# 
# library(profvis)
# 
# profvis(f())
# 
# # Improve speed score.split ---
# 
# fold <- function(n1, W1area) {
# 
#     val <- ifelse(W1area==0, 0,
#                   ifelse(n1 > 1, n1 * log((n1 - 1) / W1area), 0))
# 
#   return(val)
# }
# 
# 
# ffor <- function(n1, W1area) {
#   val <- rep(0, length(n1))
#   for (i in 1:length(n1)) {
#     val[i] <- ifelse(W1area[i] == 0, 0,
#                      ifelse(n1[i] > 1, n1[i] * log((n1[i] - 1) / W1area[i]), 0))
#   }
#   return(val)
# }
# 
# 
# fnew <- function(n1, W1area) {
# 
#   val <- rep(0, length(n1))
#   val[W1area>0 & n1>1] <- n1[W1area>0 & n1>1] * log((n1[W1area>0 & n1>1] - 1) / W1area[W1area>0 & n1>1])
# 
#   return(val)
# }
# 
# n1 = c(1,8,6,0,9,0)
# W1area = c(1,8,6,0,9,0)
# 
# fold(n1=n1, W1area=W1area)
# ffor(n1=n1, W1area=W1area)
# fnew(n1=n1, W1area=W1area)
# 
# microbenchmark::microbenchmark(fold(n1=n1, W1area=W1area),
#                                ffor(n1=n1, W1area=W1area),
#                                fnew(n1=n1, W1area=W1area), times=1e5)
