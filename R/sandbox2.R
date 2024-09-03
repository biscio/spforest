# library(microbenchmark)
# library(profvis)
# A<-profvis(splitcell2(
#      X=X,
#      valpts = valpts,
#      vecval = vecval,
#      usecovariates = usecovariates,
#      dimcov = dimcov,
#      covrangex = covrangex,
#      covrangey = covrangey,
#      areapixel = areapixel,
#      threshold = 100
#    ))
# A
#
# microbenchmark(RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = lapply(beisoilres, FUN=function(i){
#     as.im(i, dimyx=c(10,20))
#   }),
#   Ntree = 10,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 1
# ))
#
#
# A<-RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = beisoilres,
#   Ntree = 10,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 1
# )
#
# plot(A)
#
# microbenchmark(RforestPP2(
#   X = spatstat.data::bei,
#   listcovariates = lapply(beisoilres, FUN=function(i){
#     as.im(i, dimyx=c(10,20))
#   }),
#   Ntree = 10,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 1
# ),
# RforestPP(
#   X = spatstat.data::bei,
#   listcovariates = lapply(beisoilres, FUN=function(i){
#     as.im(i, dimyx=c(10,20))
#   }),
#   Ntree = 10,
#   minpts = 100,
#   mtry = 1,
#   p = 0,
#   cores_trees = 1
# ))
# #
# # library(profvis)
# #
# # profvis(RforestPP2(
# #   X = spatstat.data::bei,
# #   listcovariates = lapply(beisoilres, FUN=function(i){
# #     as.im(i, dimyx=c(10,20))
# #   }),
# #   Ntree = 1,
# #   minpts = 100,
# #   mtry = 1,
# #   p = 0,
# #   cores_trees = 1
# # ))
# #
# # profvis(splitcell2(
# #   X=X,
# #   valpts = valpts,
# #   vecval = vecval,
# #   usecovariates = usecovariates,
# #   dimcov = dimcov,
# #   covrangex = covrangex,
# #   covrangey = covrangey,
# #   areapixel = areapixel,
# #   threshold = 100
# # ))
# #
# # # ftree <- function(){
# # #   treerec(
# # #     X = spatstat.data::bei,
# # #     listcovariates = Rsandbox::beisoilres,
# # #     mtry = 1,
# # #     minpts = 500
# # #   )
# # # }
# # #
# # # gtree <- function(){
# # #   intensitytree(
# # #     X = spatstat.data::bei,
# # #     listcovariates = Rsandbox::beisoilres,
# # #     mtry = 1,
# # #     minpts = 500
# # #   )
# # # }
# # #
# # # A<-ftree()
# # # B <- gtree()
# # # plot(A)
# # # plot(B)
# # # test<-A$im-B$im ### Small difference why ??
# # # max(abs(test))
# # #
# # # sapply(intensity_tree2, FUN = function(i) {
# # #   i$right_daughter
# # # })
# # # sapply(intensity_tree2, FUN = function(i) {
# # #   i$left_daughter
# # # })
# # # sapply(intensity_tree2, FUN = function(i) {
# # #   i$status
# # # })
# # # sapply(intensity_tree2, FUN = function(i) {
# # #   i$already_split
# # # })
# # # sapply(B$tree, FUN = function(i) {
# # #   i$intensity_pred
# # # })
# # # sapply(A$tree, FUN = function(i) {
# # #   i$intensity_pred
# # # })
# # #
# # # microbenchmark(ftree(), gtree(), times = 5)
# # #
# # # library(profvis)
# # # profvis(ftree())
# # # profvis(gtree())
# # #
# # #
# # #
# # # f <- function(){
# # #   RforestPP(
# # #     X = spatstat.data::bei,
# # #     listcovariates = Rsandbox::beisoilres,
# # #     Ntree = 10,
# # #     minpts = 500,
# # #     mtry = 1,
# # #     cores_trees = 5
# # #   )
# # # }
# # #
# # # g <- function(){
# # #   RforestPP2(
# # #     X = spatstat.data::bei,
# # #     listcovariates = Rsandbox::beisoilres,
# # #     Ntree = 10,
# # #     minpts = 500,
# # #     mtry = 1,
# # #     cores_trees = 5
# # #   )
# # # }
# # #
# # # library(microbenchmark)
# # # microbenchmark(f(), g(), times = 5)
# # #
