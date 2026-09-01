# library(spatstat)
# 
# set.seed(909)
# contcovar <- spatstat.random::rGRFexpo(
#   W = square(1),
#   mu = 0,
#   var = 0.5,
#   scale = 0.2,
#   nsim = 10
# ) |> solapply(exp)
# plot(contcovar)
# 
# Znum <- as.im(
#   matrix(sample(1:20, size = 9), nrow = 3),
#   W = square(1)
# )
# plot(Znum)
# 
# 
# Z <- matrix(factor(letters[sample(9)]), nrow = 3)
# allcov <- append(contcovar, list(Z))
# names(allcov) <- c(paste0("contcovar", 1:10), "Znum")
# W <- commonGrid(contcovar[[1]], Z)
# allcov[[11]] <- as.im(allcov[[11]], W = W)
# 
# trend <- allcov[[11]] * contcovar[[10]] + contcovar[[1]]
# trueint <- 500 * trend / mean(trend)
# plot(trueint)
# 
# X <- rpoispp(lambda = trueint)
# 
# 
# foo <- spforest(
#   X = X,
#   listcovariates = allcov,
#   Ntree = 50,
#   minpts = 20,
#   p = 0,
#   mtry = 7,
#   parallel = F
# )
# mean(abs(as.im(foo) - trueint))
# 
# test <- ppm(X ~ ., covariates = allcov)
# mean(abs(predict(test) - trueint))
# 
# par(mfrow = c(1, 3))
# plot(foo)
# plot(trueint)
# points(X, pch = 16, cex = 0.5, col = 0)
# plot(as.im(foo) - trueint)
# 
# imfcttoreal(Z = Znum, X = X)
# 
# # Other test ----
# set.seed(909)
# contcovar <- spatstat.random::rGRFexpo(
#   W = square(1),
#   mu = 0,
#   var = 0.5,
#   scale = 0.2,
#   nsim = 10
# ) |> solapply(exp)
# 
# Znum <- as.im(
#   matrix(c(20, 100, 0), nrow = 3),
#   W = square(1)
# )
# plot(Znum)
# 
# W <- commonGrid(contcovar[[1]], Znum)
# newZnum <- as.im(Znum, W = W)
# allcov <- append(contcovar, list(newZnum))
# names(allcov) <- c(paste0("contcovar", 1:10), "Znum")
# 
# trend <- newZnum / 30 + contcovar[[1]]
# trueint <- 500 * trend / mean(trend)
# plot(trueint)
# 
# X <- rpoispp(lambda = trueint)
# plot(X)
# 
# foo <- spforest(
#   X = X,
#   listcovariates = allcov,
#   Ntree = 150,
#   minpts = 40,
#   p = 0,
#   mtry = 7,
#   parallel = F
# )
# mean(abs(as.im(foo) - trueint))
# vip(foo)
# 
# test <- ppm(X ~ ., covariates = allcov)
# mean(abs(predict(test) - trueint))
# 
# par(mfrow = c(1, 3))
# plot(foo)
# plot(trueint)
# points(X, pch = 16, cex = 0.5, col = 0)
# plot(as.im(foo) - trueint)
# 
# 
# # more non linearity
# 
# W <- cut.im(beisoilres$Al,
#   breaks =
#     quantile(beisoilres$Al),
#   include.lowest = T
# )
# 
# levels(W) <- letters[1:4]
# plot(W)
# 
# newbei <- beisoilres
# newbei$Al <- W
# 
# facbei <- lapply(beisoilres, FUN = function(i) {
#   W <- cut.im(i,
#     breaks = quantile(i),
#     include.lowest = T
#   )
# 
#   levels(W) <- letters[1:4]
#   return(W)
# })
# 
# plot(as.anylist(facbei))
# 
# foo <- spforest(X = bei, 
#          listcovariates = facbei,
#          Ntree = 100,
#          mtry = 10,
#          minpts = 40,
#          parallel = F)
# 
# plot(log(as.im(foo)))
# 
# ## Other pp gneeration
# Zfac <- facbei$Al
# map <- c(a = 1, b = 15, c = 0, d= 10)
# 
# # convert factor pixels to numeric values via their labels
# vnum <- matrix(
#   map[as.character(Zfac$v)],
#   nrow = nrow(Zfac$v),
#   ncol = ncol(Zfac$v)
# )
# 
# # rebuild numeric image on the same grid
# Znum <- im(vnum,
#            xcol = Zfac$xcol,
#            yrow = Zfac$yrow,
#            unitname = unitname(Zfac))
# mean(Znum)
# mean(beisoilres$Cu)
# 
# newcov <- append(beisoilres[-3], list(Znum))
# names(newcov)[15] <- "Al"
# names(newcov)
# 
# covnorm <- lapply(newcov, FUN = function(i) {
#   (i - min(i)) / (max(i)-min(i))
# })
# 
# trueint <- 500 * (covnorm$Al+covnorm$grad) / mean(covnorm$Al+covnorm$grad) / area(Znum)
# plot(trueint)
# 
# X = rpoispp(lambda = trueint)
# plot(X)
# 
# 
# 
# foo <- spforest(X = X, 
#                 listcovariates = newcov,
#                 Ntree = 100,
#                 mtry = 10,
#                 minpts = 40,
#                 parallel = F)
# 
# vip(foo)
# 
# cm <- colourmap(
#   col = default.image.colours(),
#   range = c(0, max(range(as.im(foo))[2], range(trueint)[2]))
# )
# 
# par(mfrow = c(1, 2))
# plot(foo, col=cm)
# plot(trueint, col=cm)
# points(X, pch = 16, cex = 0.5, col = 0)
# # plot(abs(as.im(foo)-trueint))
