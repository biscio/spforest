rm(list=ls())

beisoilres01 <- lapply(beisoilres,
                       FUN = function(x) {
                         (x - min(x)) / (max(x) - min(x))
                       }
)
Z <- beisoilres01$grad * beisoilres01$Mn + beisoilres01$Al
ztrue <- 1000 * Z / integral(Z)
X <- rpoispp(lambda = ztrue, nsim = 1)



newcov <- beisoilres01[sample(15)]
# Issue, invert 15 and 2. The position in the vector appears to impact the importance
newcov <- beisoilres01[c(2, 5, 8, 10, 1, 4, 7, 15, 6, 9, 12, 14, 13, 3, 11)]
A <- spforest(
  X = X,
  listcovariates = newcov,
  score = "lcv2",
  mtry = 0.8,
  minpts = 25,
  randmtry = TRUE,
  p = 0,
  Ntree = 50,
  cores = 10,
  threshold = 50
)

newcov2 <- beisoilres01[c(15, 5, 8, 10, 1, 4, 7, 2, 6, 9, 12, 14, 13, 3, 11)]
B <- spforest(
  X = X,
  listcovariates = newcov2,
  score = "lcv2",
  mtry = 0.8,
  minpts = 25,
  randmtry = TRUE,
  p = 0,
  Ntree = 50,
  cores = 10,
  threshold = 50
)


vipval <- rep(NA, 15)
vipval2 <- rep(NA, 15)
for (j in 1:15) {
  vipval[j] <- mean(importance(A, id_cov = j, viptype = 3
  ))
  vipval2[j] <- mean(importance(B, id_cov = j, viptype = 3
  ))
}
par(mfrow=c(2,1))
barplot(vipval, names.arg = names(newcov))
barplot(vipval2, names.arg = names(newcov2))


foo <- NULL
foo2 <- NULL
for (i in 1:50) {
  foo[[i]] <- summary(A$trees[[i]])[,4] |> unlist() |> table()
  foo2[[i]] <- summary(B$trees[[i]])[,4] |> unlist() |> table()
}

do.call(rbind,foo) |> colMeans()
do.call(rbind,foo2) |> colMeans()




vipval <- sapply(
  X = 1:length(newcov),
  FUN = function(i) {
    importance(A,
               id_cov = i, 
               cores = 1, 
               viptype = 3
    )
  }
)

barplot(colMeans(vipval), names.arg = names(newcov))

# plot(A)
#
# plot(Z)

smoothing(x=valpts[[1]]) |> plot(ylim=c(0.20,0.7))
for (i in 2:15) {
  lines(smoothing(x=valpts[[i]]), col=i)
}



valpts <- lapply(newcov2,
                 FUN = function(i) {
                   i[X]
                 }
)

sapply(valpts, mean)
