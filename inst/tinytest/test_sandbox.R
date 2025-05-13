

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
newcov <- beisoilres01[c(15, 5, 8, 10, 1, 4, 7, 2, 6, 9, 12, 14, 13, 3, 11)]
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


vipval <- rep(NA, 15)
for (j in 1:15) {
  vipval[j] <- mean(importance(A, id_cov = j, cores = 1, viptype = 3
  ))
}
barplot(vipval, names.arg = names(newcov))




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

foo <- NULL
for (i in 1:50) {
  foo[[i]] <- summary(A$trees[[i]])[,4] |> unlist() |> table()
}

do.call(rbind,foo) |> colMeans()


# plot(A)
#
# plot(Z)
