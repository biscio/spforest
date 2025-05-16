beisoilres01 <- lapply(beisoilres,
  FUN = function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
)
Z <- beisoilres01$grad * beisoilres01$Mn + beisoilres01$Al
ztrue <- 1000 * Z / integral(Z)
X <- rpoispp(lambda = ztrue, nsim = 1)

newcov  <- beisoilres01[c(2, 5, 8, 10, 1, 4, 7, 15, 6, 9, 12, 14, 13, 3, 11)]
newcov2 <- beisoilres01[c(15, 5, 8, 10, 1, 4, 7, 2, 6, 9, 12, 14, 13, 3, 11)]
score <- "lcv2"
mtry <- 0.8
minpts <- 25
randmtry <- TRUE
p <- 0
Ntree <- 50
cores <- 10
threshold <- 50

ptintree <- sample.int(n = X$n, size = X$n, replace = T)
X0  <- X[ptintree] # Should I use unique(X[ptintree]) ??

inforest <- T

vecval <- lapply(newcov, FUN = function(i) {
  if (!spatstat.geom::is.im(i)) {
    stop("Elements of listcovar must be an im object")
  }
  c(as.matrix.im(i))
})
vecval2 <- lapply(newcov2, FUN = function(i) {
  if (!spatstat.geom::is.im(i)) {
    stop("Elements of listcovar must be an im object")
  }
  c(as.matrix.im(i))
})

dimcov <- newcov[[1]]$dim

# Start tesscovtree

valpts <- lapply(newcov,
  FUN = function(i) {
    i[X0]
  }
)
valpts2 <- lapply(newcov2,
                 FUN = function(i) {
                   i[X0]
                 }
)


foo <- NULL
foo2 <- NULL
for (j in 1:10000) {
  res.split <- splitcell(
    X = X0,
    valpts = valpts,
    vecval = vecval,
    usecovariates = rand_covar(
      listcovariates = newcov,
      mtry = mtry,
      randmtry = randmtry
    ),
    areapixel = areapixel,
    score = score,
    threshold = threshold
  )
  foo[j] <- names(newcov)[res.split$split_var]
  
  res.split2 <- splitcell(
    X = X0,
    valpts = valpts2,
    vecval = vecval2,
    usecovariates = rand_covar(
      listcovariates = newcov2,
      mtry = mtry,
      randmtry = randmtry
    ),
    areapixel = areapixel,
    score = score,
    threshold = threshold
  )
  foo2[j] <- names(newcov2)[res.split2$split_var]
  
  if(j%%100==0){print(j)}
}

table(foo)
table(foo2)

# There is a problem there, the name should appear the same number of time on average. 


foo3 <- NULL
for (i in 1:10000) {
  foo3[[i]] <- rand_covar(
    listcovariates = newcov,
    mtry = 0.8,
    randmtry = T
  )
}
foo3<-do.call(rbind, foo3)
colMeans(foo3)

