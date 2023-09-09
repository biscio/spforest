library(profvis)

Xtest

target.points

profvis(tesstree(X=Xtest, lambda=100, 
                 target.points = target.points, test.connected = T))



a<-tessforest(X=Xtest, lambda=NULL, Ntree=100, mc.cores = 5, 
            at = NULL, test.connected = T)
plot(a)


library(spatstat)

X=rpoispp(100)

plot(X[1],pch=20)
points(X[2])
points(X[c(3,4)])

points(X[c(3,3)])

bootid <- rmultinom(n=1, size=npoints(X), prob=rep(1/npoints(X), npoints(X)))
sum(bootid)

X[c(2,2,2,3,3)]

a=NULL
for (i in 1:X$n) {
  a=c(a, sample.int(n=X$n, size = 1))
}
a

Y<-X[sample.int(n=X$n, size = X$n, replace=T)]
plot(X)
points(X[sample.int(n=X$n, size = X$n, replace=T)], pch=20,cex=0.5, col=4)
ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)


forest <- RforestPP(
  X = spatstat.data::bei,
  listcovariates = list(
    grad = spatstat.data::bei.extra$grad,
    elev = spatstat.data::bei.extra$elev
  ),
  score = "lcv2",
  p = 0,
  Ntree = 3,
  threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
  cores_trees = 1,
  mtry = 1 / 3,
  tol = Inf,
  minpts = 50,
  minsplitq = 0.5,
  maxsplitq = 0.5
)
plot(forest)
plot(forest$trees[[3]])
outputoob <- OOBscr.spforest(forest=forest, cores=1)
outputoob

X<-bei

X <- forest$X

X[ptintree] %>% unique()





a[(a %in% torm)]

unique(ptintree)

library(spatstat)
duplicated.ppp(X[ptintree])

