forest <- RforestPP(X = spatstat.data::bei[seq(1, 3604, by = 50)],
                    listcovariates = list(
                        grad = spatstat.data::bei.extra$grad,
                        elev = spatstat.data::bei.extra$elev
                    ),
                    score = "lcv2",
                    p = 1,
                    Ntree = 3,
                    threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
                    cores_trees = 1,
                    mtry = 1 / 3,
                    tol = Inf,
                    minpts = 100,
                    minsplitq = 0.5,
                    maxsplitq = 0.5)

forest$trees[[2]]$namecov


# OOBscr.spforest <- function(forest, cores = 1) {
X <- forest$X # this is alway the root

# Put listcov back in the sptree object, required in predict.sptree
for (i in 1:length(forest$trees)) {
    forest$trees[[i]]$listcov <- forest$listcov
}

### TODO: handle special case length(forest$trees)==1

OOBscr <- parallel::mclapply(1:length(forest$trees), FUN = function(i) {

    # vector of same length as number of pts in X
    OOBpts <- (forest$pt_intree[[i]] != 1)
    # If all false then there is no pts in OOB sample and nothing to do
    if (all(!OOBpts)) {
        return(OOBval)
    }
    OOBval <- rep(NA, X$n)
    OOBval



    # OOB sample
    Xout <- X[OOBpts]
    Xout

    ### OOB prediction
    pts_pred_OOB <- predict.sptree(object = forest$trees[[i]],
                                   newdata = Xout)

    OOBval[OOBpts] <- pts_pred_OOB


    # OOB score
    return(OOBval)
}, mc.cores = cores)

OOBscr

logterm <- log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE))
if (all(is.na(logterm))) {
    output <- NA
} else {
output <- sum(logterm, na.rm = TRUE)
}
# Return the average error of all the trees
return(output)


simplify2array(OOBscr[[1]])
do.call(cbind, OOBscr)

A<-rowMeans(do.call(cbind, OOBscr), na.rm = TRUE)

sum(log(A), na.rm=T)


Reduce('+', OOBscr)
library(microbenchmark)

microbenchmark(do.call(cbind, OOBscr), times = 1e6)
microbenchmark(simplify2array(OOBscr), times = 1e6)


#### Predict function

if (missing(newdata) || is.null(newdata)) {
    X <- object$X
} else if (!spatstat.geom::is.ppp(newdata)) {
    X <- spatstat.geom::ppp(x = newdata[1], y = newdata[2], window = object$X$window)
} else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n==0) {
        return(NULL)
    }
    X <- newdata
}


Zfun <- lapply(object$listcov, spatstat.geom::as.function.im)
ptxy <- cbind(X$x, X$y)
valsplits <- lapply(Zfun, FUN = function(j) {j(X)})


output <- lapply(1:nrow(ptxy),
       FUN = function(i, ...) {
           node <- object$tree[[1]]

           while (node$status == 1) {

               child <- data.table::fifelse(valsplits[[node$split_var]][i] < node$split_val,
                                            node$left_daughter,
                                            node$right_daughter
               )

               node <- object$tree[[child]]
           }

           return(node$intensity_pred)
       }
)
unlist(output)

