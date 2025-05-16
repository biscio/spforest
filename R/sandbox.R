
beisoilres01 <- lapply(beisoilres,
                       FUN = function(x) {
                         (x - min(x)) / (max(x) - min(x))
                       }
)
Z <- beisoilres01$grad * beisoilres01$Mn + beisoilres01$Al
ztrue <- 1000 * Z / integral(Z)
X <- rpoispp(lambda = ztrue, nsim = 1)


X = X
listcovariates = newcov
# listcovariates = newcov2
score = "lcv2"
mtry = 0.8
minpts = 25
randmtry = TRUE
p = 0
Ntree = 50
cores = 10
threshold = 50


# Start tessforest

nbcov <- length(listcovariates)
namescov <- names(listcovariates)

if (!do.call(spatstat.geom::compatible.im, unname(listcovariates))) {
  listcovariates <- do.call(
    spatstat.geom::harmonise.im,
    listcovariates
  )
  warning("The im objects in listcovariates have been
    harmonised with the function harmonise.im.")
}

# Used in all trees
areapixel <- listcovariates[[1]]$xstep * listcovariates[[1]]$ystep
vecval <- lapply(listcovariates, FUN = function(i) {
  if (!spatstat.geom::is.im(i)) {
    stop("Elements of listcovar must be an im object")
  }
  c(as.matrix.im(i))
})
dimcov <- listcovariates[[1]]$dim
covrangex <- listcovariates[[1]]$xrange
covrangey <- listcovariates[[1]]$yrange

# Call to tesscovtree: inloop or lapply for several trees

if (p == 0) { # bootstrap case, with replacement
  ptintree <- sample.int(n = X$n, size = X$n, replace = T)
  Xintree <- X[ptintree] # Should I use unique(X[ptintree]) ??
} else {
  ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)
  safety <- 1
  while (sum(ptintree) == 0 & safety <= 1e3) {
    ptintree <- stats::rbinom(n = X$n, size = 1, prob = p)
    safety <- safety + 1
  }
  if (safety > 1e3) {
    stop("Check your parameters, there is no point sampled in the trees")
  }
  Xintree <- X[ptintree == 1]
}

X = Xintree
inforest = T
# tree <- tesscovtree(
#   X = Xintree,
#   vecval = vecval,
#   areapixel = areapixel,
#   dimcov = dimcov,
#   covrangex = covrangex,
#   covrangey = covrangey,
#   listcovariates = listcovariates,
#   minpts = minpts,
#   mtry = mtry,
#   randmtry = randmtry,
#   score = score,
#   threshold = threshold,
#   inforest = T
# )

# Start tesscovtree

valpts <- lapply(listcovariates,
                 FUN = function(i) {
                   i[X]
                 }
)

root <- list(
  nodeID = 1,  nodeCov = vecval,   nodeValpts = valpts,
  left_daughter = NA,   right_daughter = NA,
  nX = spatstat.geom::npoints(X),
  split_var = NA,   split_val = NA,  status = 1,
  intensity_pred = spatstat.geom::npoints(X) / spatstat.geom::area(X$window),
  already_split = FALSE,  whystop = NULL,
  scrsplit = NA,  scrdcr = NA
)



intensity_tree <- list(root)

# Initialise the while loop
k <- 0
knew <- 1

while (k != knew) {
  k <- length(intensity_tree)
  
  already_split_node <- sapply(intensity_tree,
                               FUN = function(j) {
                                 j$already_split
                               }
  )
  
  for (i in (1:k)[!already_split_node]) {

    if (intensity_tree[[i]]$nX <= minpts) {
      res.split <- "Not enough points to attempt to split"
    } else {
      # Split the cell, if the split is valid under the chosen parameters
     
      foo <- NULL
      for (j in 1:1000) {
        res.split <- splitcell(
          X = X,
          valpts = intensity_tree[[i]]$nodeValpts,
          vecval = intensity_tree[[i]]$nodeCov,
          usecovariates = rand_covar(
            listcovariates = listcovariates,
            mtry = mtry,
            randmtry = randmtry
          ),
          areapixel = areapixel,
          score = score,
          threshold = threshold
        )
       foo[j] <- names(listcovariates)[res.split$split_var]
      }
      table(foo)
      
    }
    
    if (is.character(res.split)) {
      intensity_tree[[i]]$status <- 0
      intensity_tree[[i]]$already_split <- TRUE
      intensity_tree[[i]]$whystop <- res.split
      intensity_tree[[i]]$scrsplit <- NA 
    } else {
      # Update the parent
      intensity_tree[[i]]$left_daughter <- knew + 1
      intensity_tree[[i]]$right_daughter <- knew + 2
      intensity_tree[[i]]$split_var <- res.split$split_var
      intensity_tree[[i]]$split_val <- res.split$split_val
      intensity_tree[[i]]$already_split <- TRUE
      intensity_tree[[i]]$scrsplit <- res.split$scrsplit 
      intensity_tree[[i]]$scrdcr <- res.split$scrdcr
      
      # Define the children
      areasub <- (areapixel * sum(!is.na(res.split$sublevels[[res.split$split_var]])))
      areasup <- (areapixel * sum(!is.na(res.split$suplevels[[res.split$split_var]])))
      
      childleft <- list(
        nodeID = knew + 1,
        nodeCov = res.split$sublevels,
        nodeValpts = res.split$valptssub,
        nX = res.split$nsub,
        left_daughter = NA,
        right_daughter = NA,
        split_var = NA,
        split_val = NA,
        status = 1,
        intensity_pred = res.split$nsub / areasub,
        already_split = FALSE,
        whystop = NULL,
        scrsplit = NA,
        scrdcr = NA
      )
      
      childright <- list(
        nodeID = knew + 2,
        nodeCov = res.split$suplevels,
        nodeValpts = res.split$valptssup,
        nX = res.split$nsup,
        left_daughter = NA,
        right_daughter = NA,
        split_var = NA,
        split_val = NA,
        status = 1,
        intensity_pred = res.split$nsup / areasup,
        already_split = FALSE,
        whystop = NULL,
        scrsplit = NA,
        scrdcr = NA
      )
      # append the children
      intensity_tree <- append(
        intensity_tree,
        list(childleft, childright)
      )
    }
    knew <- length(intensity_tree)
  }
}