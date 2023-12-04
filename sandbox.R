rm(list = ls())
library(rworldmap)
library(rworldxtra)
library(raster)
library(spatstat)
library(maptools) # Soon deprecated
library(sf)
library(sp)
library(parallel)
library(Rsandbox)
spatstat.options(npixel = 512)

# Dans draft de fred
# {Trout_UK2}
# {UK_Trout_ppl}
# {UK_Trout_tree}
# {UK_Trout_RF}
# {UK_Trout_ppl_clip}
# {UK_Trout_RF_clip}



library(terra)
library(geodata)
library(Rsandbox)

uk <- gadm(
  country = "GBR", level = 0, resolution = 2,
  path = "maps/"
)
irl <- gadm(
  country = "IRL", level = 0, resolution = 1,
  path = "maps/"
)
spatstat.options(npixel = 1024)

split(uk)

uksp <- as(uk, "Spatial")
p1 <- slot(uksp, "polygons")
v1 <- lapply(p1, function(z) {
  SpatialPolygons(list(z))
})
UK <- spatstat.geom::as.owin(v1[[1]])

plot(UK)

UK$bdry %>% length()
A <- connected(UK, eps = 0.005)
tmp <- spatstat.geom::tiles(as.tess(A))

allarea <- sapply(tmp, area.owin)
thres <- quantile(allarea, probs = 0.75)

newtess <- tmp[allarea > thres]

listW <- lapply(newtess, as.owin)

B <- do.call(union.owin, listW)
polyW <- as.polygonal(B)
plot(B[roi.islands])
plot(polyW[roi.islands])
as.owin(newtess)
Wuk <- polyW


plot(A)


subw <- lapply(UK$bdry, extr_owin)

a <- sapply(subw, area.owin)

thres <- quantile(a, probs = 0.5)

plot(UK)

for (i in 1:sum(a > thres)) {
  temp <- subw[a > thres]
  plot(temp[[i]], col = 2, add = T)
}

# Get the map of UK ----

## high create an error on my laptop
res_world_map <- "low" # Other solution high
if (!(res_world_map %in% c("low", "high"))) {
  stop("res_world_map must be set to low or high as a string.")
}


world <- rworldmap::getMap(resolution = res_world_map)
## clip in long lag coordinate a region
clipper_europe <- as(raster::extent(-10, 32, 30, 71.9), "SpatialPolygons")
world_clip <- raster::intersect(world, clipper_europe)

# Conversion to owin object ----

## A nice code from Adrian on stackoverflow
p <- slot(world_clip, "polygons")
v <- lapply(p, function(z) {
  SpatialPolygons(list(z))
})
winlist <- lapply(v, spatstat.geom::as.owin)


# Union of all the countries
europe <- do.call(spatstat.geom::union.owin, winlist)

# Loading the dataset and convert it to ppp ----
brown_trout <- trout

## Set a more convenient name for columns
names(brown_trout)[4] <- "lat"
names(brown_trout)[5] <- "lon"

## Set the point pattern as ppp
fulldata <- ppp(
  x = brown_trout$lon,
  y = brown_trout$lat,
  window = europe
)


# Delimiting the data on the UK ----

## If using the low resolution of the worldmap
if (res_world_map == "low") {
  UK <- union.owin(winlist[[20]], winlist[[26]]) # low res
} else {
  UK <- union.owin(winlist[[20]], winlist[[26]], winlist[[27]])
}

plot(UK)
subdata <- subset(fulldata, UK)
Window(subdata) <- UK

xmin <- -8
xmax <- -5
ymin <- 56.5
ymax <- 58.7
roi.islands <- owin(c(xmin, xmax), c(ymin, ymax))


# Plot I: The data on the all UK
par(mar = c(0, 0, 0, 0))
plot(subdata, pch = 20, main = "", box = TRUE)
plot(roi.islands, add = T, lty = 2, border = 2, lwd = 2)
text(-8.8, 60, "Data", cex = 1.5)
# dev.print(pdf, "Trout_UK2.pdf", width = 5, height = 5)

# Setting the parameters for kernel intensity estim ----
spatstat.options(npixel = 512)
res.ppl <- density(subdata, sigma = bw.ppl) # bandwidth ppl

co <- colourmap(getSpatstatVariable("DefaultImageColours")[1:256],
  range = range(res.ppl)
)

# Plot II: Kernel intensity estimate on all UK
par(mar = c(0, 0, 0, 1))
plot(res.ppl,
  main = "",
  ribbon = TRUE,
  ribsep = 0.05,
  workaround = FALSE,
  ribargs = list(cex.axis = 1.5)
)
text(-8.8, 60, "Kernel", cex = 2)
# dev.print(pdf, "UK_Trout_ppl.pdf", width = 8, height = 7)

# Plot III: ZOOM of kernel intensity estimate on all UK
plot(res.ppl,
  clipwin = roi.islands,
  col = co,
  main = "",
  ribbon = TRUE,
  ribsep = 0.05,
  ribargs = list(cex.axis = 1.5)
)
text(-7.4, 58.5, "Kernel (zoom)", cex = 2.4)
# dev.print(pdf, "UK_Trout_ppl_clip.pdf", width = 10, height = 7)


# Set parameters for tree intensity estimate ----
X <- subdata
lambda <- 409
wind <- Window(X)
x.image <- res.ppl$xcol
y.image <- res.ppl$yrow
allpoints <- as.ppp(
  expand.grid(x.image, y.image),
  owin(wind$xrange, wind$yrange)
)
target.points <- subset(allpoints, wind)


spatstat.options(npixel = 256)
# Tree intensity estimate computation
set.seed(1)
res.tree <- tesstree(
  X = subdata,
  lambda = lambda,
  target.points = target.points,
  test.connected = TRUE
)

marks(allpoints) <- NA
a <- inside.owin(allpoints, w = wind)
marks(allpoints)[a] <- res.tree[[1]]
tree.im <- as.im(
  t(matrix(
    marks(allpoints),
    length(x.image), length(x.image)
  )),
  W = owin(wind$xrange, wind$yrange)
)
plot(tree.im,
  col = co,
  main = "",
  ribbon = TRUE,
  ribsep = 0.05,
  ribargs = list(cex.axis = 1.5)
)
plot(res.tree[[2]], add = T, border = "white")
text(-9, 60, "Tree", cex = 2)
# dev.print(pdf, "UK_Trout_Tree.pdf", width = 8, height = 7)

# for testing: zoom/clip of the tree estimate on the island
plot(tree.im,
  clipwin = roi.islands,
  col = co,
  main = "",
  ribbon = TRUE,
  ribsep = 0.05,
  ribargs = list(cex.axis = 1.5)
)
plot(res.tree[[2]], add = T, border = "black", lwd = 4)



# RF (full + clip)
allpoints <- as.ppp(
  expand.grid(x.image, y.image),
  owin(wind$xrange, wind$yrange)
)
target.points <- subset(allpoints, wind)

res <- tessforest(subdata,
  lambda = 409,
  mc.cores = 2,
  N = 100,
  test.connected = TRUE,
  at = target.points
) # lambda=409, consistent with lambda chosen on the ROI

marks(allpoints) <- NA
a <- inside.owin(allpoints, w = wind)
marks(allpoints)[a] <- res
N <- length(x.image)

res.im <- as.im(t(matrix(marks(allpoints), N, N)),
  W = owin(wind$xrange, wind$yrange)
)


plot(res.im,
  col = co, main = "", ribbon = TRUE,
  ribsep = 0.05, ribargs = list(cex.axis = 1.5)
)
text(-9, 60, "RF", cex = 2)
# dev.print(pdf, "UK_Trout_RF.pdf", width = 8, height = 7)

plot(res.im,
  clipwin = roi.islands, col = co, main = "",
  ribbon = TRUE, ribsep = 0.05, ribargs = list(cex.axis = 1.5)
)
points(subdata[roi.islands], col = 1, pch = 20, cex = 2)
text(-7.5, 58.5, "RF (zoom)", cex = 2.4)
dev.print(pdf, "test.pdf", width = 10, height = 7)


############# 3
#############
#############
#############


test <- lapply(1:nn, FUN = function(i) {
  spatstat.geom::as.polygonal(tile.tmp[[i]])$bdry
})

concomp <- sapply(test, length)



plot(tmp)
for (i in which(concomp > 1)) {
  plot(as.polygonal(tile.tmp[[i]]), add = T, col = 2)
}

test[[15]]

dim(gDistance(spts,
  firstSpatialPoly,
  byid = TRUE
))


apply(apply(gDistance(spts,
  firstSpatialPoly,
  byid = TRUE
), 2, min), 2, which.min)

apply(sf::st_distance(
  st_as_sf(spts),
  st_as_sf(firstSpatialPoly),
  by_element = F
), 1, min)



if (testing) {
  argmtry <- c(0.5, 0.7)
  argminpts <- c(70, 80)
  Ntree <- 5
}
argmtry <- seq(0.1, 1, by = 0.1)
argminpts <- c(seq(10, 60, by = 5))
Ntree <- 300

argu <- expand.grid(argmtry, argminpts)
dim(argu)


f <- function(object, newdata, ...) {
  # Test if the covariates are im object
  whichcovim <- unlist(lapply(object$listcov, spatstat.geom::is.im))
  if (!all(whichcovim)) {
    stop("It appears that in predict.sptree, the covariables of the
             tree are not spatstat im objects")
  }

  # Handles the newdata to be in the correct form
  if (missing(newdata) || is.null(newdata)) {
    X <- object$X
  } else if (!spatstat.geom::is.ppp(newdata)) {
    X <- spatstat.geom::ppp(x = newdata[1], y = newdata[2], window = object$X$window)
  } else if (spatstat.geom::is.ppp(newdata)) {
    if (newdata$n == 0) {
      return(NULL)
    }
    X <- newdata
  }

  Zfun <- lapply(object$listcov, spatstat.geom::as.function.im)
  ptxy <- cbind(X$x, X$y)
  valsplits <- lapply(Zfun, FUN = function(j) {
    j(X)
  }) # FIXME What I am doing there ????

  A<-lapply(object$tree, FUN=function(i) {
    c(i$status, 
      i$split_var, 
      i$split_val, 
      i$intensity_pred,
      i$left_daughter, 
      i$right_daughter)}
  )
  
  B <- do.call(rbind, A)
  
  output <- lapply(1:nrow(ptxy),
    FUN = function(i, ...) {
      node <- B[1,]

      while (node[1] == 1) {
        if (valsplits[[node[2]]][i] < node[3]) {
          child <- node[5]
        } else {
          child <- node[6]
        }
        # child <- data.table::fifelse(
        #   valsplits[[node[2]]][i] < node[3],
        #   node[5],
        #   node[6]
        # )

        node <- B[child,]
      }

      return(node[4])
    }
  )

  return(unlist(output))
}


microbenchmark(f(object = arbre, 
                 newdata = spatstat.data::bei), times=500)
microbenchmark(predict.sptree(object = arbre, 
                              newdata = spatstat.data::bei, times=500))
profvis(f(object = arbre, 
          newdata =testX))

A<-f(object = arbre, newdata = spatstat.data::bei)

B<-predict.sptree(object = arbre, newdata = spatstat.data::bei)

max(A-B)

names(B) <- c("status", "split_var", 
              "left_daughter", "right_daughter" )



microbenchmark(do.call(cbind, A))
microbenchmark(do.call(rbind, A))

library(microbenchmark)
cbind(arbre$tree)




library(profvis)
arbre <- treerec(
  X = spatstat.data::bei,
  threshold = 1000,
  score = "lcv2",
  listcovariates = beisoilres,
  mtry = 1,
  tol = Inf,
  minpts = 50
)
arbre$listcov <- beisoilres
profvis(predict.sptree(object = arbre, 
                       newdata = rshift.ppp(spatstat.data::bei, radius=20)))

profvis(f(object = arbre,
          newdata = rshift.ppp(spatstat.data::bei, radius=20)))

testX<-rshift.ppp(spatstat.data::bei, radius=20)

A<-f(object = arbre, 
     newdata = testX)

B<-predict.sptree(object = arbre, 
                  newdata = testX)

max(abs(A-B))


g <- function(forest, cores = 1) {
  X <- forest$X # this is always the root
  
  # Put listcov back in the sptree object, required in predict.sptree
  for (i in 1:length(forest$trees)) {
    forest$trees[[i]]$listcov <- forest$listcov
  }
  
  OOBscr <- parallel::mclapply(1:length(forest$trees), FUN = function(i) {
    OOBval <- rep(NA, X$n)
    
    if (forest$p == 0) {
      torm <- unique(forest$pt_intree[[i]])
      if (length(torm) == X$n) { # If all points are drawn in the bootstrap,
        # then nothing do do.
        return(OOBval)
      }
      OOBpts <- 1:X$n
      OOBpts <- OOBpts[!(OOBpts %in% torm)]
    } else {
      # vector of same length as number of pts in X
      OOBpts <- (forest$pt_intree[[i]] != 1)
      if (all(!OOBpts)) { # If no points in OOBpts, then nothing do do.
        return(OOBval)
      }
    }
    
    # OOB sample
    Xout <- X[OOBpts]
    
    ### OOB prediction
    pts_pred_OOB <- f(
      object = forest$trees[[i]],
      newdata = Xout
    )
    
    OOBval[OOBpts] <- pts_pred_OOB
    
    
    # OOB score
    return(OOBval)
  }, mc.cores = cores)
  
  logterm <- log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE))
  if (all(is.na(logterm))) {
    output <- NA
  } else {
    output <- sum(logterm, na.rm = TRUE)
  }
  # output <- sum(log(rowMeans(do.call(cbind, OOBscr), na.rm = TRUE)), na.rm = TRUE)
  
  # Return the average error of all the trees
  return(output)
}

forest <- RforestPP(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  score = "lcv2",
  p = 1 / 2,
  Ntree = 10,
  threshold = spatstat.geom::area(spatstat.data::bei) / 2^4,
  cores_trees = 1,
  mtry = 1 / 3,
  tol = Inf,
  minpts = 50,
  minsplitq = 0.5,
  maxsplitq = 0.5
)

microbenchmark(OOBscr.spforest(forest=forest, cores=1), times = 100)
microbenchmark(g(forest=forest, cores=1), times = 100)

profvis(g(forest=forest, cores=1))

A<-g(forest=forest, cores=2)
B<-OOBscr.spforest(forest=forest, cores=2)

A-B
