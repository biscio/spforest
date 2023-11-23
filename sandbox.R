rm(list = ls())
library(rworldmap)
# library(rworldxtra)
library(raster)
library(spatstat)
# library(maptools) # Soon deprecated
library(sf)
library(parallel)
library(Rsandbox)
# spatstat.options(npixel = 512)

# Dans draft de fred
# {Trout_UK2}
# {UK_Trout_ppl}
# {UK_Trout_tree}
# {UK_Trout_RF}
# {UK_Trout_ppl_clip}
# {UK_Trout_RF_clip}

# Get the map of UK ----

## high create an error on my laptop
res_world_map <- "low" # Other solution high
if (!(res_world_map %in% c("low", "high"))) {
  stop("res_world_map must be set to low or high as a string.")
}


world <- rworldmap::getMap(resolution = res_world_map)
## clip in long lag coordinate a region
clipper_europe <- as(raster::extent(-10, 32, 30, 72), "SpatialPolygons")
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
spatstat.options(npixel = 256)
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
plot(res.tree[[2]], add = T, border = "black", lwd=4)



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


#############3
#############
#############
#############


test<- lapply(1:nn, FUN=function(i){
  spatstat.geom::as.polygonal(tile.tmp[[i]])$bdry
})

concomp <- sapply(test, length)



plot(tmp)
for (i in which(concomp>1)) {
  plot(as.polygonal(tile.tmp[[i]]), add=T, col=2)
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

