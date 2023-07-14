#' Intensity tesselation tree
#'
#' @param X
#' @param lambda
#' @param target.points
#' @param test.connected
#'
#' @return
#' @export
#'
#' @examples
tesstree <- function(X,
                     lambda,
                     target.points,
                     test.connected = TRUE) {
  # lambda is the nb of points (not the intensity)
  wind <- Window(X)
  enclose.rect <- spatstat.geom::owin(wind$xrange, wind$yrange)

  tol <- 2 / sqrt(lambda) # in order to add points outside the window
  # Alternative to avoid a polygonal windows. 
  tess.points <- spatstat.geom::runifrect(
    lambda,
    owin(xrange = enclose.rect$xrange + c(-tol, tol), 
         yrange = enclose.rect$yrange + c(-tol, tol))
  )
  # tess.points <- spatstat.geom::runifrect(
  #   lambda,
  #   spatstat.geom::dilation(
  #     enclose.rect,
  #     tol
  #   )
  # ) 
  # simulate the dummy points

  del <- spatstat.geom::dirichlet(tess.points) # associated tesselation
  tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows

  if (test.connected) {
    tmp <- spatstat.geom::connected(tmp)
  } # split unconnected tiles
  # tmp<-connected(tmp,dimyx=c(512,512)) #supply dimyx to improve precision (better way to deal with connected?)
  delarea <- spatstat.geom::tile.areas(tmp) # collect the areas


  # nn of target.points in the dummy points
  # corresponds to the index of tile in tmp if no split has been done by connected()
  ind <- spatstat.geom::nncross(target.points,
    tess.points,
    what = "which"
  )

  # check if some cells were splitted by connected()
  if (test.connected) {
    indpb <- NULL
    for (i in 1:target.points$n) {
      if (!as.character(ind[i]) %in% names(delarea)) {
        indpb <- c(indpb, i)
      }
    }

    # repair ind for the splitted tiles due to connected()
    if (!is.null(indpb)) {
      fun2 <- function(v, a) names(a)[as.logical(v)]
      tile.tmp <- spatstat.geom::tiles(tmp)
      a <- lapply(tile.tmp, FUN = function(jj) {
        spatstat.geom::inside.owin(target.points[indpb], w = jj)
      })
      aa <- do.call(cbind, a)
      # fun <- function(w, X) spatstat.geom::inside.owin(X, w = w)
      # a <- lapply(tile.tmp, fun, target.points[indpb])
      # aa <- as.data.frame(a)
      
      newnames <- apply(aa, 1, fun2, a)
      ind[indpb] <- newnames
    }
  }
  ind <- as.character(ind)
  ind.vec <- spatstat.geom::nncross(X,
    tess.points,
    what = "which"
  )

  if (test.connected) {
    indpb <- NULL
    for (i in 1:X$n) {
      if (!as.character(ind.vec[i]) %in% names(delarea)) {
        indpb <- c(indpb, i)
      }
    }
    if (!is.null(indpb)) {
      tile.tmp <- tiles(tmp)
      Y <- X[indpb]
      spts <- sp::SpatialPoints(list(Y$x, Y$y))
      nn <- length(tile.tmp)
      dist.tile <- NULL
      for (i in 1:nn) {
        b <- spatstat.geom::as.polygonal(tile.tmp[[i]])
        bb <- list()
        for (k in 1:length(b$bdry)) {
          if (k == 1) {
            bb[[k]] <- sp::Polygon(matrix(unlist(b$bdry[[k]]), ncol = 2),
              hole = FALSE
            )
          } else {
            bb[[k]] <- sp::Polygon(matrix(unlist(b$bdry[[k]]), ncol = 2),
              hole = TRUE
            )
          }
        }
        firstPoly <- sp::Polygons(bb, ID = names(delarea)[i])
        firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
        # rgeos::gDistance will be deprecated in 2023.
        # Alternative sf::st_distance() 	terra::distance()
        # dist.tile <- rbind(dist.tile, apply(gDistance(spts,
        #                                               firstSpatialPoly,
        #                                               byid = TRUE), 2, min))
        dist.tile <- rbind(dist.tile, apply(sf::st_distance(
          st_as_sf(spts),
          st_as_sf(firstSpatialPoly),
          by_element = F
        ), 2, min))
      }
      index <- apply(dist.tile, 2, which.min)
      ind.vec[indpb] <- names(delarea)[index]
    }
  }
  ind.vec <- as.factor(ind.vec)

  noms <- names(delarea)
  card <- as.table(rep(0, length(noms)))
  names(card) <- noms
  card0 <- table(ind.vec) # number of points of X in each tile
  card[names(card0)] <- card0

  return(card[ind] / delarea[ind])
}


#' Intensity tesselation forest
#'
#' @param X
#' @param lambda
#' @param Ntrees
#' @param at
#' @param mc.cores
#' @param test.connected
#'
#' @return
#' @export
#'
#' @examples
tessforest <- function(X,
                       lambda = NULL,
                       Ntrees = 100,
                       at = NULL,
                       mc.cores = 1,
                       test.connected = TRUE) {
  # at : a ppp where the intensity is estimated (by default on a 128x128 image)
  if (is.null(at)) {
    wind <- spatstat.geom::Window(X)
    N <- 128 # size of image
    x.image <- seq(wind$xrange[1], wind$xrange[2], length.out = N)
    y.image <- seq(wind$yrange[1], wind$yrange[2], length.out = N)
    allpoints <- spatstat.geom::as.ppp(
      expand.grid(x.image, y.image),
      owin(wind$xrange, wind$yrange)
    )
    target.points <- spatstat.geom::subset.ppp(allpoints, wind)
  } else {
    target.points <- at
    N <- at$n
  }

  if (is.null(lambda)) {
    lambda <- floor(mean(c(
      grDevices::nclass.FD(X$x),
      grDevices::nclass.FD(X$y)
    ))^2)
  }

  fun <- function(i) {
    tesstree(
      X,
      lambda,
      target.points,
      test.connected
    )
  }
  tmp <- unlist(mclapply(1:Ntrees, fun, mc.cores = mc.cores))
  res <- rowMeans(matrix(tmp, nrow = length(tmp) / Ntrees)) # TODO: use columns ?

  if (is.null(at)) {
    marks(allpoints) <- NA
    a <- spatstat.geom::inside.owin(allpoints, w = wind)
    marks(allpoints)[a] <- res
    return(as.im(t(matrix(marks(allpoints), N, N)),
      W = owin(wind$xrange, wind$yrange)
    ))
  } else {
    return(res)
  }
}
