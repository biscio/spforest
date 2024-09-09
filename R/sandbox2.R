# Z <- spforest(
#   X = bei,
#   lambda = 100,
#   dimyx = c(101, 201),
#   test.connected = FALSE,
#   Ntree = 5,
#   cores = 1
# )
#
# Z2 <- spforest(
#   X = spatstat.data::bei,
#   listcovariates = spatstat.data::bei.extra,
#   Ntree = 3,
#   minpts = 100,
#   mtry = 1 / 3,
#   p = 0,
#   cores = 1
# )
#
# plot(Z2)

# library(microbenchmark)
#
# fnew <- function() {
#   tessforest3(
#     X = bei,
#     lambda = 100,
#     dimyx = c(128, 128),
#     test.connected = FALSE,
#     Ntree = 5,
#     cores = 1
#   )
# }
#
# fold <- function() {
#   tessforest(
#     X = bei,
#     lambda = 100,
#     Ntree = 5,
#     test.connected = FALSE
#   )
# }
#
# microbenchmark(fnew(), fold(), times = 2)
#
# profvis::profvis(tessforest2(
#   X = bei,
#   lambda = 100,
#   dimyx = c(128, 128),
#   test.connected = FALSE,
#   Ntree = 5,
#   cores = 1
# ))
#
#
# myinsideowin <- function(x, y, w) {
#   if ((missing(y) || is.null(y)) && all(c("x", "y") %in% names(x))) {
#     y <- x$y
#     x <- x$x
#   }
#   # w <- as.owin(w)
#   if (length(x) == 0) {
#     return(logical(0))
#   }
#   xr <- w$xrange
#   yr <- w$yrange
#   eps <- sqrt(.Machine$double.eps)
#   frameok <- (x >= xr[1L] - eps) & (x <= xr[2L] + eps) & (y >=
#     yr[1L] - eps) & (y <= yr[2L] + eps)
#   if (!any(frameok)) {
#     return(frameok)
#   }
#   ok <- frameok
#   switch(w$type,
#     rectangle = {
#       return(ok)
#     },
#     polygonal = {
#       framesize <- max(diff(xr), diff(yr))
#       if (framesize > 1e+06 || framesize < 1e-06) {
#         scalefac <- as.numeric(framesize) / 100
#         w <- as.polygonal(rescale(w, scalefac))
#         x <- x / scalefac
#         y <- y / scalefac
#       }
#       xy <- list(x = x, y = y)
#       bdry <- w$bdry
#       total <- numeric(length(x))
#       on.bdry <- rep.int(FALSE, length(x))
#       for (i in seq_along(bdry)) {
#         score <- inside.xypolygon(xy, bdry[[i]], test01 = FALSE)
#         total <- total + score
#         on.bdry <- on.bdry | attr(score, "on.boundary")
#       }
#       total[on.bdry] <- 1
#       uhoh <- (total * (1 - total) != 0)
#       if (any(uhoh)) {
#         nuh <- sum(uhoh)
#         warning(
#           paste(
#             "point-in-polygon test had difficulty with",
#             nuh, ngettext(nuh, "point", "points"), "(total score not 0 or 1)"
#           ),
#           call. = FALSE
#         )
#         total[uhoh] <- 0
#       }
#       return(ok & (total != 0))
#     },
#     mask = {
#       xf <- x[frameok]
#       yf <- y[frameok]
#       loc <- nearest.raster.point(xf, yf, w)
#       okf <- (w$m)[cbind(loc$row, loc$col)]
#       ok[frameok] <- okf
#       return(ok)
#     },
#     stop("unrecognised window type", w$type)
#   )
# }
#
# mytileid <- function(x, y, Z) {
#   stopifnot(is.tess(Z))
#   if ((missing(y) || is.null(y)) && all(c("x", "y") %in% names(x))) {
#     y <- x$y
#     x <- x$x
#   }
#   stopifnot(length(x) == length(y))
#   switch(Z$type,
#     rect = {
#       jx <- findInterval(x, Z$xgrid, rightmost.closed = TRUE)
#       iy <- findInterval(y, Z$ygrid, rightmost.closed = TRUE)
#       nrows <- length(Z$ygrid) - 1
#       ncols <- length(Z$xgrid) - 1
#       iy[iy < 1 | iy > nrows] <- NA
#       jx[jx < 1 | jx > ncols] <- NA
#       jcol <- jx
#       irow <- nrows - iy + 1
#       ktile <- jcol + ncols * (irow - 1)
#       m <- factor(ktile, levels = seq_len(nrows * ncols))
#       ij <- expand.grid(j = seq_len(ncols), i = seq_len(nrows))
#       levels(m) <- paste("Tile row ", ij$i, ", col ", ij$j,
#         sep = ""
#       )
#     },
#     tiled = {
#       n <- length(x)
#       todo <- seq_len(n)
#       nt <- length(Z$tiles)
#       m <- integer(n)
#       for (i in 1:nt) {
#         ti <- Z$tiles[[i]]
#         hit <- myinsideowin(x[todo], y[todo], ti)
#         if (any(hit)) {
#           m[todo[hit]] <- i
#           todo <- todo[!hit]
#         }
#         if (length(todo) == 0) break
#       }
#       m[m == 0] <- NA
#       nama <- names(Z$tiles)
#       lev <- seq_len(nt)
#       lab <- if (!is.null(nama) && all(nzchar(nama))) {
#         nama
#       } else {
#         paste(
#           "Tile",
#           lev
#         )
#       }
#       m <- factor(m, levels = lev, labels = lab)
#     },
#     image = {
#       Zim <- Z$image
#       m <- lookup.im(Zim, x, y, naok = TRUE)
#       if (anyNA(m)) {
#         isna <- is.na(m)
#         rc <- nearest.valid.pixel(x[isna], y[isna], Zim,
#           nsearch = 2
#         )
#         m[isna] <- Zim$v[cbind(rc$row, rc$col)]
#       }
#     }
#   )
#   return(m)
# }
#
#
#
#
#
#
# # # https://cp-algorithms.com/geometry/point-in-convex-polygon.html
# #
# # ## TODO: handle edges and vertices cases.
# # ## TODO: vectorise
# #
# # dd <- dirichlet(cells)
# #
# # polyg <- dd$tiles[[4]]$bdry[[1]]
# # polyg$x
# # polyg$y
# # plot(polyg$x, polyg$y, col = c(1, 2, 3, 4), pch = 16)
# #
# # orient2d <- function(a1, a2, a3) {
# #   v1 <- a1 - a3
# #   v2 <- a2 - a3
# #   return(v1[1] * v2[2] - v1[2] * v2[1])
# # }
# #
# # reorder_po <- function(po) {
# #   id <- which(po$x == min(po$x, na.rm = TRUE))
# #   if (length(id)>1) {
# #     subid <- which(po$y[id] == min(po$y[id], na.rm = TRUE))
# #     id <- id[subid]
# #   }
# #
# #   px <- c(po$x[id:length(po$x)], c(po$x[1:(id - 1)]))
# #   py <- c(po$y[id:length(po$y)], c(po$y[1:(id - 1)]))
# #   return(list(x = px, y = py))
# # }
# #
# # orient2d(a1 = p1, a2 = pt, a3 = p0)
# #
# # po2 <- reorder_po(polyg)
# #
# # p0 <- c(po2$x[1], po2$y[1])
# # p1 <- c(po2$x[2], po2$y[2])
# # p2 <- c(po2$x[3], po2$y[3])
# # pt <- c(0.75, 0.04)
# # pt2 <- c(0.695, 0.04)
# #
# # plot(polyg)
# # points(p0[1], p0[2], pch = 16, col = 1)
# # points(p1[1], p1[2], pch = 2, col = 2, lwd = 2)
# # points(p2[1], p2[2], pch = 3, col = 4, lwd = 2)
# # points(pt[1], pt[2], pch = 4, col = 6, lwd = 2)
# # points(pt2[1], pt2[2], pch = 4, col = 7, lwd = 2)
# #
# # orient2d(a1 = p1, a2 = p2, a3 = p0)
# # orient2d(a1 = p1, a2 = pt, a3 = p0)
# #
# # orient2d(a1 = p2, a2 = pt2, a3 = p0)
# # areatri(p0, p1, p2) -
# #   areatri(p0, p1, pt) -
# #   areatri(p1, p2, pt) -
# #   areatri(p0, p2, pt)
# #
# # abs(orient2d(p0, p1, p2)) -
# #   abs(orient2d(p0, p1, pt)) -
# #   abs(orient2d(p1, p2, pt)) -
# #   abs(orient2d(p0, p2, pt))
# #
# # intri <- function(pt, a1, a2, a3, eps=1e-10){
# #   x <- abs(orient2d(a1, a2, a3)) -
# #     abs(orient2d(a1, a2, pt)) -
# #     abs(orient2d(a2, a3, pt)) -
# #     abs(orient2d(a1, a3, pt))
# #   return(abs(x) < eps)
# # }
# #
# # inpoly <- function(pt, polyg, eps = 1e-10) {
# #   po <- reorder_po(polyg)
# #   n <- length(po$x)
# #   p1 <- c(po$x[1], po$y[1])
# #   p2 <- c(po$x[2], po$y[2])
# #   pn <- c(po$x[n], po$y[n])
# #
# #   if (orient2d(a1 = p2, a2 = pt , a3 = p1)<0) {
# #     return(FALSE)
# #   }
# #
# #   if (orient2d(a1 = pt, a2 = pn , a3 = p1)<0) {
# #     return(FALSE)
# #   }
# #
# #   z <- NULL
# #   for (i in 2:(n-1)) {
# #     b1 <-  c(po$x[i], po$y[i])
# #     b2 <-  c(po$x[i+1], po$y[i+1])
# #     z[i-1] = intri(pt=pt, a1=p1, a2=b1, a3=b2, eps=1e-10)
# #   }
# #
# #   return(sum(z)>0)
# # }
# #
# # inpoly(pt=pt, polyg=polyg, eps = 1e-10)
# #
# # polyg2 <- dd$tiles[[18]]$bdry[[1]]
# # plot(polyg2)
# # pt <- c(0.45,0.45)
# # points(pt[1], pt[2], pch=16)
# #
# # po2 <- reorder_po(polyg2)
# # for (i in 1:length(polyg2$x)) {
# #   points(po2$x[i], po2$y[i], col=i, pch=i, lwd=i)
# # }
# #
# # inpoly(pt=pt, polyg=polyg2, eps = 1e-10)
# #
