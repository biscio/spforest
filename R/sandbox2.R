insideconvex <- function(p, polyg) {
  orient <- NULL
  for (i in 1:(dim(polyg)[1] - 1)) {
    v1 <- polyg[i + 1, ] - polyg[i, ]
    v2 <- p - polyg[i, ]
    orient[i] <- v1[1] * v2[2] - v1[2] * v2[1]
  }
  return(all(orient > 0) | all(orient<0))
}

X <- runifpoint(5)
dd <- dirichlet(cells)
plot(dd); points(X)
tileindex(X$x, X$y, dd)

mytileid <- function(X, dd){
  m1 <- NULL
  m2 <- NULL
  for (i in 1:length(dd$tiles)) {
      polyg <- cbind(dd$tiles[[i]]$bdry[[1]]$x, dd$tiles[[i]]$bdry[[1]]$y)
      A <- NULL
      B <- NULL
    for (j in 1:npoints(X)) {
      A[j] <- insideconvex(c(X$x[j], X$y[j]), polyg = polyg)
      B[j] <- is_point_inside(point = c(X$x[j], X$y[j]), polyg = polyg)
    }
    m1 <- cbind(m1, A)
    m2 <- cbind(m2, B)
  }
  
  output <- NULL
  for (i in 1:npoints(X)) {
    output[i] <- which(m2[i, ])
  }
  rowsum(m1) ### ISSUES m1, the insideconvex is not working
  
}



library(profvis)
X <- runifpoint(3)
dd <- dirichlet(cells)

test<- cut(X, B) |> marks()
which(test==22)
polyg <- cbind(B$tiles[[22]]$bdry[[1]]$x, B$tiles[[22]]$bdry[[1]]$y)

is_point_inside(point = c(X$x[30], X$y[30]), polyg = polyg)
is_point_inside(point = c(X$x[31], X$y[31]), polyg = polyg)

A<-NULL
B<-NULL
for (i in 1:X$n) {
  A[i] <- is_point_inside(point = c(X$x[i], X$y[i]), polyg = polyg)
  B[i] <-  f(p = c(X$x[i], X$y[i]), polyg = polyg)
}

all(A == B)
B


is_point_inside <-  function(point, polyg) {
    p <- as.numeric(point)
    # library(mgcv)
    # return(in.out(as.matrix(polyg), p))
    is.vertex <- sum(apply(polyg, 1, function(x) all(x == p)))
    if (is.vertex == 1) {
      return(FALSE)
    }
    sx1 <- sum(!p[1] <= polyg[, 1])
    sx2 <- sum(!p[1] >= polyg[, 1])
    if (sx1 == 0 | sx2 == 0) {
      return(FALSE)
    }
    sy1 <- sum(!p[2] <= polyg[, 2])
    sy2 <- sum(!p[2] >= polyg[, 2])
    if (sy1 == 0 | sy2 == 0) {
      return(FALSE)
    }
    px <- polyg[, 1]
    py <- polyg[, 2]
    px <- c(px, px[1])
    py <- c(py, py[1])
    segments <- (p[2] - py) / (p[1] - px)
    condition <- segments[-1] == segments[-length(segments)]
    is.segment <- sum(condition)
    if (is.segment > 0) {
      c1 <- px[-1][condition]
      c2 <- px[-length(px)][condition]
      if (c1 < c2) {
        if (c1 < p[1] & p[1] < c2) {
          return(FALSE)
        }
      } else {
        if (c2 < p[1] & p[1] < c1) {
          return(FALSE)
        }
      }
    }
    xcross <- polyg[, 1] + (px[-1] - px[-length(px)]) * (p[2] - polyg[, 2]) / (py[-1] - py[-length(py)])
    xcross <- c(xcross[length(xcross)], xcross[-length(xcross)])
    px2 <- polyg[, 1]
    px2 <- c(px2[length(px2)], px2[-length(px2)])
    px1 <- px[-length(px)]
    xcross2 <- xcross[p[1] <= xcross]
    px1 <- px1[p[1] <= xcross]
    px2 <- px2[p[1] <= xcross]
    n.intersections1 <- sum(px1 <= xcross2 & xcross2 <= px2)
    n.intersections2 <- sum(px2 <= xcross2 & xcross2 <= px1)
    n.int <- n.intersections1 + n.intersections2
    if (any(p[2] == polyg[, 2][p[1] < polyg[, 1]])) n.int <- n.int - 1
    if (n.int == 1) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

plot(B)
points(X)
points(p[1], p[2], pch = 16, cex = 0.6, col = 5)
points(polyg[1, 1], polyg[1, 2])
points(polyg[2, 1], polyg[2, 2], col = 2)


####



f <- function() {
  is_point_inside(point = c(X$x[1], X$y[1]), polyg = polyg)
}
g <- function() {
  inside.owin(X$x[1], X$y[1], B$tiles[[1]])
}
microbenchmark(f(), g())
is_point_inside(point = c(X$x[1], X$y[1]), polyg = polyg)

test <- function() {
  function(x, y, Z) {
    stopifnot(is.tess(Z))
    if ((missing(y) || is.null(y)) && all(c("x", "y") %in% names(x))) {
      y <- x$y
      x <- x$x
    }
    stopifnot(length(x) == length(y))
    switch(Z$type,
      rect = {
        jx <- findInterval(x, Z$xgrid, rightmost.closed = TRUE)
        iy <- findInterval(y, Z$ygrid, rightmost.closed = TRUE)
        nrows <- length(Z$ygrid) - 1
        ncols <- length(Z$xgrid) - 1
        iy[iy < 1 | iy > nrows] <- NA
        jx[jx < 1 | jx > ncols] <- NA
        jcol <- jx
        irow <- nrows - iy + 1
        ktile <- jcol + ncols * (irow - 1)
        m <- factor(ktile, levels = seq_len(nrows * ncols))
        ij <- expand.grid(j = seq_len(ncols), i = seq_len(nrows))
        levels(m) <- paste("Tile row ", ij$i, ", col ", ij$j,
          sep = ""
        )
      },
      tiled = {
        n <- length(x)
        todo <- seq_len(n)
        nt <- length(Z$tiles)
        m <- integer(n)
        for (i in 1:nt) {
          ti <- Z$tiles[[i]]
          hit <- inside.owin(x[todo], y[todo], ti)
          if (any(hit)) {
            m[todo[hit]] <- i
            todo <- todo[!hit]
          }
          if (length(todo) == 0) break
        }
        m[m == 0] <- NA
        nama <- names(Z$tiles)
        lev <- seq_len(nt)
        lab <- if (!is.null(nama) && all(nzchar(nama))) {
          nama
        } else {
          paste(
            "Tile",
            lev
          )
        }
        m <- factor(m, levels = lev, labels = lab)
      },
      image = {
        Zim <- Z$image
        m <- lookup.im(Zim, x, y, naok = TRUE)
        if (anyNA(m)) {
          isna <- is.na(m)
          rc <- nearest.valid.pixel(x[isna], y[isna], Zim,
            nsearch = 2
          )
          m[isna] <- Zim$v[cbind(rc$row, rc$col)]
        }
      }
    )
    return(m)
  }
}

Z <- B
x <- X$x
y <- X$y



test2 <- function(x, y, Z) {
  stopifnot(is.tess(Z))
  if ((missing(y) || is.null(y)) && all(c("x", "y") %in% names(x))) {
    y <- x$y
    x <- x$x
  }
  stopifnot(length(x) == length(y))

  n <- length(x)
  todo <- seq_len(n)
  nt <- length(Z$tiles)
  m <- integer(n)
  for (i in 1:nt) {
    hit <- inside.owin(x[todo], y[todo], Z$tiles[[i]])
    if (any(hit)) {
      m[todo[hit]] <- i
      todo <- todo[!hit]
    }
    if (length(todo) == 0) break
  }
  m[m == 0] <- NA
  nama <- names(Z$tiles)
  lev <- seq_len(nt)
  lab <- if (!is.null(nama) && all(nzchar(nama))) {
    nama
  } else {
    paste(
      "Tile",
      lev
    )
  }
  m <- factor(m, levels = lev, labels = lab)


  return(m)
}
