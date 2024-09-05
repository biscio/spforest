finttiles <- function(x,y, tile){
  po <- tile$bdry[[1]]
  # Extracting side segments from the poly
  x0 <- y0 <- x1 <- y1 <- numeric(0)
  ni <- length(po$x)
  nxt <- c(2:ni, 1)
  x0 <- c(x0, po$x)
  y0 <- c(y0, po$y)
  x1 <- c(x1, po$x[nxt])
  y1 <- c(y1, po$y[nxt])
  
  hit <- NULL
  
  for (j in seq_along(x)) {
    if (sum(po$x == x[j]) > 0 | sum(po$y ==y[j]) > 0) {
      warning("The point is a vertex of the polygon")
      hit[j] <- FALSE
    }
    
    cnt <- 0
    for (i in seq_along(x0)) {
      xint <- x0[i] + (y[j] - y0[i]) * (x1[i] - x0[i]) / (y1[i] - y0[i])
      if (((y[j] < y0[i]) != (y[j] < y1[i])) & (x[j] < xint)) {
        cnt <- cnt + 1
      }
    }
    hit[j] <- (cnt %% 2 == 1)
  }
  
  return(hit)
}
mytileindex <-  function (x, y, Z) {
  stopifnot(is.tess(Z))
  if ((missing(y) || is.null(y)) && all(c("x", "y") %in% names(x))) {
    y <- x$y
    x <- x$x
  }
  stopifnot(length(x) == length(y))
  switch(Z$type, rect = {
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
                       sep = "")
  }, tiled = {
    n <- length(x)
    todo <- seq_len(n)
    nt <- length(Z$tiles)
    m <- integer(n)
    for (i in 1:nt) {
      ti <- Z$tiles[[i]]
      hit <- finttiles(x[todo], y[todo], ti)
      if (any(hit)) {
        m[todo[hit]] <- i
        todo <- todo[!hit]
      }
      if (length(todo) == 0) break
    }
    m[m == 0] <- NA
    nama <- names(Z$tiles)
    lev <- seq_len(nt)
    lab <- if (!is.null(nama) && all(nzchar(nama))) nama else paste("Tile", 
                                                                    lev)
    m <- factor(m, levels = lev, labels = lab)
  }, image = {
    Zim <- Z$image
    m <- lookup.im(Zim, x, y, naok = TRUE)
    if (anyNA(m)) {
      isna <- is.na(m)
      rc <- nearest.valid.pixel(x[isna], y[isna], Zim, 
                                nsearch = 2)
      m[isna] <- Zim$v[cbind(rc$row, rc$col)]
    }
  })
  return(m)
}


tesstree3 <- function(X,
                      lambda = 100,
                      dimyx = c(128, 128),
                      test.connected = FALSE) {
  # lambda is the nb of points (not the intensity)
  wind <- Window(X)
  enclose.rect <- spatstat.geom::owin(wind$xrange, wind$yrange)
  
  tol <- 2 / sqrt(lambda) # in order to add points outside the window
  # Alternative to avoid a polygonal windows.
  tess.points <- spatstat.geom::runifrect(
    lambda,
    owin(
      xrange = enclose.rect$xrange + c(-tol, tol),
      yrange = enclose.rect$yrange + c(-tol, tol)
    )
  )
  
  del <- spatstat.geom::dirichlet(tess.points) # associated tessellation
  tmp <- spatstat.geom::intersect.tess(del, wind) # intersected with the windows
  if (test.connected) {
    tmp <- spatstat.geom::connected(tmp)
  }
  
  delarea <- spatstat.geom::tile.areas(tmp) # collect the areas
  
  # mX <- marks(cut(X, tmp)) ## the alternative is very slightly quicker
  mX <- g(X$x, X$y, tmp)
  
  ptintess <- c()
  if (test.connected) {
    for (i in levels(tmp$image)) {
      ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
    }
  } else {
    for (i in names(tmp$tiles)) {
      ptintess <- c(ptintess, sum(mX == i, na.rm = TRUE))
    }
  }
  return(as.im(tmp, values = ptintess / delarea, dimyx = dimyx))
}
