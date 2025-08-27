
#' Change resolution of image
#'
#' @param Z An im object from spatstat
#' @param dimyx The new values of dimyx
#' @param xstep If dimyx not specified, the x-size of the pixel 
#' @param ystep  If dimyx not specified, the y-size of the pixel 
#'
#' @returns An im object
#' @export
#'
#' @examples
#' beismooth <- lapply(beisoilres, imsmooth, dimyx = c(256, 512))
#' plot(as.anylist(beismooth))
imsmooth <- function(Z, dimyx = NULL, xstep, ystep = xstep) {
  if (is.null(dimyx)) {
    x <- seq(Z$xrange[1] + xstep / 2,
             Z$xrange[2] - xstep / 2,
             by = xstep
    )
    y <- seq(Z$yrange[1] + ystep / 2,
             Z$yrange[2] - ystep / 2,
             by = ystep
    )
  } else {
    x <- seq(Z$xrange[1] + (diff(Z$xrange) / dim(Z)[2]) / 2,
             Z$xrange[2] - (diff(Z$xrange) / dim(Z)[2]) / 2,
             length.out = dimyx[2]
    )
    y <- seq(Z$yrange[1] + (diff(Z$yrange) / dim(Z)[1]) / 2,
             Z$yrange[2] - (diff(Z$yrange) / dim(Z)[1]) / 2,
             length.out = dimyx[1]
    )
  }
  mygrid <- expand.grid(x=x, y=y) 
  mygrid <- mygrid[!duplicated.data.frame(mygrid),]
  
  Bv <- interp.im(Z, x = mygrid$x, y = mygrid$y, bilinear = T)
  Bmat <- matrix(Bv, nrow = length(x), ncol = length(y))
  
  return(im(t(Bmat), xrange = Z$xrange, yrange = Z$yrange))
}
