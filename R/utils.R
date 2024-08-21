
#' Find the small pixel area in list of images
#'
#' @param x an im object
#'
#' @return A number
#' @export 
#'
#' @examples
smallest_pixelarea <- function(x) {
  allarea <- sapply(x, FUN = function(i) {
    y <- unclass(i)[c("xstep", "ystep")]
    pixelarea <- y$xstep * y$ystep
  })

  return(min(allarea))
}
