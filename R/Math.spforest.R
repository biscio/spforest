#' Dispatch standard math functions
#'
#' @param x An object of class \code{spforestmesh}
#' @param ... Additional arguments passed to the function
#'
#' @returns An object of class \code{spforestmesh}
#' @export
#' 
#' @examples
#' res <- pptomesh(X=spatstat.data::bei,
#' elev= spatstat.data::bei.extra$elev)
#' forest <- spforest(X = res)
#' log(forest+exp(-8))
Math.spforestmesh <- function(x, ...) {
  m <- do.call(.Generic, list(x$tridensity, ...))
  rslt <- list(
    tridensity = m,
    mesh = x$mesh,
    pp = x$pp
  )
  class(rslt) <- "spforestmesh"
  return(rslt)
}



#' Standard binary operations Ops
#'
#' @param e1 An object of class \code{spforestmesh}
#' @param e2 Either a numeric or an object of class \code{spforestmesh}
#'
#' @returns An object of class \code{spforestmesh}
#' @export
#'
#' @examples 
#' res <- pptomesh(X=spatstat.data::bei,
#' elev= spatstat.data::bei.extra$elev)
#' forest <- spforest(X = res)
#' forest + forest
Ops.spforestmesh <- function(e1, e2 = NULL) {
  # Do test compatibility objects

  if (is.null(e2)) {
    #' unary operation
    return(e1)
  }

  if (class(e2) != "spforestmesh") {
    m <- do.call(.Generic, list(e1$tridensity, e2))
  } else {
    m <- do.call(.Generic, list(e1$tridensity, e2$tridensity))
  }

  rslt <- list(
    tridensity = m,
    mesh = e1$mesh,
    pp = e1$pp
  )
  class(rslt) <- "spforestmesh"
  return(rslt)
}
