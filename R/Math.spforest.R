#' Dispatch standard math functions
#'
#' @param x
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
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
#' @param e1
#' @param e2
#'
#' @returns
#' @export
#'
#' @examples
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
