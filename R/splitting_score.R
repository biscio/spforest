#' Score of a split
#'
#' @description
#' It is used to compute the score.
#'
#'
#' @param n1 Number of points in left child
#' @param n2 Number of points in right child
#' @param W1area Area of left child
#' @param W2area Area of right child
#' @param score The score to use
#'
#' @return A numeric number
#'
#' @details
#' Provide a lot of infor blabla
#'
#' @export
#'
#' @examples
#' score.split(5, 5, 5, 5)
score.split <- function(n1, n2, W1area, W2area, score = "lcv") {
  stopifnot(score %in% c("lcv", "lcv2", "ent", "star", "ise", "isecv"))

  # Return value depending on the score considered
  ## likelihood cross validation
  ## QUESTION: quoi mettre si n1 ou n2 = 0 ?? It is -infinite for now
  if (score == "lcv") {
    val <- ifelse(n1 > 1, n1 * log((n1 - 1) / W1area), -Inf) +
      ifelse(n2 > 1, n2 * log((n2 - 1) / W2area), -Inf)
    return(val)
  }

  # After some simulations it does not appear to change the result
  if (score == "lcv2") {
    val <- ifelse(n1 > 1, n1 * log((n1 - 1) / W1area), 0) +
      ifelse(n2 > 1, n2 * log((n2 - 1) / W2area), 0)
    return(val)
  }

  ## entropie de Poisson (idem lcv sans le -1)
  if (score == "ent") {
    val <- ifelse(n1 >= 1, n1 * log(n1 / W1area), 0) +
      ifelse(n2 >= 1, n2 * log(n2 / W2area), 0)
    return(val)
  }

  ## star discrepancy
  if (score == "star") {
    n <- n1 + n2
    val <- abs(n1 / n - W1area / (W1area + W2area))
    return(val)
  }

  ## integrated square errors (for density estimation)
  ## As I maximise the score in the code, I actually put minus the ISE.
  if (score == "ise") {
    val <- -n1^2 / W1area - n2^2 / W2area
    return(val)
  }

  ## CV integrated square errors (for density estimation)
  if (score == "isecv") {
    val <- -(n1^2 - 2 * n1) / W1area - (n2^2 - 2 * n2) / W2area
    return(val)
  }

  return(val)
}



#' Score parent
#'
#' @param X A ppp object for the parent
#' @param score The score used
#'
#' @return A number
#' @export
#'
#' @examples
#' score.pp(X = spatstat.random::rpoispp(10))
score.pp <- function(X, score = "lcv") {
  stopifnot(score %in% c("lcv", "lcv2", "ent", "star", "ise", "isecv"))
  stopifnot(spatstat.geom::is.ppp(X))

  n <- spatstat.geom::npoints(X)
  W <- spatstat.geom::area(X)

  if (score == "lcv") {
    val <- ifelse(n > 1, n * log((n - 1) / W), -Inf)
  }

  if (score == "lcv2") {
    val <- ifelse(n > 1, n * log((n - 1) / W), 0)
  }

  ## entropie de Poisson (idem lcv sans le -1)
  if (score == "ent") {
    val <- ifelse(n > 1, n * log(n / W), -Inf)
  }

  ## integrated square errors (for density estimation)
  ## As I maximise the score in the code, I actually put minus the ISE.
  if (score == "ise") {
    val <- -n^2 / W
  }

  ## CV integrated square errors (for density estimation)
  if (score == "isecv") {
    val <- -(n^2 - 2 * n) / W
  }
  return(val)
}
