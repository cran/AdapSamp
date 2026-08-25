#' @keywords internal
#' @importFrom pracma fderiv
#' @importFrom stats runif integrate D
#' @importFrom utils tail
NULL

#' Adaptive Rejection Sampling Algorithm
#'
#' rARS generates a sequence of random numbers using the adaptive rejection sampling algorithm.
#'
#' @param n Desired sample size (positive integer).
#' @param formula Kernel of the target density, as a character string (e.g. \code{"exp(-x^2/2)"}).
#' @param min,max Domain of the target distribution, including \code{-Inf} and \code{Inf}.
#' @param sp Supporting set; a numeric vector of initial abscissae.
#' @return A numeric vector of length \code{n} containing samples from the target distribution.
#' @export
#' @author Dong Zhang
#'
#' @examples
#'
#' # Example 1: Standard normal distribution
#' set.seed(123)
#' x1 <- rARS(100, "exp(-x^2/2)", -Inf, Inf, c(-2, 2))
#' hist(x1, breaks = 20, probability = TRUE)
#'
#' # Example 2: Truncated normal distribution
#' \dontrun{
#' x2 <- rARS(100, "exp(-x^2/2)", -2.1, 2.1, c(-2, 2))
#' }
#'
#' # Example 3: Exponential distribution with rate=3
#' \dontrun{
#' x4 <- rARS(100, "exp(-3*x)", 0, Inf, c(2, 3, 100))
#' }
#'
#' # Example 4: Beta distribution with alpha=3 and beta=4
#' \dontrun{
#' x5 <- rARS(100, "x^2*(1-x)^3", 0, 1, c(0.4, 0.6))
#' }
#'
#' # Example 5: Gamma distribution with alpha=5 and lambda=2
#' \dontrun{
#' x6 <- rARS(100, "x^(5-1)*exp(-2*x)", 0, Inf, c(1, 10))
#' }
#'
rARS <- function(n, formula, min = -Inf, max = Inf, sp) {
  ## --- input validation ---
  sp <- sort(sp)
  if (!is.character(formula))
    stop("'formula' must be a character string.")
  if (n <= 0)
    stop("Sample size 'n' must be a positive integer.")
  if (min >= max)
    stop("'min' must be strictly less than 'max'.")
  if (length(sp) < 2)
    stop("Supporting set 'sp' must contain at least 2 points.")

  ## --- density & negative log-density ---
  p <- function(x) eval(parse(text = formula))
  V <- function(x) {
    px <- p(x)
    ifelse(px > 0, -log(px), Inf)
  }

  x_final <- numeric(n)

  for (j in 1:n) {
    Support <- sp
    u <- 0
    compareprop <- -1

    while (u > compareprop) {
      tangent <- fderiv(V, Support, 1)   # numerical derivatives from pracma

      ## cross-points of adjacent tangent lines
      k <- length(Support)
      crosspoint <- numeric(k + 1)
      crosspoint[1] <- min
      crosspoint[k + 1] <- max
      crossvalue <- numeric(k - 1)

      for (i in 1:(k - 1)) {
        A <- matrix(c(tangent[i], -1, tangent[i + 1], -1), nrow = 2, byrow = TRUE)
        b <- c(tangent[i] * Support[i] - V(Support)[i],
               tangent[i + 1] * Support[i + 1] - V(Support)[i + 1])
        sol <- solve(A, b)          # solve once, reuse result
        crosspoint[i + 1] <- sol[1]
        crossvalue[i] <- sol[2]
      }

      ## per-segment integrals of the exponential envelope
      IntSum <- numeric(k)
      for (i in 1:k) {
        t_i <- tangent[i]
        if (abs(t_i) < 1e-12) {
          ## degenerate (flat) segment: uniform envelope
          IntSum[i] <- crosspoint[i + 1] - crosspoint[i]
        } else {
          expfun <- function(x) exp(-t_i * (x - Support[i]) - V(Support)[i])
          IntSum[i] <- integrate(expfun, crosspoint[i], crosspoint[i + 1])$value
        }
      }

      ## inverse-CDF sampling from the envelope
      rdm <- runif(1)
      cum <- c(0, cumsum(IntSum / sum(IntSum)))
      idx <- which(rdm < cumsum(IntSum / sum(IntSum)))[1]

      t_idx <- tangent[idx]
      if (abs(t_idx) < 1e-12) {
        ## uniform segment
        x_star <- runif(1, crosspoint[idx], crosspoint[idx + 1])
      } else {
        x_star <-
          log((rdm - cum[idx] +
                 exp(t_idx * Support[idx] - V(Support)[idx]) *
                 exp(-t_idx * crosspoint[idx]) / sum(IntSum) / (-t_idx)) *
                sum(IntSum) * (-t_idx) /
                exp(t_idx * Support[idx] - V(Support)[idx])) / (-t_idx)
      }

      ## accept / reject
      u <- runif(1)
      compareprop <-
        p(x_star) / exp(-t_idx * (x_star - Support[idx]) - V(Support)[idx])

      ## adapt: add accepted candidate to the support set
      Support <- sort(c(Support, x_star))
    }
    x_final[j] <- x_star
  }
  x_final
}
