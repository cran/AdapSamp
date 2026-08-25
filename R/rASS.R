#' Adaptive Slice Sampling Algorithm With Stepping-Out Procedures
#'
#' rASS generates a sequence of random numbers by the adaptive slice sampling
#' algorithm with stepping-out procedures (Neal, 2003).
#'
#' @param n Desired sample size (positive integer).
#' @param x0 Initial value for the Markov chain.
#' @param formula Target density function \code{p(x)}, as a character string.
#' @param w Length of the initial coverage interval.
#' @param max_steps Maximum number of stepping-out iterations (default 1000).
#' @return A numeric vector of length \code{n} containing samples (excludes the initial value \code{x0}).
#' @references Neal R M. Slice sampling - Rejoinder[J]. Annals of Statistics, 2003, 31(3):758-767.
#' @author Dong Zhang
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rASS(100, -1, "1.114283*exp(-(4-x^2)^2)", 3)
#' plot(density(x))
#'
rASS <- function(n, x0 = 0, formula, w = 3, max_steps = 1000) {
  ## --- input validation ---
  if (!is.character(formula))
    stop("'formula' must be a character string.")
  if (n <= 0)
    stop("Sample size 'n' must be a positive integer.")
  if (w <= 0)
    stop("Coverage interval 'w' must be positive.")
  if (max_steps <= 0)
    stop("'max_steps' must be a positive integer.")

  f <- function(x) eval(parse(text = formula))

  ## pre-allocate result vector (length n, excluding x0)
  x_final <- numeric(n)
  x_final[1] <- x0

  for (i in 1:n) {
    ## vertical slice: sample height uniformly in (0, f(x_current)]
    fx_cur <- f(x_final[i])
    if (!is.finite(fx_cur) || fx_cur <= 0) {
      ## fallback: stay at current point
      if (i < n) x_final[i + 1] <- x_final[i]
      next
    }
    slice_h <- runif(1, 0, fx_cur)

    ## stepping-out: expand (left, right) until both ends are below slice
    left <- x_final[i] - runif(1, 0, w)
    right <- left + w
    step_count <- 0
    while (!((f(left) < slice_h) && (f(right) < slice_h))) {
      left <- left - w
      right <- right + w
      step_count <- step_count + 1
      if (step_count > max_steps) break
    }

    ## shrinkage: sample uniformly in (left, right), shrink until acceptance
    x <- runif(1, left, right)
    shrink_count <- 0
    while (f(x) < slice_h) {
      if (x > x_final[i]) {
        right <- x
      } else {
        left <- x
      }
      x <- runif(1, left, right)
      shrink_count <- shrink_count + 1
      if (shrink_count > max_steps) break
    }

    ## accept (use >= to handle floating-point edge case)
    if (f(x) >= slice_h) {
      x_final[i + 1] <- x
    } else {
      ## fallback after max_steps: keep current point
      x_final[i + 1] <- x_final[i]
    }
  }

  ## return only the n samples (drop x0)
  x_final[-1]
}
