#' Concave-Convex Adaptive Rejection Sampling Algorithm
#'
#' rCCARS generates a sequence of random numbers by the concave-convex adaptive
#' rejection sampling algorithm from target distributions with bounded domain.
#'
#' @param n Desired sample size (positive integer).
#' @param cvformula Convex part of \code{-log(p(x))}, as a character string.
#' @param ccformula Concave part of \code{-log(p(x))}, as a character string.
#' @param min,max Bounded domain limits (must be finite).
#' @param sp Supporting set; a numeric vector of initial abscissae.
#' @return A numeric vector of length \code{n} containing samples from the target distribution.
#' @author Dong Zhang
#' @references Teh Y W. Concave-Convex Adaptive Rejection Sampling[J]. Journal of Computational & Graphical Statistics, 2011, 20(3):670-691.
#' @export
#'
#' @details Strictly speaking, the concave-convex adaptive rejection sampling
#'   algorithm can generate samples from target distributions with bounded
#'   domains. For distributions with unbounded domain, rCCARS can also be used
#'   for approximate sampling by truncating the domain.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rCCARS(100, "x+x^-1", "2*log(x)", 0.001, 100, 1)
#' hist(x, breaks = 20, probability = TRUE)
#' }
#'
rCCARS <- function(n, cvformula, ccformula, min, max, sp) {
  ## --- input validation ---
  if (!is.character(cvformula) || !is.character(ccformula))
    stop("'cvformula' and 'ccformula' must be character strings.")
  if (n <= 0)
    stop("Sample size 'n' must be a positive integer.")
  if (min >= max)
    stop("'min' must be strictly less than 'max'.")
  if (!is.finite(min) || !is.finite(max))
    stop("rCCARS requires bounded domain: 'min' and 'max' must be finite.")
  if (length(sp) < 1)
    stop("Supporting set 'sp' must contain at least 1 point.")

  p <- function(x) {
    eval(parse(text = paste0("exp(-(", cvformula, ")-(", ccformula, "))")))
  }

  x_final <- numeric(n)

  for (k in 1:n) {
    support <- sp
    xrange <- c(min, max)
    u <- 0
    prop <- -1

    while (u > prop) {
      allpt <- sort(c(xrange, support))

      ## --- convex part: tangent-line upper envelope ---
      convex <- function(x) eval(parse(text = cvformula))
      drv1or <- function(x) eval(D(parse(text = cvformula), "x"))
      der <- drv1or(allpt)

      crossx <- numeric(length(allpt) - 1)
      crossy <- numeric(length(allpt) - 1)
      for (i in 1:(length(allpt) - 1)) {
        A <- matrix(c(der[i], -1, der[i + 1], -1), nrow = 2, byrow = TRUE)
        b <- c(der[i] * allpt[i] - convex(allpt)[i],
               der[i + 1] * allpt[i + 1] - convex(allpt)[i + 1])
        sol <- solve(A, b)
        crossx[i] <- sol[1]
        crossy[i] <- sol[2]
      }

      ## merge cross-points with support points, sorted
      merged <- data.frame(X = c(crossx, allpt), Y = c(crossy, convex(allpt)))
      merged <- merged[order(merged$X), ]
      ## remove duplicate X (zero-length intervals)
      dup <- duplicated(merged$X)
      if (any(dup)) merged <- merged[!dup, ]
      xconvex <- merged$X
      yconvex <- merged$Y

      n_seg_convex <- length(xconvex) - 1
      tan1 <- numeric(n_seg_convex)
      int1 <- numeric(n_seg_convex)
      for (i in 1:n_seg_convex) {
        dx <- xconvex[i + 1] - xconvex[i]
        dy <- yconvex[i + 1] - yconvex[i]
        tan1[i] <- dy / dx
        int1[i] <- dy / dx * (-xconvex[i]) + yconvex[i]
      }

      ## --- concave part: secant (chord) lower envelope ---
      concave <- function(x) eval(parse(text = ccformula))
      n_seg_concave <- length(allpt) - 1
      tan2_raw <- numeric(n_seg_concave)
      int2_raw <- numeric(n_seg_concave)
      for (i in 1:n_seg_concave) {
        dx <- allpt[i + 1] - allpt[i]
        tan2_raw[i] <- (concave(allpt[i + 1]) - concave(allpt[i])) / dx
        int2_raw[i] <- -tan2_raw[i] * allpt[i] + concave(allpt[i])
      }

      ## replicate concave segments to align with convex segments
      xconcave <- rep(allpt[1:(length(allpt) - 1)], each = 2)
      yconcave <- concave(xconcave)
      tan2 <- rep(tan2_raw, each = 2)
      int2 <- rep(int2_raw, each = 2)

      ## --- integrate each combined envelope segment ---
      IntSum <- numeric(length(tan2))
      for (i in 1:length(tan2)) {
        lo <- xconvex[i]
        hi <- xconvex[i + 1]
        if (hi <= lo) {
          IntSum[i] <- 0
          next
        }
        slope <- tan1[i] + tan2[i]
        intercept <- int1[i] + int2[i]
        fun <- function(x) exp(-(slope * x + intercept))
        IntSum[i] <- integrate(fun, lo, hi)$value
      }

      ## --- inverse-CDF sampling ---
      total_int <- sum(IntSum)
      if (total_int <= 0) break
      cum <- c(0, cumsum(IntSum / total_int))
      rdm <- runif(1)
      idx <- which(rdm < cumsum(IntSum / total_int))[1]
      if (is.na(idx)) idx <- 1

      slope <- tan1[idx] + tan2[idx]
      intercept <- int1[idx] + int2[idx]

      if (abs(slope) < 1e-12) {
        ## degenerate (flat) segment: uniform sampling
        x_star <- runif(1, xconvex[idx], xconvex[idx + 1])
      } else {
        val <- (rdm - cum[idx]) * (-slope) * total_int *
               exp(intercept) + exp(-slope * xconvex[idx])
        if (val > 0) {
          x_star <- log(val) / (-slope)
        } else {
          ## fallback: uniform in segment
          x_star <- runif(1, xconvex[idx], xconvex[idx + 1])
        }
      }

      ## guard against out-of-domain
      if (!is.finite(x_star) || x_star < min || x_star > max) {
        x_star <- runif(1, xconvex[idx], xconvex[idx + 1])
      }

      ## accept / reject
      u <- runif(1)
      prop <- p(x_star) / exp(-(slope * x_star + intercept))

      support <- sort(c(x_star, support))
    }
    x_final[k] <- x_star
  }
  x_final
}
