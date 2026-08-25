#' Modified Adaptive Rejection Sampling Algorithm
#'
#' rMARS generates a sequence of random numbers using the modified adaptive
#' rejection sampling algorithm, which extends ARS to non-log-concave densities.
#'
#' @param n Desired sample size (positive integer).
#' @param formula Kernel of the target distribution, as a character string.
#' @param min,max Domain of the target distribution, including \code{-Inf} and \code{Inf}.
#' @param sp Supporting set; a numeric vector of initial abscissae.
#' @param infp Inflexion set; a numeric vector of inflexion points of \code{-log(p(x))}.
#' @param m A small numeric parameter for judging concavity and convexity near inflexion points.
#' @return A numeric vector of length \code{n} containing samples from the target distribution.
#' @author Dong Zhang
#' @references Martino L, Miguez J. A generalization of the adaptive rejection sampling algorithm[J]. Statistics & Computing, 2011, 21(4):633-647.
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rMARS(100, "exp(-(4-x^2)^2)", -Inf, Inf,
#'            c(-2.5, 0, 2.5), c(-2/sqrt(3), 2/sqrt(3)))
#' hist(x, probability = TRUE, xlim = c(-3, 3), ylim = c(0, 1.2), breaks = 20)
#' lines(density(x, bw = 0.05), col = "blue")
#'
rMARS <- function(n, formula, min = -Inf, max = Inf, sp, infp, m = 1e-4) {
  ## --- input validation ---
  sp <- sort(sp)
  infp <- sort(infp)
  if (!is.character(formula))
    stop("'formula' must be a character string.")
  if (n <= 0)
    stop("Sample size 'n' must be a positive integer.")
  if (min >= max)
    stop("'min' must be strictly less than 'max'.")
  if (length(sp) < 1)
    stop("Supporting set 'sp' must contain at least 1 point.")
  if (length(infp) < 1)
    stop("Inflexion set 'infp' must contain at least 1 point.")

  p <- function(x) eval(parse(text = formula))
  V <- function(x) {
    px <- p(x)
    ifelse(px > 0, -log(px), Inf)
  }

  ## --- symbolic derivatives (computed once) ---
  neglog_formula <- paste0("-log(", formula, ")")
  deriv1 <- function(x) eval(D(parse(text = neglog_formula), "x"))
  deriv2 <- function(x) eval(D(D(parse(text = neglog_formula), "x"), "x"))

  ## --- internal helpers for tangent (tu) / secant (ao) envelope segments ---
  ## Each returns list(tg, int, crp, crv) for a given row index x in pandl.

  ## tangent-based envelope: cross-points from adjacent tangent lines
  make_tangent_seg <- function(usepoint, n_pts) {
    tg <- deriv1(usepoint)
    int <- V(usepoint) - tg * usepoint
    crp <- numeric(n_pts)
    crv <- numeric(n_pts)
    for (i in 1:n_pts) {
      A <- matrix(c(tg[i], -1, tg[i + 1], -1), nrow = 2, byrow = TRUE)
      b <- -c(int[i], int[i + 1])
      sol <- solve(A, b)
      crp[i] <- sol[1]
      crv[i] <- sol[2]
    }
    list(tg = tg, int = int, crp = crp, crv = crv)
  }

  ## secant-based envelope: linear interpolation between endpoints
  make_secant_seg <- function(usepoint, n_pts, drop_last = FALSE) {
    crp <- usepoint
    crv <- V(usepoint)
    tg <- numeric(n_pts)
    int <- numeric(n_pts)
    for (i in 1:n_pts) {
      tg[i] <- (V(usepoint[i + 1]) - V(usepoint[i])) /
               (usepoint[i + 1] - usepoint[i])
      int[i] <- V(usepoint[i]) - tg[i] * usepoint[i]
    }
    if (drop_last) {
      crp <- crp[-length(crp)]
      crv <- crv[-length(crv)]
    }
    list(tg = tg, int = int, crp = crp, crv = crv)
  }

  ## dispatch by region type
  envelope <- function(x) {
    rtype <- pandl[x, 1][[1]]
    pts <- pandl[x, 2][[1]]
    n_t <- pandl[x, 3][[1]]
    n_s <- pandl[x, 4][[1]]

    switch(rtype,
      "ltuInf" = make_tangent_seg(c(pts, infp[1]), n_t),
      "ltuFin" = make_tangent_seg(c(min, pts, infp[1]), n_t),
      "laoFin" = {
        r <- make_secant_seg(c(min, pts, infp[1]), n_t)
        r$crp <- r$crp[-1]
        r$crv <- r$crv[-1]
        r
      },
      "rtuInf" = make_tangent_seg(c(tail(infp, 1), pts), n_t),
      "rtuFin" = make_tangent_seg(c(tail(infp, 1), pts, max), n_t),
      "raoFin" = {
        r <- make_secant_seg(c(tail(infp, 1), pts, max), n_t)
        r$crp <- r$crp[-length(r$crp)]
        r$crv <- r$crv[-length(r$crv)]
        r
      },
      "tu" = make_tangent_seg(c(infp[x - 1], pts, infp[x]), n_t),
      "ao" = make_secant_seg(c(infp[x - 1], pts, infp[x]), n_s),
      stop("Unknown region type: ", rtype)
    )
  }

  x_final <- numeric(n)

  for (N in 1:n) {
    u <- 0
    rate <- -1

    while (u >= rate) {
      ## --- classify each region as tangent (concave) or secant (convex) ---
      corc <- character(length(infp) + 1)

      ## right boundary
      d2_right <- deriv2(infp[length(infp)] + m)
      if (is.finite(d2_right) && d2_right > 0) {
        corc[length(infp) + 1] <- if (is.infinite(max)) "rtuInf" else "rtuFin"
      } else if (is.finite(d2_right) && d2_right < 0) {
        corc[length(infp) + 1] <- if (!is.infinite(max)) "raoFin"
                                  else "rtuInf"
      } else {
        ## d2 == 0 or non-finite: default to tangent as fallback
        corc[length(infp) + 1] <- if (is.infinite(max)) "rtuInf" else "rtuFin"
      }

      ## left boundary
      d2_left <- deriv2(infp[1] - m)
      if (is.finite(d2_left) && d2_left > 0) {
        corc[1] <- if (is.infinite(min)) "ltuInf" else "ltuFin"
      } else if (is.finite(d2_left) && d2_left < 0) {
        corc[1] <- if (!is.infinite(min)) "laoFin" else "ltuInf"
      } else {
        corc[1] <- if (is.infinite(min)) "ltuInf" else "ltuFin"
      }

      ## interior regions
      if (length(infp) > 1) {
        for (i in 2:length(infp)) {
          d2_mid <- deriv2(infp[i - 1] + m)
          if (is.finite(d2_mid) && d2_mid < 0) {
            corc[i] <- "ao"
          } else {
            corc[i] <- "tu"
          }
        }
      }

      ## --- partition support points by inflexion intervals ---
      parsp <- vector("list", length(infp) + 1)
      parsp[[1]] <- sp[sp < infp[1]]
      parsp[[length(infp) + 1]] <- sp[sp > infp[length(infp)]]
      if (length(infp) > 1) {
        for (i in 2:length(infp))
          parsp[[i]] <- sp[sp > infp[i - 1] & sp < infp[i]]
      }

      ## --- build pandl table: type, points, tangent-count, secant-count ---
      pandl <- cbind(corc, parsp,
                     pt = numeric(length(corc)),
                     lns = numeric(length(corc)))
      for (i in 1:nrow(pandl)) {
        npts <- length(pandl[i, 2][[1]])
        rtype <- pandl[i, 1][[1]]
        if (rtype == "ao") {
          pandl[i, 3][[1]] <- npts
          pandl[i, 4][[1]] <- npts + 1
        } else if (rtype == "tu") {
          pandl[i, 3][[1]] <- npts + 1
          pandl[i, 4][[1]] <- npts + 2
        } else if (rtype %in% c("ltuInf", "rtuInf")) {
          pandl[i, 3][[1]] <- npts
          pandl[i, 4][[1]] <- npts + 1
        } else if (rtype %in% c("ltuFin", "rtuFin")) {
          pandl[i, 3][[1]] <- npts + 1
          pandl[i, 4][[1]] <- npts + 2
        } else if (rtype %in% c("laoFin", "raoFin")) {
          pandl[i, 3][[1]] <- npts + 1
          pandl[i, 4][[1]] <- npts + 1
        }
      }

      ## --- compute envelope segments for all regions ---
      tg_total <- c()
      int_total <- c()
      crp_total <- c()
      crv_total <- c()

      for (i in 1:nrow(pandl)) {
        tsf <- envelope(i)
        tg_total <- c(tg_total, tsf$tg)
        int_total <- c(int_total, tsf$int)
        crp_total <- c(crp_total, tsf$crp)
        crv_total <- c(crv_total, tsf$crv)
      }

      ## --- assemble sorted (x, V(x)) table for integration limits ---
      fdtfram <- rbind(cbind(crp_total, crv_total),
                       c(min, V(min)), c(max, V(max)))
      fdtfram <- rbind(fdtfram,
                       matrix(c(infp, V(infp)), nrow = length(infp), byrow = FALSE))
      fdtfram <- fdtfram[order(fdtfram[, 1]), ]

      ## --- integrate each envelope segment ---
      intsum <- numeric(length(tg_total))
      for (i in 1:length(tg_total)) {
        lo <- fdtfram[i, 1]
        hi <- fdtfram[i + 1, 1]
        if (is.finite(lo) && is.finite(hi) && hi > lo) {
          integ <- function(x) exp(-(tg_total[i] * x + int_total[i]))
          intsum[i] <- integrate(integ, lo, hi)$value
        } else if (is.finite(hi) && (!is.finite(lo) || lo < hi)) {
          ## unbounded left
          integ <- function(x) exp(-(tg_total[i] * x + int_total[i]))
          intsum[i] <- integrate(integ, -Inf, hi)$value
        } else if (is.finite(lo) && (!is.finite(hi) || hi > lo)) {
          ## unbounded right
          integ <- function(x) exp(-(tg_total[i] * x + int_total[i]))
          intsum[i] <- integrate(integ, lo, Inf)$value
        } else {
          intsum[i] <- 0
        }
      }

      ## --- inverse-CDF sampling ---
      rdm <- runif(1)
      total_int <- sum(intsum)
      if (total_int <= 0) break
      cum <- c(0, cumsum(intsum / total_int))
      idx <- which(rdm < cumsum(intsum / total_int))[1]
      if (is.na(idx)) idx <- 1

      t_idx <- tg_total[idx]
      ## proper if/else instead of ifelse
      if (idx > 1) {
        x_star <- (log(-(rdm - cum[idx]) * total_int * t_idx +
                       exp(-t_idx * fdtfram[idx, 1] - int_total[idx])) +
                   int_total[idx]) / (-t_idx)
      } else {
        x_star <- (log(-rdm * total_int * t_idx +
                       exp(-t_idx * fdtfram[2, 1] - int_total[1])) +
                   int_total[1]) / (-t_idx)
      }

      ## guard against NaN / Inf from degenerate tangent
      if (!is.finite(x_star)) {
        x_star <- runif(1, fdtfram[idx, 1], fdtfram[idx + 1, 1])
      }

      ## accept / reject
      u <- runif(1)
      rate <- p(x_star) / exp(-t_idx * x_star - int_total[idx])

      sp <- sort(c(sp, x_star))
    }
    x_final[N] <- x_star
  }
  x_final
}
