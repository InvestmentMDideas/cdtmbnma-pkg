## Internal: pull a named parameter block as a draws-by-index matrix ----------
.draw_mat <- function(draws, base, k) {
  if (k == 0L) return(matrix(0, nrow(draws), 0))
  cols <- sprintf("%s[%d]", base, seq_len(k))
  miss <- setdiff(cols, names(draws))
  if (length(miss)) stop("Fitted draws are missing parameter column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  as.matrix(draws[, cols, drop = FALSE])
}

.scalar_draw <- function(draws, name) {
  if (!name %in% names(draws)) stop("Fitted draws are missing parameter column: ", name, call. = FALSE)
  as.numeric(draws[[name]])
}

#' Predict relative effects at component-dose combinations
#'
#' Computes the posterior distribution of the relative effect against the
#' all-zero reference at one or more component-dose combinations, including dose
#' pairs or combinations that no trial observed. The random effect is marginalised
#' at zero, so the prediction is the population-average structural effect rather
#' than a new-study predictive interval. The effect is on the outcome scale for a
#' continuous outcome and on the log-odds scale for a binary outcome.
#'
#' @param object A fitted `cdtmbnma` object.
#' @param newdata Data frame of dose combinations, carrying the component columns
#'   used in fitting.
#' @param probs Two or more quantiles for credible intervals.
#' @param ... Unused.
#'
#' @return A data frame with the supplied doses, posterior mean, posterior
#'   standard deviation, and requested quantiles for each row of `newdata`.
#' @export
predict.cdtmbnma <- function(object, newdata, probs = c(0.025, 0.975), ...) {
  if (!inherits(object, "cdtmbnma")) stop("'object' must be a cdtmbnma fit.", call. = FALSE)
  if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.", call. = FALSE)
  if (!length(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1)) stop("'probs' must be between 0 and 1.", call. = FALSE)

  d <- object$data
  comp <- d$components
  C <- length(comp)
  miss <- setdiff(comp, names(newdata))
  if (length(miss)) stop("newdata is missing component column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  Dp <- as.matrix(newdata[, comp, drop = FALSE])
  storage.mode(Dp) <- "double"
  if (anyNA(Dp) || any(Dp < 0)) stop("newdata component doses must be non-negative and non-missing.", call. = FALSE)
  if (nrow(Dp) == 0L) stop("'newdata' must contain at least one row.", call. = FALSE)

  draws <- object$draws
  emax <- .draw_mat(draws, "emax", C)
  ED50 <- .draw_mat(draws, "ED50", C)
  itype <- object$spec$interaction
  P <- nrow(d$pairs)
  dstar <- d$dstar
  eta <- if (itype == "bilinear") .draw_mat(draws, "eta", P) else NULL
  INT <- if (itype == "gpdi") .draw_mat(draws, "INT", P) else NULL
  kappa <- if (itype == "gpdi") .draw_mat(draws, "kappa", C) else NULL

  out <- vector("list", nrow(Dp))
  for (g in seq_len(nrow(Dp))) {
    dose <- Dp[g, ]
    lp <- rep(0, nrow(draws))
    for (c in seq_len(C)) {
      lp <- lp + emax[, c] * dose[c] / (ED50[, c] + dose[c])
    }
    if (itype == "bilinear" && P > 0L) {
      for (p in seq_len(P)) {
        a <- d$pairs[p, 1]
        b <- d$pairs[p, 2]
        lp <- lp + eta[, p] * (dose[a] / dstar[a]) * (dose[b] / dstar[b])
      }
    } else if (itype == "gpdi" && P > 0L) {
      for (p in seq_len(P)) {
        a <- d$pairs[p, 1]
        b <- d$pairs[p, 2]
        lp <- lp + INT[, p] * (dose[a] / (kappa[, a] + dose[a])) *
          (dose[b] / (kappa[, b] + dose[b]))
      }
    }
    qs <- stats::quantile(lp, probs = probs, names = TRUE)
    out[[g]] <- c(mean = mean(lp), sd = stats::sd(lp), qs)
  }
  res <- as.data.frame(do.call(rbind, out), check.names = FALSE)
  cbind(newdata[, comp, drop = FALSE], res, row.names = NULL)
}

#' Predict long-term effects from a two-component time-course fit
#'
#' Computes the posterior distribution of the long-term, asymptotic relative
#' effect against the all-zero reference for the stage-one longitudinal model.
#'
#' @param object A fitted `cdt_timefit` object.
#' @param newdata Data frame carrying the two component-dose columns used in
#'   fitting.
#' @param probs Quantiles for credible intervals.
#' @param ... Unused.
#'
#' @return A data frame with the supplied doses, posterior mean, posterior
#'   standard deviation, and requested quantiles.
#' @export
predict.cdt_timefit <- function(object, newdata, probs = c(0.025, 0.975), ...) {
  if (!inherits(object, "cdt_timefit")) stop("'object' must be a cdt_timefit fit.", call. = FALSE)
  if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.", call. = FALSE)
  comp <- object$data$components
  miss <- setdiff(comp, names(newdata))
  if (length(miss)) stop("newdata is missing component column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  Dp <- as.matrix(newdata[, comp, drop = FALSE])
  storage.mode(Dp) <- "double"
  if (anyNA(Dp) || any(Dp < 0)) stop("newdata component doses must be non-negative and non-missing.", call. = FALSE)
  if (nrow(Dp) == 0L) stop("'newdata' must contain at least one row.", call. = FALSE)

  draws <- object$draws
  a1 <- .scalar_draw(draws, "a1_EmaxA")
  a2 <- .scalar_draw(draws, "a2_EmaxB")
  ED50A <- .scalar_draw(draws, "ED50A")
  ED50B <- .scalar_draw(draws, "ED50B")
  eta <- .scalar_draw(draws, "eta")
  dAstar <- object$data$dstar[[1]]
  dBstar <- object$data$dstar[[2]]

  out <- vector("list", nrow(Dp))
  for (g in seq_len(nrow(Dp))) {
    dA <- Dp[g, 1]
    dB <- Dp[g, 2]
    eff <- a1 * dA / (ED50A + dA) +
      a2 * dB / (ED50B + dB) +
      eta * (dA / dAstar) * (dB / dBstar)
    qs <- stats::quantile(eff, probs = probs, names = TRUE)
    out[[g]] <- c(mean = mean(eff), sd = stats::sd(eff), qs)
  }
  res <- as.data.frame(do.call(rbind, out), check.names = FALSE)
  cbind(newdata[, comp, drop = FALSE], res, row.names = NULL)
}
