.cdt_draw_summary <- function(draws, vars) {
  s <- posterior::summarise_draws(
    draws,
    "mean", "sd",
    ~ stats::quantile(.x, probs = c(0.025, 0.975), names = TRUE),
    "rhat", "ess_bulk"
  )
  s <- s[s$variable %in% vars, , drop = FALSE]
  as.data.frame(s)
}

#' @export
print.cdtmbnma <- function(x, ...) {
  d <- x$data
  cat("<cdtmbnma fit>\n")
  cat(sprintf("  %d arms, %d studies, %d components\n",
              d$n_arms, d$n_studies, length(d$components)))
  cat(sprintf("  outcome: %s   interaction: %s   backend: %s\n",
              d$outcome, x$spec$interaction, x$backend))
  cat("\nComponent effects (posterior means and 95% credible intervals):\n")
  print(stats::coef(x), row.names = FALSE)
  invisible(x)
}

#' Posterior summary of a cdtmbnma fit
#'
#' @param object A fitted `cdtmbnma` object.
#' @param ... Unused.
#' @return A data frame summarising the structural parameters, including
#'   convergence diagnostics where available.
#' @export
summary.cdtmbnma <- function(object, ...) {
  if (!inherits(object, "cdtmbnma")) stop("'object' must be a cdtmbnma fit.", call. = FALSE)
  d <- object$data
  C <- length(d$components)
  P <- nrow(d$pairs)
  itype <- object$spec$interaction
  vars <- c(sprintf("emax[%d]", seq_len(C)),
            sprintf("ED50[%d]", seq_len(C)),
            if (itype == "bilinear" && P) sprintf("eta[%d]", seq_len(P)),
            if (itype == "gpdi" && P) sprintf("INT[%d]", seq_len(P)),
            if (itype == "gpdi") sprintf("kappa[%d]", seq_len(C)),
            "omega")
  s <- .cdt_draw_summary(object$draws, vars)

  lab <- s$variable
  for (cc in seq_len(C)) {
    lab <- sub(sprintf("^emax\\[%d\\]$", cc), paste0("emax: ", d$components[cc]), lab)
    lab <- sub(sprintf("^ED50\\[%d\\]$", cc), paste0("ED50: ", d$components[cc]), lab)
    lab <- sub(sprintf("^kappa\\[%d\\]$", cc), paste0("kappa: ", d$components[cc]), lab)
  }
  if (P) {
    for (p in seq_len(P)) {
      nm <- paste(d$components[d$pairs[p, ]], collapse = " x ")
      lab <- sub(sprintf("^eta\\[%d\\]$", p), paste0("interaction: ", nm), lab)
      lab <- sub(sprintf("^INT\\[%d\\]$", p), paste0("interaction: ", nm), lab)
    }
  }
  s$variable <- lab
  s
}

#' Extract component-effect coefficients
#'
#' @param object A fitted `cdtmbnma` object.
#' @param ... Unused.
#' @return A data frame with posterior means and 95% credible intervals for
#'   component Emax values and posterior means for ED50 values.
#' @export
coef.cdtmbnma <- function(object, ...) {
  if (!inherits(object, "cdtmbnma")) stop("'object' must be a cdtmbnma fit.", call. = FALSE)
  d <- object$data
  C <- length(d$components)
  em <- .draw_mat(object$draws, "emax", C)
  ed <- .draw_mat(object$draws, "ED50", C)
  data.frame(
    component = d$components,
    emax = colMeans(em),
    emax_q2.5 = apply(em, 2, stats::quantile, 0.025),
    emax_q97.5 = apply(em, 2, stats::quantile, 0.975),
    ED50 = colMeans(ed),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Plot component dose-response curves
#'
#' Draws each component's marginal dose-response curve, varying the named
#' component while holding the others at zero. The curve is shown with a posterior
#' credible band on the model scale.
#'
#' @param x A fitted `cdtmbnma` object.
#' @param ngrid Number of dose points per curve.
#' @param probs Two quantiles for the credible band.
#' @param ... Passed to [plot()].
#' @return Invisibly, a named list of plotted grids.
#' @export
plot.cdtmbnma <- function(x, ngrid = 60, probs = c(0.025, 0.975), ...) {
  if (!inherits(x, "cdtmbnma")) stop("'x' must be a cdtmbnma fit.", call. = FALSE)
  if (ngrid < 2L) stop("'ngrid' must be at least 2.", call. = FALSE)
  if (length(probs) != 2L || anyNA(probs) || any(probs <= 0 | probs >= 1)) stop("'probs' must contain two probabilities between 0 and 1.", call. = FALSE)
  d <- x$data
  comp <- d$components
  C <- length(comp)
  Dobs <- d$standata$D
  draws <- x$draws
  emax <- .draw_mat(draws, "emax", C)
  ED50 <- .draw_mat(draws, "ED50", C)

  nr <- ceiling(sqrt(C))
  nc <- ceiling(C / nr)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op), add = TRUE)
  grids <- vector("list", C)
  for (cc in seq_len(C)) {
    dmax <- max(Dobs[, cc])
    if (dmax == 0) dmax <- 1
    g <- seq(0, dmax, length.out = ngrid)
    M <- sapply(g, function(dose) emax[, cc] * dose / (ED50[, cc] + dose))
    m <- colMeans(M)
    lo <- apply(M, 2, stats::quantile, probs[1])
    hi <- apply(M, 2, stats::quantile, probs[2])
    graphics::plot(g, m, type = "n", ylim = range(lo, hi, 0),
                   xlab = paste0(comp[cc], " dose"), ylab = "effect",
                   main = comp[cc], ...)
    graphics::polygon(c(g, rev(g)), c(lo, rev(hi)),
                      col = grDevices::adjustcolor("steelblue", 0.2), border = NA)
    graphics::lines(g, m, col = "steelblue", lwd = 2)
    graphics::abline(h = 0, col = "grey60", lty = 2)
    grids[[cc]] <- data.frame(dose = g, mean = m, lo = lo, hi = hi)
  }
  names(grids) <- comp
  invisible(grids)
}

#' @export
print.cdt_timefit <- function(x, ...) {
  d <- x$data
  cat("<cdt_timefit>\n")
  cat(sprintf("  %d observations, %d arms, %d studies\n", d$n_obs, d$n_arms, d$n_studies))
  cat(sprintf("  components: %s\n", paste(d$components, collapse = ", ")))
  cat(sprintf("  backend:    %s\n", x$backend))
  cat("\nLong-term component effects (posterior means and 95% credible intervals):\n")
  print(stats::coef(x), row.names = FALSE)
  invisible(x)
}

#' Posterior summary of a two-component time-course fit
#'
#' @param object A fitted `cdt_timefit` object.
#' @param ... Unused.
#' @return A data frame summarising key structural parameters.
#' @export
summary.cdt_timefit <- function(object, ...) {
  if (!inherits(object, "cdt_timefit")) stop("'object' must be a cdt_timefit fit.", call. = FALSE)
  vars <- c("a1_EmaxA", "a2_EmaxB", "ED50A", "ED50B", "eta",
            "b1_rateA", "b2_rateB", "omega_E", "omega_k")
  s <- .cdt_draw_summary(object$draws, vars)
  labs <- c(a1_EmaxA = paste0("emax: ", object$data$components[1]),
            a2_EmaxB = paste0("emax: ", object$data$components[2]),
            ED50A = paste0("ED50: ", object$data$components[1]),
            ED50B = paste0("ED50: ", object$data$components[2]),
            eta = paste0("interaction: ", paste(object$data$components, collapse = " x ")),
            b1_rateA = paste0("log-rate dose effect: ", object$data$components[1]),
            b2_rateB = paste0("log-rate dose effect: ", object$data$components[2]),
            omega_E = "random-effect SD: asymptote",
            omega_k = "random-effect SD: log-rate")
  s$variable <- unname(labs[s$variable])
  s
}

#' Extract long-term coefficients from a two-component time-course fit
#'
#' @param object A fitted `cdt_timefit` object.
#' @param ... Unused.
#' @return A data frame with Emax and ED50 summaries for the two components.
#' @export
coef.cdt_timefit <- function(object, ...) {
  if (!inherits(object, "cdt_timefit")) stop("'object' must be a cdt_timefit fit.", call. = FALSE)
  em <- cbind(.scalar_draw(object$draws, "a1_EmaxA"), .scalar_draw(object$draws, "a2_EmaxB"))
  ed <- cbind(.scalar_draw(object$draws, "ED50A"), .scalar_draw(object$draws, "ED50B"))
  data.frame(
    component = object$data$components,
    emax = colMeans(em),
    emax_q2.5 = apply(em, 2, stats::quantile, 0.025),
    emax_q97.5 = apply(em, 2, stats::quantile, 0.975),
    ED50 = colMeans(ed),
    row.names = NULL,
    check.names = FALSE
  )
}
