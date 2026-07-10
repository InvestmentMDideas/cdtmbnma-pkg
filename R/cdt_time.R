#' Assemble a two-component longitudinal component-dose network
#'
#' Builds data for the stage-one component-dose-time model. The input is a long
#' data frame with one row per study arm per follow-up time. The model currently
#' supports two components, a continuous normal likelihood with known standard
#' errors, an exponential time-course, additive component effects on log onset
#' rate, and a bilinear dose-dependent interaction on the long-term asymptote.
#'
#' @param data A data frame with one row per study-arm-time observation.
#' @param study Name of the study column.
#' @param arm Name of the arm column.
#' @param time Name of the follow-up time column, on a consistent positive scale
#'   such as weeks.
#' @param components Character vector of exactly two component dose columns.
#' @param y,se Outcome mean and standard error columns.
#' @param sd,n Optional standard deviation and sample-size columns, used as
#'   `se = sd / sqrt(n)` when `se` is absent.
#' @param ref Optional row-level or arm-level column flagging the within-study
#'   reference arm. If supplied, every row of an arm may carry the same value, or
#'   at least one row of the reference arm may be `TRUE`; exactly one arm per
#'   study must be flagged. When absent, the all-zero-dose arm is used, or the
#'   lowest-total-dose arm if no all-zero arm exists.
#' @param dstar Optional named numeric vector of dose-normalisation references
#'   for the bilinear interaction.
#'
#' @return An object of class `cdt_time_data`.
#' @export
cdt_time_data <- function(data, study, arm, time, components,
                          y, se = NULL, sd = NULL, n = NULL,
                          ref = NULL, dstar = NULL) {
  .cdt_require_data_frame(data)
  if (!is.character(components) || length(components) != 2L) {
    stop("'components' must name exactly two component dose columns for cdt_time_data().", call. = FALSE)
  }
  if (anyDuplicated(components)) stop("'components' contains duplicate names.", call. = FALSE)
  .cdt_require_columns(data, c(study, arm, time, components, y), "data")

  tvec <- as.double(data[[time]])
  ybar <- as.double(data[[y]])
  if (!is.null(se) && se %in% names(data)) {
    sevec <- as.double(data[[se]])
  } else if (!is.null(sd) && !is.null(n) && all(c(sd, n) %in% names(data))) {
    nn <- as.double(data[[n]])
    if (anyNA(nn) || any(nn <= 0)) stop("Sample sizes in 'n' must be positive.", call. = FALSE)
    sevec <- as.double(data[[sd]]) / sqrt(nn)
  } else {
    stop("Provide 'se', or both 'sd' and 'n'.", call. = FALSE)
  }
  if (anyNA(tvec) || any(tvec < 0)) stop("Follow-up times must be non-negative and non-missing.", call. = FALSE)
  if (anyNA(ybar) || anyNA(sevec) || any(sevec <= 0)) stop("Outcome means must be non-missing and standard errors positive.", call. = FALSE)

  Drow <- as.matrix(data[, components, drop = FALSE])
  storage.mode(Drow) <- "double"
  if (anyNA(Drow) || any(Drow < 0)) stop("Component doses must be non-negative and non-missing.", call. = FALSE)

  study_chr <- as.character(data[[study]])
  study_levels <- unique(study_chr)
  study_fac <- factor(study_chr, levels = study_levels)
  study_idx_row <- as.integer(study_fac)

  arm_key <- paste(study_chr, as.character(data[[arm]]), sep = "\r")
  arm_levels <- unique(arm_key)
  arm_fac <- factor(arm_key, levels = arm_levels)
  row_arm <- as.integer(arm_fac)
  n_arms <- length(arm_levels)
  n_obs <- nrow(data)
  n_studies <- length(study_levels)

  arm_sid <- integer(n_arms)
  arm_dose <- matrix(NA_real_, n_arms, 2L)
  arm_label <- character(n_arms)
  for (a in seq_len(n_arms)) {
    rows <- which(row_arm == a)
    arm_sid[a] <- study_idx_row[rows[1]]
    arm_label[a] <- as.character(data[[arm]][rows[1]])
    u <- unique(as.data.frame(Drow[rows, , drop = FALSE]))
    if (nrow(u) != 1L) stop("Component doses must be constant within each study arm.", call. = FALSE)
    arm_dose[a, ] <- as.numeric(u[1, ])
  }
  colnames(arm_dose) <- components

  arm_df <- data.frame(.study = study_levels[arm_sid], arm = arm_label,
                       arm_dose, check.names = FALSE)
  if (!is.null(ref) && ref %in% names(data)) {
    raw_ref <- as.logical(data[[ref]])
    if (anyNA(raw_ref)) stop("Reference column contains values that cannot be interpreted as TRUE/FALSE.", call. = FALSE)
    arm_ref <- vapply(seq_len(n_arms), function(a) any(raw_ref[row_arm == a]), logical(1))
    arm_df$.ref <- arm_ref
    is_ref <- .cdt_reference_arms(arm_df, arm_sid, study_levels, arm_dose, ref = ".ref")
  } else {
    is_ref <- .cdt_reference_arms(arm_df, arm_sid, study_levels, arm_dose, ref = NULL)
  }

  active_idx <- which(is_ref == 0L)
  if (!length(active_idx)) stop("At least one non-reference arm is required.", call. = FALSE)

  if (is.null(dstar)) {
    dstar <- apply(arm_dose, 2, max)
    dstar[dstar == 0] <- 1
  } else {
    if (is.null(names(dstar))) stop("'dstar' must be a named numeric vector.", call. = FALSE)
    dstar <- as.double(dstar[components])
    if (anyNA(dstar)) stop("'dstar' must name both components.", call. = FALSE)
    if (any(dstar <= 0)) stop("Every value in 'dstar' must be positive.", call. = FALSE)
  }
  names(dstar) <- components

  standata <- list(
    n_obs = as.integer(n_obs),
    n_arms = as.integer(n_arms),
    n_studies = as.integer(n_studies),
    n_active = as.integer(length(active_idx)),
    row_arm = as.integer(row_arm),
    time = tvec,
    ybar = ybar,
    se = sevec,
    arm_sid = as.integer(arm_sid),
    arm_dA = as.array(arm_dose[, 1]),
    arm_dB = as.array(arm_dose[, 2]),
    active_idx = as.array(as.integer(active_idx)),
    dAstar = as.numeric(dstar[1]),
    dBstar = as.numeric(dstar[2])
  )

  structure(list(
    standata = standata,
    components = components,
    studies = study_levels,
    arms = arm_label,
    dstar = dstar,
    is_ref = is_ref,
    n_obs = n_obs,
    n_arms = n_arms,
    n_studies = n_studies
  ), class = "cdt_time_data")
}

#' @export
print.cdt_time_data <- function(x, ...) {
  cat("<cdt_time_data>\n")
  cat(sprintf("  observations: %d\n", x$n_obs))
  cat(sprintf("  arms:         %d across %d studies\n", x$n_arms, x$n_studies))
  cat(sprintf("  components:   %s\n", paste(x$components, collapse = ", ")))
  cat("  model:        two-component exponential time-course with bilinear interaction\n")
  invisible(x)
}

#' Fit the two-component component-dose-time model
#'
#' Fits the stage-one longitudinal model in `inst/stan/cdt_stage1.stan`. This is
#' a two-component model for repeated arm means over follow-up time. The long-term
#' asymptote carries the component dose-response curves and bilinear interaction;
#' the onset rate carries additive component dose effects.
#'
#' @param data A [cdt_time_data()] object.
#' @param priors A named list from [cdt_time_priors()].
#' @param newdata Optional data frame of two-component dose combinations for
#'   generated-quantities predictions of the long-term relative effect.
#' @param backend One of `"auto"`, `"cmdstanr"`, or `"rstan"`.
#' @param chains,iter_warmup,iter_sampling Sampler settings.
#' @param adapt_delta Target acceptance probability.
#' @param seed Random seed.
#' @param refresh Console refresh interval; `0` silences progress.
#' @param ... Passed to the backend sampler.
#'
#' @return An object of class `cdt_timefit`.
#' @export
cdt_time_fit <- function(data, priors = cdt_time_priors(), newdata = NULL,
                         backend = "auto", chains = 4,
                         iter_warmup = 1000, iter_sampling = 1000,
                         adapt_delta = 0.95, seed = 1, refresh = 0, ...) {
  if (!inherits(data, "cdt_time_data")) stop("'data' must be a cdt_time_data object.", call. = FALSE)
  comp <- data$components
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.", call. = FALSE)
    miss <- setdiff(comp, names(newdata))
    if (length(miss)) stop("newdata is missing component column(s): ", paste(miss, collapse = ", "), call. = FALSE)
    Dp <- as.matrix(newdata[, comp, drop = FALSE])
    storage.mode(Dp) <- "double"
    if (anyNA(Dp) || any(Dp < 0)) stop("newdata component doses must be non-negative and non-missing.", call. = FALSE)
    n_pred <- nrow(Dp)
    if (n_pred == 0L) stop("'newdata' must contain at least one row.", call. = FALSE)
  } else {
    Dp <- matrix(0, 1, 2)
    n_pred <- 1L
  }

  priors <- .cdt_validate_priors(priors, cdt_time_priors())
  standata <- c(data$standata, priors,
                list(n_pred = as.integer(n_pred),
                     pred_dA = as.array(Dp[, 1]),
                     pred_dB = as.array(Dp[, 2])))

  be <- .cdt_backend(backend)
  mod <- .cdt_model(be, "cdt_stage1.stan")
  fit <- .cdt_sample(mod, standata, be, chains, iter_warmup, iter_sampling,
                     adapt_delta, seed, refresh, ...)
  draws <- .cdt_as_draws(fit, be)
  structure(list(draws = draws, fit = fit, data = data, backend = be,
                 spec = list(priors = priors), newdata = newdata),
            class = "cdt_timefit")
}
