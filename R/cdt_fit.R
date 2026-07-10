.cdt_cache <- new.env(parent = emptyenv())

.cdt_has_cmdstan <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  ver <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                  error = function(e) NULL)
  !is.null(ver)
}

.cdt_backend <- function(backend = c("auto", "cmdstanr", "rstan")) {
  backend <- match.arg(backend)
  has_cmd <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rst <- requireNamespace("rstan", quietly = TRUE)
  if (backend == "cmdstanr" && !has_cmd) stop("cmdstanr is not installed.", call. = FALSE)
  if (backend == "rstan" && !has_rst) stop("rstan is not installed.", call. = FALSE)
  if (backend == "auto") {
    # The cmdstanr namespace can load without a CmdStan installation behind it,
    # so 'auto' must confirm CmdStan before preferring that backend.
    backend <- if (.cdt_has_cmdstan()) "cmdstanr" else if (has_rst) "rstan" else
      stop("Install 'cmdstanr' with CmdStan or install 'rstan' to fit models.", call. = FALSE)
  }
  backend
}

.cdt_model <- function(backend, stan_name = "cdtmbnma.stan") {
  key <- paste("model", backend, stan_name, sep = "_")
  if (!is.null(.cdt_cache[[key]])) return(.cdt_cache[[key]])
  stan_file <- system.file("stan", stan_name, package = "cdtmbnma")
  if (stan_file == "") stop("Stan model file not found in the installed package: ", stan_name, call. = FALSE)
  mod <- if (backend == "cmdstanr") {
    cmdstanr::cmdstan_model(stan_file)
  } else {
    rstan::stan_model(file = stan_file)
  }
  .cdt_cache[[key]] <- mod
  mod
}

.cdt_validate_priors <- function(priors, defaults) {
  if (!is.list(priors)) stop("'priors' must be a named list.", call. = FALSE)
  priors <- utils::modifyList(defaults, priors)
  need <- names(defaults)
  miss <- setdiff(need, names(priors))
  if (length(miss)) stop("Missing prior value(s): ", paste(miss, collapse = ", "), call. = FALSE)
  priors <- priors[need]
  priors[] <- lapply(priors, function(x) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) stop("Prior values must be finite scalars.", call. = FALSE)
    as.numeric(x)
  })
  priors
}

.cdt_sample <- function(mod, standata, backend, chains, iter_warmup, iter_sampling,
                        adapt_delta, seed, refresh, ...) {
  if (backend == "cmdstanr") {
    mod$sample(data = standata, chains = chains,
               parallel_chains = chains,
               iter_warmup = iter_warmup, iter_sampling = iter_sampling,
               adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
  } else {
    rstan::sampling(mod, data = standata, chains = chains,
                    warmup = iter_warmup, iter = iter_warmup + iter_sampling,
                    control = list(adapt_delta = adapt_delta),
                    seed = seed, refresh = refresh, ...)
  }
}

.cdt_as_draws <- function(fit, backend) {
  posterior::as_draws_df(if (backend == "cmdstanr") fit$draws() else as.array(fit))
}

#' Fit a component dose-response network meta-analysis
#'
#' Compiles the Stan model on first use and samples the posterior. The component
#' and study structure comes from [cdt_data()]. The interaction surface is chosen
#' here.
#'
#' @param data A [cdt_data()] object.
#' @param interaction Interaction surface: `"bilinear"`, `"none"`, or
#'   `"gpdi"`. `"none"` is the additive component model. `"bilinear"` uses
#'   one parameter per co-occurring component pair. `"gpdi"` uses a saturating
#'   general pharmacodynamic interaction surface.
#' @param priors A named list from [cdt_priors()].
#' @param newdata Optional data frame of component-dose combinations at which
#'   Stan should calculate posterior predictions in generated quantities. The
#'   exported [predict.cdtmbnma()] method can also calculate predictions after
#'   fitting.
#' @param backend One of `"auto"`, `"cmdstanr"`, or `"rstan"`.
#' @param chains,iter_warmup,iter_sampling Sampler settings.
#' @param adapt_delta Target acceptance probability passed to the Stan backend.
#' @param seed Random seed.
#' @param refresh Console refresh interval; `0` silences progress.
#' @param ... Passed to the backend sampler.
#'
#' @return An object of class `cdtmbnma`.
#' @export
cdt_fit <- function(data, interaction = c("bilinear", "none", "gpdi"),
                    priors = cdt_priors(), newdata = NULL,
                    backend = "auto", chains = 4,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, seed = 1, refresh = 0, ...) {
  if (!inherits(data, "cdt_data")) stop("'data' must be a cdt_data object.", call. = FALSE)
  interaction <- match.arg(interaction)
  if (data$standata$P == 0L && interaction != "none") {
    warning("No component pairs co-occur, so the requested interaction is not identifiable; fitting the additive model.", call. = FALSE)
    interaction <- "none"
  }
  icode <- c(none = 0L, bilinear = 1L, gpdi = 2L)[[interaction]]

  C <- data$standata$C
  if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.", call. = FALSE)
    miss <- setdiff(data$components, names(newdata))
    if (length(miss)) stop("newdata is missing component column(s): ", paste(miss, collapse = ", "), call. = FALSE)
    Dp <- as.matrix(newdata[, data$components, drop = FALSE])
    storage.mode(Dp) <- "double"
    if (anyNA(Dp) || any(Dp < 0)) stop("newdata component doses must be non-negative and non-missing.", call. = FALSE)
    G <- nrow(Dp)
    if (G == 0L) stop("'newdata' must contain at least one row.", call. = FALSE)
  } else {
    Dp <- matrix(0, 1, C)
    G <- 1L
  }

  priors <- .cdt_validate_priors(priors, cdt_priors())
  stan_base <- data$standata
  if (stan_base$P == 0L) {
    stan_base$P <- 1L
    stan_base$pair_a <- as.array(1L)
    stan_base$pair_b <- as.array(1L)
  }
  standata <- c(stan_base, priors,
                list(interaction = as.integer(icode), G = as.integer(G), Dpred = Dp))

  be <- .cdt_backend(backend)
  mod <- .cdt_model(be, "cdtmbnma.stan")
  fit <- .cdt_sample(mod, standata, be, chains, iter_warmup, iter_sampling,
                     adapt_delta, seed, refresh, ...)

  draws <- .cdt_as_draws(fit, be)
  structure(list(
    draws = draws,
    fit = fit,
    data = data,
    backend = be,
    spec = list(interaction = interaction, priors = priors),
    newdata = newdata
  ), class = "cdtmbnma")
}

#' Component dose-response network meta-analysis: one-call interface
#'
#' Convenience wrapper that builds the design with [cdt_data()] and fits it with
#' [cdt_fit()].
#'
#' @inheritParams cdt_data
#' @param interaction,priors,newdata,backend,chains,iter_warmup,iter_sampling,adapt_delta,seed,refresh,... Passed to [cdt_fit()].
#' @return An object of class `cdtmbnma`.
#' @export
cdtmbnma <- function(data, study, components,
                     outcome = c("continuous", "binary"),
                     y = NULL, se = NULL, sd = NULL, n = NULL,
                     events = NULL, n_binary = NULL,
                     ref = NULL, dstar = NULL, interactions = NULL,
                     interaction = c("bilinear", "none", "gpdi"),
                     priors = cdt_priors(), newdata = NULL, backend = "auto",
                     chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                     adapt_delta = 0.95, seed = 1, refresh = 0, ...) {
  d <- cdt_data(data, study = study, components = components,
                outcome = match.arg(outcome), y = y, se = se, sd = sd, n = n,
                events = events, n_binary = n_binary, ref = ref,
                dstar = dstar, interactions = interactions)
  cdt_fit(d, interaction = match.arg(interaction), priors = priors,
          newdata = newdata, backend = backend, chains = chains,
          iter_warmup = iter_warmup, iter_sampling = iter_sampling,
          adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
}
