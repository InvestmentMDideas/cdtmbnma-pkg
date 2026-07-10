#' Prior settings for single-timepoint cdtmbnma models
#'
#' Weakly informative defaults on the model's structural parameters. Effects
#' live on the outcome scale for continuous outcomes and on the log-odds scale
#' for binary outcomes, so scale the component-effect and interaction priors to
#' the outcome.
#'
#' @param emax_sd Prior standard deviation for each component's Emax, the
#'   maximal effect.
#' @param logED50_mean,logED50_sd Prior mean and standard deviation for each
#'   component's log half-maximal dose. Also used for the saturating-surface
#'   half-doses in the `gpdi` model.
#' @param int_sd Prior standard deviation for the interaction parameter, either
#'   bilinear `eta` or saturating `INT`, centred at additivity.
#' @param ref_sd Prior standard deviation for the free per-study reference level.
#' @param omega_sd Scale of the half-normal prior on the random-effect standard
#'   deviation.
#'
#' @return A named list of prior hyperparameters.
#' @export
cdt_priors <- function(emax_sd = 10, logED50_mean = log(50), logED50_sd = 1,
                       int_sd = 5, ref_sd = 10, omega_sd = 2) {
  vals <- list(prior_emax_sd = emax_sd,
               prior_logED50_mean = logED50_mean,
               prior_logED50_sd = logED50_sd,
               prior_int_sd = int_sd,
               prior_ref_sd = ref_sd,
               prior_omega_sd = omega_sd)
  for (nm in names(vals)) {
    x <- vals[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) stop("Prior values must be finite numeric scalars.", call. = FALSE)
  }
  if (emax_sd <= 0 || logED50_sd <= 0 || int_sd <= 0 || ref_sd <= 0 || omega_sd <= 0) {
    stop("Prior scale parameters must be positive.", call. = FALSE)
  }
  vals
}

#' Prior settings for the two-component time-course model
#'
#' These priors are passed to the stage-one longitudinal Stan model. The model
#' places an exponential time-course on the arm mean and a bilinear
#' dose-dependent interaction on the long-term asymptote.
#'
#' @param ref_E_sd Prior standard deviation for study-specific reference
#'   asymptotes.
#' @param ref_lograte_mean,ref_lograte_sd Prior mean and standard deviation for
#'   study-specific reference log onset rates.
#' @param emax_sd Prior standard deviation for each component's Emax.
#' @param logED50_mean,logED50_sd Prior mean and standard deviation for the log
#'   half-maximal dose parameters.
#' @param eta_sd Prior standard deviation for the bilinear interaction.
#' @param rate_sd Prior standard deviation for component log-rate dose effects.
#' @param omega_E_sd,omega_k_sd Half-normal prior scales for random-effect
#'   standard deviations on the asymptote and log-rate.
#'
#' @return A named list of prior hyperparameters.
#' @export
cdt_time_priors <- function(ref_E_sd = 5,
                            ref_lograte_mean = -3,
                            ref_lograte_sd = 1,
                            emax_sd = 10,
                            logED50_mean = 0,
                            logED50_sd = 1,
                            eta_sd = 5,
                            rate_sd = 1,
                            omega_E_sd = 2,
                            omega_k_sd = 0.5) {
  vals <- list(prior_ref_E_sd = ref_E_sd,
               prior_ref_lograte_mean = ref_lograte_mean,
               prior_ref_lograte_sd = ref_lograte_sd,
               prior_emax_sd = emax_sd,
               prior_logED50_mean = logED50_mean,
               prior_logED50_sd = logED50_sd,
               prior_eta_sd = eta_sd,
               prior_rate_sd = rate_sd,
               prior_omega_E_sd = omega_E_sd,
               prior_omega_k_sd = omega_k_sd)
  for (nm in names(vals)) {
    x <- vals[[nm]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) stop("Prior values must be finite numeric scalars.", call. = FALSE)
  }
  scale_names <- setdiff(names(vals), c("prior_ref_lograte_mean", "prior_logED50_mean"))
  if (any(vapply(vals[scale_names], function(x) x <= 0, logical(1)))) {
    stop("Prior scale parameters must be positive.", call. = FALSE)
  }
  vals
}
