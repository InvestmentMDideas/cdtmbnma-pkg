#' cdtmbnma: component, dose, and time network meta-analysis
#'
#' `cdtmbnma` fits Bayesian component model-based network meta-analysis models
#' in which treatment arms are decomposed into component-specific doses. The main
#' model estimates component Emax dose-response curves and pairwise
#' dose-dependent interaction surfaces for single-timepoint arm-level networks.
#' Continuous and binary outcomes are supported. A second, stage-one interface
#' fits a two-component exponential time-course model for repeated arm means.
#'
#' The usual single-timepoint workflow is [cdt_data()] to build the design, then
#' [cdt_fit()] to sample, or [cdtmbnma()] to do both in one call. Use
#' `summary()`, `coef()`, `plot()`, and `predict()` on the result. Estimation
#' uses Stan through `cmdstanr` or `rstan`.
#'
#' @section Interaction surfaces:
#' * `none`: additive component effects.
#' * `bilinear`: one parameter per pair; the interaction grows with the product
#'   of the normalised component doses.
#' * `gpdi`: a saturating surface adapted from the general pharmacodynamic
#'   interaction model; it is bilinear at low dose and plateaus as both doses
#'   increase.
#'
#' @keywords internal
"_PACKAGE"

.cdt_read_extdata <- function(filename) {
  f <- system.file("extdata", filename, package = "cdtmbnma")
  if (f == "") stop("Package data file not found: ", filename, call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

#' Read an installed package data file
#'
#' @param name Filename under `inst/extdata`.
#' @return A data frame.
#' @export
cdt_extdata <- function(name) {
  if (!is.character(name) || length(name) != 1L) stop("'name' must be a single filename.", call. = FALSE)
  .cdt_read_extdata(name)
}

#' Sacubitril and valsartan blood-pressure dose plane
#'
#' Reads the eight-week sitting systolic blood pressure dose plane used in the
#' package vignette. The data contain a valsartan dose-ranging axis, a fixed-ratio
#' diagonal, and off-diagonal arms that vary sacubitril at fixed valsartan.
#'
#' @return A data frame with columns `study`, `arm`, `d_sac`, `d_val`, `n`, `y`,
#'   `sd`, `se`, and `sd_source`.
#' @export
sacval_example <- function() {
  .cdt_read_extdata("sacval_msSBP.csv")
}

#' Amlodipine-anchored hypertension factorial extraction template
#'
#' Reads the package's antihypertensive factorial extraction form or its data
#' dictionary. The template is intended as a structured data-entry scaffold for a
#' calcium-channel blocker by angiotensin receptor blocker component-dose
#' network.
#'
#' @param type Either `"extraction"` for the extraction form or `"dictionary"`
#'   for the data dictionary.
#' @return A data frame.
#' @export
antihtn_factorial_template <- function(type = c("extraction", "dictionary")) {
  type <- match.arg(type)
  file <- if (type == "extraction") "antihtn_factorial_extraction_form.csv" else "antihtn_factorial_data_dictionary.csv"
  .cdt_read_extdata(file)
}

#' COPD BGF triple-therapy extraction template
#'
#' Reads the package's COPD budesonide/glycopyrronium/formoterol extraction form
#' or its data dictionary. The template provides a structured scaffold for a
#' three-component respiratory network with lung-function and exacerbation
#' outcomes.
#'
#' @param type Either `"extraction"` for the extraction form or `"dictionary"`
#'   for the data dictionary.
#' @return A data frame.
#' @export
copd_bgf_template <- function(type = c("extraction", "dictionary")) {
  type <- match.arg(type)
  file <- if (type == "extraction") "copd_bgf_extraction_form.csv" else "copd_bgf_data_dictionary.csv"
  .cdt_read_extdata(file)
}
