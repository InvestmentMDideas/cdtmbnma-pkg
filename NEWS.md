# cdtmbnma 0.2.1

- Quotes software names and adds method references with DOIs in the `DESCRIPTION`
  `Description` field, per CRAN review.
- Replaces `\dontrun{}` with `\donttest{}` in the `cdt_fit()` and `cdtmbnma()`
  examples, guarded by `requireNamespace("rstan")`, and runs the design-building
  part of the `cdt_fit()` example unconditionally.
- `backend = "auto"` no longer selects `cmdstanr` when its namespace loads but
  no CmdStan installation is available; it now falls back to `rstan`. The
  vignette applies the same check.

# cdtmbnma 0.2.0

- Adds complete manual pages and package metadata for a full source release.
- Adds `cdt_time_data()` and `cdt_time_fit()` for the two-component stage-one
  longitudinal model in `inst/stan/cdt_stage1.stan`.
- Adds `cdt_time_priors()` and S3 `print()`, `summary()`, `coef()`, and
  `predict()` methods for time-course fits.
- Adds structured example-data readers for the sacubitril/valsartan dose plane,
  antihypertensive factorial extraction template, and COPD BGF extraction
  template.
- Hardens input validation for component doses, references, interaction pairs,
  priors, and prediction data.
- Updates Stan generated-quantities prediction grids to avoid zero-length
  prediction data.

# cdtmbnma 0.1.0

- Initial source release with single-timepoint component dose-response network
  meta-analysis, continuous and binary likelihoods, bilinear and saturating
  interaction surfaces, Stan models, tests, and a sacubitril/valsartan example.
