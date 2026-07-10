# cdtmbnma

`cdtmbnma` is an R source package for Bayesian component, dose, and time
model-based network meta-analysis with dose-dependent interaction surfaces.

The core model fits arm-level component-dose networks. Each treatment arm is
represented by a dose vector, where zero means the component is absent. The model
estimates component Emax dose-response curves and, when the data identify them,
pairwise dose-dependent interaction surfaces. It supports continuous outcomes
with known arm standard errors and binary outcomes through a binomial-logit
likelihood.

The package also includes a stage-one longitudinal model for repeated arm means:
a two-component exponential time-course model with a bilinear interaction on the
long-term asymptote.

## Contents

- `cdt_data()` prepares single-timepoint arm-level data for Stan.
- `cdt_fit()` fits the single-timepoint model.
- `cdtmbnma()` is a one-call wrapper around `cdt_data()` and `cdt_fit()`.
- `cdt_time_data()` and `cdt_time_fit()` prepare and fit the two-component
  longitudinal model.
- `predict()`, `summary()`, `coef()`, and `plot()` methods are provided.
- `sacval_example()` reads the sacubitril/valsartan blood-pressure dose plane.
- `antihtn_factorial_template()` and `copd_bgf_template()` read structured
  extraction templates included with the package.

## Installation

A Stan backend is required to fit models. `cmdstanr` is recommended.

```r
install.packages("posterior")
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

install.packages("remotes")
remotes::install_local("cdtmbnma_0.2.1.tar.gz", build_vignettes = TRUE)
```

`rstan` can be used instead of `cmdstanr` by installing `rstan` and calling
`cdt_fit(..., backend = "rstan")` or `cdt_time_fit(..., backend = "rstan")`.

## Quick start: single-timepoint dose plane

```r
library(cdtmbnma)

sv <- sacval_example()

sv_design <- cdt_data(
  sv,
  study = "study",
  components = c("d_sac", "d_val"),
  outcome = "continuous",
  y = "y",
  se = "se",
  dstar = c(d_sac = 200, d_val = 320)
)

fit <- cdt_fit(
  sv_design,
  interaction = "bilinear",
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 1
)

summary(fit)
coef(fit)
plot(fit)
predict(fit, data.frame(d_sac = 100, d_val = 160))
```

For a binary outcome, use `outcome = "binary"` and supply `events` and
`n_binary` instead of `y` and `se`.

## Interaction surfaces

For component doses `d_c` and `d_c'`, the single-timepoint model has three
options.

The additive model uses no interaction:

```r
cdt_fit(design, interaction = "none")
```

The bilinear model adds one parameter per component pair:

\[
\psi_{cc'}(d_c, d_{c'}) = \eta_{cc'}
\left(\frac{d_c}{d_c^\star}\right)
\left(\frac{d_{c'}}{d_{c'}^\star}\right).
\]

The saturating general pharmacodynamic interaction surface is:

\[
\psi_{cc'}(d_c, d_{c'}) = \mathrm{INT}_{cc'}
\frac{d_c}{\kappa_c + d_c}
\frac{d_{c'}}{\kappa_{c'} + d_{c'}}.
\]

Both surfaces are centred on additivity. They vanish whenever either component is
absent.

## Quick start: two-component time-course model

```r
long_dat <- data.frame(
  study = rep("trial1", 6),
  arm = rep(c("placebo", "A", "AB"), each = 2),
  week = rep(c(4, 8), 3),
  dose_A = rep(c(0, 10, 10), each = 2),
  dose_B = rep(c(0, 0, 20), each = 2),
  y = c(0, 0, -1, -2, -2, -3),
  se = rep(1, 6)
)

time_design <- cdt_time_data(
  long_dat,
  study = "study",
  arm = "arm",
  time = "week",
  components = c("dose_A", "dose_B"),
  y = "y",
  se = "se",
  dstar = c(dose_A = 10, dose_B = 20)
)

# fit_time <- cdt_time_fit(time_design)
```

The time-course interface is intentionally narrower than the main interface: two
components, continuous outcomes, exponential time-course, and a bilinear
interaction on the long-term asymptote.

## Notes on interpretation

The null is additivity on the modelled scale: mean difference for continuous
outcomes and log odds for binary outcomes. A favourable interaction relative to
that additive component model is a statistical interaction on the chosen scale;
it is not automatically equivalent to Bliss independence, Loewe additivity, or a
pharmacological mechanism.

The interaction surface is identifiable only for component pairs that co-occur in
at least one arm. If no component pair co-occurs, `cdt_fit()` falls back to the
additive model.
