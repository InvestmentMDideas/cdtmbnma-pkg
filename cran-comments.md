## Resubmission

This is a resubmission of cdtmbnma (now 0.2.1). Thank you for the review. The
following changes address each point raised.

* **Software names in single quotes.** The `Description` field now writes
  'Stan' in single quotes, alongside the existing 'cmdstanr' and 'rstan'. The
  `Title` field contains no software, package, or API names.

* **References describing the methods.** The `Description` field now cites the
  methods the package implements, in the requested `authors (year) <doi:...>`
  form, with no space after `doi:` and with angle brackets for auto-linking:
  - Welton et al. (2009) <doi:10.1093/aje/kwp014> — additive component network
    meta-analysis.
  - Mawdsley et al. (2016) <doi:10.1002/psp4.12091> — model-based network
    meta-analysis with dose-response.
  - Wicha et al. (2017) <doi:10.1038/s41467-017-01929-y> — the general
    pharmacodynamic interaction surface used by `interaction = "gpdi"`.
  - Pedder et al. (2019) <doi:10.1002/jrsm.1351> — the time-course model used by
    the stage-one longitudinal interface.

* **`\dontrun{}` replaced by `\donttest{}`.** The two fitting examples
  (`cdt_fit()` and `cdtmbnma()`) no longer use `\dontrun{}`; no `\dontrun{}`
  remains anywhere in the package. They use `\donttest{}` rather than running
  unwrapped because they take longer than 5 seconds: `cdt_fit()` compiles a Stan
  program with the system C++ toolchain on first use. Under
  `R CMD check --as-cran` the `\donttest` examples run in their own step and
  pass in about 61 seconds, of which roughly 45-50 seconds is that one-time
  compilation. The two examples share the compiled model within a session, so
  only the first pays the cost.

* **Examples using packages in `Suggests`.** Both fitting examples are wrapped in
  `if (requireNamespace("rstan", quietly = TRUE)) { }`. They pass
  `backend = "rstan"` explicitly so the example never depends on a CmdStan
  installation being present. The design-building portion of the `cdt_fit()`
  example was moved out of `\donttest{}` and now runs unconditionally, since it
  needs no Stan backend.

Sampler settings in the examples were reduced (`chains = 2`, 500 warmup, 500
sampling) to keep run time down.

While making these changes we also fixed a related robustness problem:
`backend = "auto"` previously selected 'cmdstanr' whenever the cmdstanr
namespace could be loaded, even with no CmdStan installation behind it, which
then failed at model compilation. It now confirms CmdStan is available before
choosing that backend, and otherwise falls back to 'rstan'. The vignette applies
the same check.

## Test environments

* Local: Ubuntu 24.04, R 4.6.0 — `R CMD check --as-cran` (with PDF manual)
* R-hub v2 (GitHub Actions), R-devel on Linux, Windows and macOS — `Status: OK`
  on all three, with no errors, warnings or notes
* win-builder: R-devel and R-release — submitted

On every R-hub platform the `\donttest` examples were executed
(`checking examples with --run-donttest ... OK`, 71-76 seconds) and raised no
timing note.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  Maintainer: 'Tyler Pitre <pitretmed@gmail.com>'

  New submission

  This is the package's first submission to CRAN.

  The same NOTE reports 'cmdstanr' under "Suggests or Enhances not in mainstream
  repositories", and confirms it is available through the
  Additional_repositories field (see below).

  The incoming check may also flag possibly misspelled words in DESCRIPTION:
  "asymptote", "bilinear" and "timepoint" are correctly spelled statistical
  terms; "cmdstanr" and "rstan" are package names; "et", "al" and the surnames
  "Welton", "Mawdsley", "Wicha" and "Pedder" come from the added references.

A second NOTE appears only locally: "checking HTML version of manual ...
Skipping checking HTML validation: no command 'tidy' found." That is a missing
local tool, not a package issue.

## Suggested non-mainstream package

The package suggests 'cmdstanr', which is not on a mainstream repository. It is
available via the Additional_repositories field
(https://stan-dev.r-universe.dev). Model fitting can use either 'cmdstanr' (with
CmdStan) or 'rstan' (on CRAN). All examples, tests, and vignette code that need a
Stan backend are conditional on a usable backend being present: the test suite
uses synthetic posterior draws, the man-page fitting examples are wrapped in
`\donttest{}` and guarded with `requireNamespace("rstan")`, and the vignette's
sampling chunks are guarded with `eval = FALSE` when no backend is available. The
package therefore installs and checks with no Stan backend available.

## Downstream dependencies

There are currently no downstream dependencies (new package).
