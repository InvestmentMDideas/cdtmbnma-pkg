test_that("priors validate", {
  expect_named(cdt_priors(), c("prior_emax_sd", "prior_logED50_mean", "prior_logED50_sd",
                               "prior_int_sd", "prior_ref_sd", "prior_omega_sd"))
  expect_error(cdt_priors(emax_sd = 0), "positive")
  expect_named(cdt_time_priors(), c("prior_ref_E_sd", "prior_ref_lograte_mean", "prior_ref_lograte_sd", "prior_emax_sd", "prior_logED50_mean", "prior_logED50_sd", "prior_eta_sd", "prior_rate_sd", "prior_omega_E_sd", "prior_omega_k_sd"))
})

test_that("predict works with posterior draws without calling Stan", {
  df <- data.frame(
    study = c("s1", "s1", "s1"),
    A = c(0, 1, 1),
    B = c(0, 0, 1),
    y = c(0, -1, -2),
    se = c(1, 1, 1)
  )
  d <- cdt_data(df, study = "study", components = c("A", "B"),
                outcome = "continuous", y = "y", se = "se",
                dstar = c(A = 1, B = 1))
  draws_df <- data.frame(
    "emax[1]" = c(-2, -3, -4),
    "emax[2]" = c(-1, -1, -1),
    "ED50[1]" = c(1, 1, 1),
    "ED50[2]" = c(1, 1, 1),
    "eta[1]" = c(0, -0.5, -1),
    omega = c(0.1, 0.2, 0.3),
    check.names = FALSE
  )
  fit <- structure(list(draws = posterior::as_draws_df(draws_df), data = d,
                        backend = "test", spec = list(interaction = "bilinear")),
                   class = "cdtmbnma")
  p <- predict(fit, data.frame(A = 1, B = 1))
  expect_equal(nrow(p), 1L)
  expect_true("mean" %in% names(p))
  expect_true(p$mean < 0)
})

test_that("stage-one time-course prediction works with posterior draws", {
  dat <- data.frame(
    study = rep(c("s1", "s1", "s1"), each = 2),
    arm = rep(c("pbo", "A", "AB"), each = 2),
    week = rep(c(4, 8), 3),
    A = rep(c(0, 10, 10), each = 2),
    B = rep(c(0, 0, 20), each = 2),
    y = c(0, 0, -1, -2, -2, -3),
    se = rep(1, 6)
  )
  d <- cdt_time_data(dat, study = "study", arm = "arm", time = "week",
                     components = c("A", "B"), y = "y", se = "se",
                     dstar = c(A = 10, B = 20))
  expect_s3_class(d, "cdt_time_data")
  draws_df <- data.frame(a1_EmaxA = c(-3, -4), a2_EmaxB = c(-1, -2),
                         ED50A = c(10, 10), ED50B = c(20, 20), eta = c(0, -1),
                         check.names = FALSE)
  fit <- structure(list(draws = posterior::as_draws_df(draws_df), data = d,
                        backend = "test"), class = "cdt_timefit")
  p <- predict(fit, data.frame(A = 10, B = 20))
  expect_equal(nrow(p), 1L)
  expect_true(p$mean < 0)
})
