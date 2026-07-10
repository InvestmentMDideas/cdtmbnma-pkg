test_that("cdt_data builds the design and detects co-occurring pairs", {
  df <- data.frame(
    study = c("s1", "s1", "s1", "s2", "s2"),
    A = c(0, 1, 0, 0, 2),
    B = c(0, 0, 1, 1, 1),
    y = c(0, -2, -1, -1.5, -3),
    se = c(1, 0.5, 0.5, 0.7, 0.6)
  )
  d <- suppressWarnings(
    cdt_data(df, study = "study", components = c("A", "B"),
             outcome = "continuous", y = "y", se = "se"))
  expect_s3_class(d, "cdt_data")
  expect_equal(d$standata$N, 5L)
  expect_equal(d$standata$S, 2L)
  expect_equal(d$standata$C, 2L)
  # A and B co-occur only in s2's last arm -> exactly one interaction pair
  expect_equal(nrow(d$pairs), 1L)
  # baseline: s1 placebo (row 1); s2 has no placebo -> its lowest-dose arm (row 4)
  expect_equal(d$standata$is_ref, c(1L, 0L, 0L, 1L, 0L))
  # dstar defaults to the max observed dose per component
  expect_equal(unname(d$dstar), c(2, 1))
})

test_that("se is derived from sd and n when se is absent", {
  df <- data.frame(study = c("s1", "s1"), A = c(0, 1), B = c(0, 0),
                   y = c(0, -2), sd = c(4, 4), n = c(16, 16))
  d <- cdt_data(df, study = "study", components = c("A", "B"),
                outcome = "continuous", y = "y", sd = "sd", n = "n")
  expect_equal(d$standata$se, c(1, 1))
})

test_that("binary outcome is encoded", {
  df <- data.frame(study = c("s1", "s1"), A = c(0, 1), B = c(0, 0),
                   r = c(5, 3), n = c(20, 20))
  d <- cdt_data(df, study = "study", components = c("A", "B"),
                outcome = "binary", events = "r", n_binary = "n",
                dstar = c(A = 1, B = 1))
  expect_equal(d$standata$outcome, 2L)
  expect_equal(d$standata$r, c(5L, 3L))
  expect_equal(unname(d$dstar), c(1, 1))
})

test_that("malformed input is rejected", {
  df <- data.frame(study = "s1", A = -1, B = 0, y = 0, se = 1)
  expect_error(cdt_data(df, study = "study", components = c("A", "B"),
                        outcome = "continuous", y = "y", se = "se"),
               "non-negative")
})
