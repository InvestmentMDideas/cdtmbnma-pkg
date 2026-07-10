test_that("installed data readers return data frames", {
  sv <- sacval_example()
  expect_s3_class(sv, "data.frame")
  expect_true(all(c("study", "d_sac", "d_val", "y", "se") %in% names(sv)))
  expect_s3_class(antihtn_factorial_template("dictionary"), "data.frame")
  expect_s3_class(copd_bgf_template("extraction"), "data.frame")
})
