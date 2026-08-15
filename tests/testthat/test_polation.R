test_that("Cubic interpolation can return NaN", {
  xc <- cubic_interpolate(
    0.9493646,
    1.04033,
    0.6949782,
    0.1898729,
    0.8138197,
    0.03981831,
    ignoreWarnings = TRUE
  )
  testthat::expect_true(is.nan(xc))
})

test_that("Cubic extrapolation warning suppression preserves its result", {
  default_result <- cubic_extrapolate(0, 1, -2, 0.5, 0.25, -1)
  ignored_result <- cubic_extrapolate(
    0,
    1,
    -2,
    0.5,
    0.25,
    -1,
    ignoreWarnings = TRUE
  )

  expect_equal(default_result, 1)
  expect_identical(default_result, ignored_result)

  expect_warning(
    default_invalid <- cubic_extrapolate(
      0.9493646,
      1.04033,
      0.6949782,
      0.1898729,
      0.8138197,
      0.03981831
    ),
    "NaNs produced"
  )
  expect_warning(
    ignored_invalid <- cubic_extrapolate(
      0.9493646,
      1.04033,
      0.6949782,
      0.1898729,
      0.8138197,
      0.03981831,
      ignoreWarnings = TRUE
    ),
    NA
  )
  expect_true(is.nan(default_invalid))
  expect_identical(default_invalid, ignored_invalid)
})
