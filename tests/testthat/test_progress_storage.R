test_that("progress preserves late columns, types and missing observations", {
  # Synthetic summaries expose field arrival/disappearance independently of
  # callback caching and cross the buffer's capacity boundaries.
  progress <- make_progress()
  for (i in 0:299) {
    row <- list(iter = i * 2L, nf = i, step = 0.5)
    if (i > 0) row$ls_reason <- "wolfe"
    if (i == 2) row$f <- 3
    progress$add(row)
  }
  frame <- progress$finish()
  expected <- data.frame(
    nf = 0:299,
    step = rep(0.5, 300),
    ls_reason = c(NA_character_, rep("wolfe", 299)),
    f = rep(NA_real_, 300),
    row.names = as.character((0:299) * 2L)
  )
  expected$f[3] <- 3
  expect_identical(frame, expected)
})

test_that("public progress preserves cadence and the final iteration", {
  fg <- list(fn = function(x) sum(x), gr = function(x) rep(1, length(x)))
  for (iterations in c(0, 1, 9, 260)) {
    result <- mize(
      0,
      fg,
      method = "SD",
      line_search = "const",
      step0 = 1,
      max_iter = iterations,
      check_conv_every = 2,
      log_every = 4,
      abs_tol = 0,
      rel_tol = NULL,
      step_tol = NULL,
      store_progress = TRUE
    )
    labels <- unique(c(seq(0, iterations, by = 4), iterations))
    expect_identical(rownames(result$progress), as.character(labels))
    expect_equal(result$progress$f, -labels)
    expect_equal(tail(result$progress$step, 1), as.numeric(iterations > 0))
    expect_false(anyDuplicated(rownames(result$progress)) > 0)
  }
})
