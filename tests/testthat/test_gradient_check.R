make_gradient_check_quadratic <- function() {
  a <- matrix(
    c(
      4, 1, 0,
      1, 3, 0.5,
      0, 0.5, 2
    ),
    nrow = 3
  )
  b <- c(1, -2, 0.5)

  fn <- function(x) {
    drop(0.5 * t(x) %*% a %*% x - t(b) %*% x)
  }
  gr <- function(x) {
    as.vector(a %*% x - b)
  }

  list(
    fn = fn,
    gr = gr,
    fg = function(x) {
      list(fn = fn(x), gr = gr(x))
    }
  )
}

test_that("check_mize_gradient reports small errors for correct gradients", {
  fg <- make_gradient_check_quadratic()
  par <- c(alpha = 0.25, beta = -1.5, gamma = 2)

  res <- check_mize_gradient(fg, par)

  expect_equal(res$method, "central")
  expect_equal(res$par, par)
  expect_equal(res$by_coord$index, seq_along(par))
  expect_equal(res$by_coord$name, names(par))
  expect_lt(res$max_abs_error, 1e-8)
  expect_lt(res$max_rel_error, 1e-8)
  expect_equal(res$worst_abs$abs_error, max(res$abs_error))
  expect_equal(res$worst_rel$rel_error, max(res$rel_error))

  expect_true(res$fg_consistency$checked)
  expect_equal(res$fg_consistency$fn$abs_error, 0)
  expect_equal(res$fg_consistency$gr$max_abs_error, 0)
})

test_that("check_mize_gradient identifies deliberately bad gradients", {
  fg <- make_gradient_check_quadratic()
  bad_fg <- fg
  bad_fg$gr <- function(x) {
    fg$gr(x) + c(0, 0.25, 0)
  }

  res <- check_mize_gradient(bad_fg, c(0.25, -1.5, 2))

  expect_gt(res$max_abs_error, 0.24)
  expect_equal(res$worst_abs$index, 2)
  expect_equal(res$worst_abs$gr - res$worst_abs$gr_fd, res$worst_abs$error)
})

test_that("check_mize_gradient supports forward differences", {
  fg <- make_gradient_check_quadratic()
  par <- c(0.25, -1.5, 2)

  res <- check_mize_gradient(fg, par, method = "forward", rel_eps = 1e-6)

  expect_equal(res$method, "forward")
  expect_lt(res$max_abs_error, 1e-4)
})

test_that("check_mize_gradient reports combined fg inconsistencies", {
  fg <- make_gradient_check_quadratic()
  inconsistent_fg <- fg
  inconsistent_fg$fg <- function(x) {
    combined <- fg$fg(x)
    combined$fn <- combined$fn + 1
    combined$gr[3] <- combined$gr[3] - 0.5
    combined
  }

  res <- check_mize_gradient(inconsistent_fg, c(0.25, -1.5, 2))

  expect_lt(res$max_abs_error, 1e-8)
  expect_true(res$fg_consistency$checked)
  expect_equal(res$fg_consistency$fn$abs_error, 1)
  expect_gt(res$fg_consistency$gr$max_abs_error, 0.49)
  expect_equal(res$fg_consistency$gr$worst_abs$index, 3)
})

test_that("check_mize_gradient validates public inputs", {
  fg <- make_gradient_check_quadratic()

  expect_error(
    check_mize_gradient(list(fn = fg$fn), c(1)),
    "fg\\$gr must be a function"
  )
  expect_error(
    check_mize_gradient(fg, matrix(1, nrow = 1)),
    "par must be a non-empty numeric vector"
  )
  expect_error(
    check_mize_gradient(fg, c(1), rel_eps = 0),
    "rel_eps must be a positive finite numeric scalar"
  )
})
