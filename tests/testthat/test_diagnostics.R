capture_mize_messages <- function(expr) {
  messages <- character()
  withCallingHandlers(
    force(expr),
    message = function(cnd) {
      messages <<- c(messages, trimws(conditionMessage(cnd)))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

test_that("opt_report prints optional progress fields", {
  step_info <- list(
    iter = 3,
    f = 1.25,
    g2n = 0.5,
    ginfn = 0.25,
    nf = 4,
    ng = 5,
    step = 0.125,
    alpha = 0.75
  )

  messages <- capture_mize_messages(
    opt_report(step_info, print_par = TRUE, par = c(1, -2))
  )

  expect_equal(
    messages,
    paste(
      "iter 3 f = 1.25 g2 = 0.5 ginf = 0.25 nf = 4 ng = 5",
      "step = 0.125 alpha = 0.75 par = 1, -2"
    )
  )
})

test_that("list_hooks reports registered handlers and dependency hooks", {
  quadratic_fg <- list(
    fn = function(x) {
      sum(x^2)
    },
    gr = function(x) {
      2 * x
    }
  )
  opt <- make_mize(
    method = "SD",
    line_search = "constant",
    step0 = 0.1,
    par = c(1, -2),
    fg = quadratic_fg,
    abs_tol = NULL,
    rel_tol = NULL,
    grad_tol = NULL
  )

  messages <- capture_mize_messages(list_hooks(opt))

  expect_equal(messages[1], "handlers")
  expect_true("hooks" %in% messages)
  expect_true("init opt" %in% messages)
  expect_true("during gradient_descent direction" %in% messages)
  expect_true("during gradient_descent step_size" %in% messages)
  expect_true("before gradient_descent: gradient" %in% messages)
})
