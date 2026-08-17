# These values were recorded from the official minimize.m Rosenbrock example.
# They are an external numerical oracle for the supported CG integration, not a
# contract for Rasmussen's private R implementation.
rasmussen_cg_oracle <- list(
  objective_values = c(
    1.00000000000000,
    0.77160942667725,
    0.58224024884105,
    0.40492742502160,
    0.32466327341368,
    0.28960411147824,
    0.07623420070067,
    0.06786211944378,
    0.03378423679313,
    0.00108990808914,
    0.00108795243321,
    0.00008974308332,
    0.00000012183819,
    0.00000000675602,
    0.00000000000000,
    0
  ),
  cumulative_evaluations = c(
    0,
    6,
    8,
    9,
    13,
    16,
    19,
    22,
    27,
    31,
    34,
    41,
    44,
    46,
    49,
    51
  ),
  gradient_squared_norms = c(
    2,
    4.853,
    7.288,
    2.231,
    3.540,
    6.203,
    3.739,
    0.454,
    3.098,
    0.0406,
    0.0454,
    0.394,
    0.00518,
    0.00308,
    1.52e-6,
    1.21e-9
  ),
  step_lengths = c(
    0,
    0.156,
    0.172,
    0.0892,
    0.0995,
    0.0860,
    0.383,
    0.00931,
    0.176,
    0.274,
    9.62e-5,
    0.0651,
    0.00697,
    6.366e-4,
    1.003e-4,
    1.167e-7
  ),
  parameters = c(1, 1)
)

run_rasmussen_cg_oracle <- function(objective) {
  optimizer <- make_opt(
    make_stages(
      gradient_stage(
        direction = cg_direction(cg_update = pr_update),
        step_size = rasmussen_ls(
          initial_step_length = "r",
          max_alpha_mult = 10
        )
      ),
      verbose = FALSE
    )
  )
  optimizer$count_res_fg <- FALSE

  opt_loop(
    optimizer,
    c(0, 0),
    objective,
    15,
    store_progress = TRUE,
    verbose = FALSE,
    grad_tol = 0,
    abs_tol = 0
  )
}

test_that("Rasmussen CG retains its minimize.m output oracle", {
  objective_without_combined <- rosenbrock_fg
  objective_without_combined$fg <- NULL
  objectives <- list(
    combined_callback = rosenbrock_fg,
    separate_callbacks = objective_without_combined
  )

  for (name in names(objectives)) {
    result <- run_rasmussen_cg_oracle(objectives[[name]])

    expect_equal(
      result$progress$nf,
      rasmussen_cg_oracle$cumulative_evaluations,
      info = name
    )
    expect_equal(
      result$progress$ng,
      rasmussen_cg_oracle$cumulative_evaluations,
      info = name
    )
    expect_equal(
      result$progress$f,
      rasmussen_cg_oracle$objective_values,
      tolerance = 1e-3,
      info = name
    )
    expect_equal(
      result$progress$g2n,
      rasmussen_cg_oracle$gradient_squared_norms,
      tolerance = 1e-3,
      info = name
    )
    expect_equal_abs(
      result$progress$step,
      rasmussen_cg_oracle$step_lengths,
      tolerance = 1e-3,
      info = name
    )
    expect_equal(
      result$par,
      rasmussen_cg_oracle$parameters,
      tolerance = 1e-3,
      info = name
    )
  }
})
