#!/usr/bin/env Rscript

source(
  file.path("scripts", "funconstrain-mize-harness.R"),
  local = environment()
)

hessian_probe_control <- function() {
  epsilon <- .Machine$double.eps
  list(
    step_factor = epsilon^(1 / 3),
    symmetry_relative_tolerance = sqrt(epsilon),
    directional_absolute_tolerance = 100 * sqrt(epsilon),
    directional_relative_tolerance = 100 * sqrt(epsilon),
    directional_roundoff_multiplier = 1
  )
}

hessian_probe_norm <- function(x) {
  if (length(x) == 0L) {
    return(0)
  }

  magnitude <- abs(x)
  scale <- max(magnitude)
  if (!is.finite(scale) || scale == 0) {
    return(scale)
  }

  scale * sqrt(sum((magnitude / scale)^2))
}

hessian_probe_directions <- function(point) {
  dimension <- length(point)
  seeds <- list(
    uniform = rep(1, dimension),
    alternating = rep(c(1, -1), length.out = dimension),
    trend = seq(-1, 1, length.out = dimension)
  )
  coordinate_scale <- pmax(1, abs(point))
  directions <- list()

  for (name in names(seeds)) {
    seed <- seeds[[name]]
    seed_norm <- hessian_probe_norm(seed)
    if (!is.finite(seed_norm) || seed_norm == 0) {
      next
    }

    unit_seed <- seed / seed_norm
    scaled_seed <- coordinate_scale * unit_seed
    direction_scale <- hessian_probe_norm(scaled_seed)
    if (!is.finite(direction_scale) || direction_scale == 0) {
      stop(
        "could not construct a finite nonzero Hessian probe direction",
        call. = FALSE
      )
    }
    direction <- scaled_seed / direction_scale
    duplicate <- any(vapply(
      directions,
      function(existing) {
        abs(sum(existing$direction * direction)) >= 1 - 1e-12
      },
      logical(1)
    ))
    if (duplicate) {
      next
    }

    directions[[name]] <- list(
      direction = direction,
      direction_scale = direction_scale
    )
  }

  directions
}

hessian_probe_callback <- function(callback, point) {
  tryCatch(
    list(ok = TRUE, value = callback(point), error = ""),
    error = function(err) {
      list(ok = FALSE, value = NULL, error = conditionMessage(err))
    }
  )
}

hessian_spectrum_control <- function() {
  list(
    absolute_tolerance = 0,
    relative_tolerance = sqrt(.Machine$double.eps)
  )
}

hessian_spectrum_symmetrize <- function(hessian) {
  transposed <- t(hessian)
  if (all(hessian == transposed)) {
    return(list(
      value = hessian,
      method = "preserved_exact_symmetry"
    ))
  }

  left <- as.vector(hessian)
  right <- as.vector(transposed)
  overflow_risk <- sign(left) == sign(right) &
    abs(left) > .Machine$double.xmax - abs(right)
  averaged <- numeric(length(left))
  safe_sum <- !overflow_risk
  averaged[safe_sum] <- (left[safe_sum] + right[safe_sum]) / 2
  averaged[overflow_risk] <-
    left[overflow_risk] / 2 + right[overflow_risk] / 2
  symmetric <- hessian
  symmetric[] <- averaged

  list(
    value = symmetric,
    method = "overflow_safe_pairwise_average"
  )
}

classify_hessian_spectrum <- function(
  hessian,
  expected_dimension,
  control = hessian_spectrum_control(),
  eigen_solver = base::eigen
) {
  required_controls <- c("absolute_tolerance", "relative_tolerance")
  if (
    !is.list(control) ||
      !all(required_controls %in% names(control)) ||
      any(!vapply(control[required_controls], is.numeric, logical(1))) ||
      any(vapply(control[required_controls], length, integer(1)) != 1L) ||
      any(!is.finite(unlist(control[required_controls]))) ||
      control$absolute_tolerance < 0 ||
      control$relative_tolerance <= 0
  ) {
    stop(
      "spectrum control tolerances must be finite numeric scalars, with absolute tolerance nonnegative and relative tolerance positive",
      call. = FALSE
    )
  }
  if (
    !is.numeric(expected_dimension) ||
      length(expected_dimension) != 1L ||
      !is.finite(expected_dimension) ||
      expected_dimension < 1 ||
      expected_dimension != as.integer(expected_dimension)
  ) {
    stop("expected_dimension must be a positive integer", call. = FALSE)
  }
  if (!is.function(eigen_solver)) {
    stop("eigen_solver must be a function", call. = FALSE)
  }

  expected_dimension <- as.integer(expected_dimension)
  hessian_shape_ok <- is.matrix(hessian) &&
    is.numeric(hessian) &&
    identical(dim(hessian), c(expected_dimension, expected_dimension))
  hessian_finite <- hessian_shape_ok && all(is.finite(hessian))
  result <- list(
    hessian_shape_ok = hessian_shape_ok,
    hessian_finite = hessian_finite,
    symmetrization = NA_character_,
    symmetrized_hessian_finite = FALSE,
    eigensolver_ok = FALSE,
    eigensolver_error = "",
    eigenvalues_finite = FALSE,
    eigenvalues = "",
    eigen_min = NA_real_,
    eigen_max = NA_real_,
    eigen_abs_min = NA_real_,
    eigen_abs_max = NA_real_,
    spectral_scale = NA_real_,
    spectral_absolute_tolerance = control$absolute_tolerance,
    spectral_relative_tolerance = control$relative_tolerance,
    spectral_absolute_tolerance_term = control$absolute_tolerance,
    spectral_relative_tolerance_term = NA_real_,
    spectral_sign_tolerance = NA_real_,
    inertia_positive = NA_integer_,
    inertia_zero = NA_integer_,
    inertia_negative = NA_integer_,
    inertia_unresolved = NA_integer_,
    has_resolved_positive_curvature = NA,
    has_resolved_negative_curvature = NA,
    singular = NA,
    singularity = "calculation_failed",
    spectral_condition_estimate = NA_real_,
    spectral_condition_definition = paste(
      "max(abs(lambda)) / min(abs(lambda));",
      "Inf when any exact zero; otherwise NA for unresolved or failed"
    ),
    spectral_classification = "calculation_failed"
  )

  if (!hessian_finite) {
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  symmetrized <- hessian_spectrum_symmetrize(hessian)
  symmetric_hessian <- symmetrized$value
  result$symmetrization <- symmetrized$method
  result$symmetrized_hessian_finite <- all(is.finite(symmetric_hessian))
  if (!result$symmetrized_hessian_finite) {
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  eigen_result <- tryCatch(
    eigen_solver(symmetric_hessian, symmetric = TRUE, only.values = TRUE),
    error = identity
  )
  if (inherits(eigen_result, "error")) {
    result$eigensolver_error <- conditionMessage(eigen_result)
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  result$eigensolver_ok <- TRUE

  eigenvalues <- eigen_result$values
  valid_eigenvalues <- is.numeric(eigenvalues) &&
    is.null(dim(eigenvalues)) &&
    length(eigenvalues) == expected_dimension
  if (!valid_eigenvalues) {
    result$eigensolver_error <- paste0(
      "eigensolver must return ",
      expected_dimension,
      " numeric values"
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  result$eigenvalues_finite <- all(is.finite(eigenvalues))
  if (!result$eigenvalues_finite) {
    result$eigensolver_error <- "eigensolver returned nonfinite values"
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  eigenvalues <- as.numeric(eigenvalues)
  absolute_eigenvalues <- abs(eigenvalues)
  spectral_scale <- max(absolute_eigenvalues)
  relative_tolerance_term <- control$relative_tolerance * spectral_scale
  sign_tolerance <- control$absolute_tolerance + relative_tolerance_term
  tolerance_values_finite <- all(is.finite(c(
    spectral_scale,
    relative_tolerance_term,
    sign_tolerance
  )))
  if (!tolerance_values_finite) {
    result$eigensolver_error <- "spectral scale or tolerance is nonfinite"
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  positive <- eigenvalues > sign_tolerance
  negative <- eigenvalues < -sign_tolerance
  exact_zero <- eigenvalues == 0
  unresolved <- !(positive | negative | exact_zero)
  inertia_positive <- sum(positive)
  inertia_zero <- sum(exact_zero)
  inertia_negative <- sum(negative)
  inertia_unresolved <- sum(unresolved)

  result$eigenvalues <- paste(
    format(eigenvalues, digits = 17, scientific = TRUE, trim = TRUE),
    collapse = ";"
  )
  result$eigen_min <- min(eigenvalues)
  result$eigen_max <- max(eigenvalues)
  result$eigen_abs_min <- min(absolute_eigenvalues)
  result$eigen_abs_max <- max(absolute_eigenvalues)
  result$spectral_scale <- spectral_scale
  result$spectral_relative_tolerance_term <- relative_tolerance_term
  result$spectral_sign_tolerance <- sign_tolerance
  result$inertia_positive <- inertia_positive
  result$inertia_zero <- inertia_zero
  result$inertia_negative <- inertia_negative
  result$inertia_unresolved <- inertia_unresolved
  result$has_resolved_positive_curvature <- inertia_positive > 0L
  result$has_resolved_negative_curvature <- inertia_negative > 0L

  if (inertia_zero > 0L) {
    result$singular <- TRUE
    result$singularity <- "exactly_singular"
    result$spectral_condition_estimate <- Inf
  } else if (inertia_unresolved > 0L) {
    result$singularity <- "numerically_unresolved"
  } else {
    result$singular <- FALSE
    result$singularity <- "resolved_nonsingular"
    result$spectral_condition_estimate <-
      result$eigen_abs_max / result$eigen_abs_min
  }

  if (inertia_unresolved > 0L) {
    result$spectral_classification <- "numerically_unresolved"
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  result$spectral_classification <- if (
    inertia_positive == expected_dimension
  ) {
    "positive_definite"
  } else if (inertia_negative == expected_dimension) {
    "negative_definite"
  } else if (inertia_zero == expected_dimension) {
    "zero"
  } else if (inertia_positive > 0L && inertia_negative > 0L) {
    "indefinite"
  } else if (inertia_positive > 0L && inertia_zero > 0L) {
    "positive_semidefinite"
  } else if (inertia_negative > 0L && inertia_zero > 0L) {
    "negative_semidefinite"
  } else {
    stop("unreachable Hessian spectrum classification", call. = FALSE)
  }

  as.data.frame(result, stringsAsFactors = FALSE)
}

probe_hessian_spectrum <- function(
  fg,
  point,
  control = hessian_spectrum_control(),
  eigen_solver = base::eigen
) {
  if (!is.list(fg) || !is.function(fg$hs)) {
    stop("fg must provide function callback fg$hs", call. = FALSE)
  }
  if (
    !is.numeric(point) ||
      !is.null(dim(point)) ||
      length(point) == 0L ||
      !all(is.finite(point))
  ) {
    stop("point must be a non-empty finite numeric vector", call. = FALSE)
  }

  hessian_result <- hessian_probe_callback(fg$hs, point)
  spectrum <- classify_hessian_spectrum(
    hessian = hessian_result$value,
    expected_dimension = length(point),
    control = control,
    eigen_solver = eigen_solver
  )
  data.frame(
    hessian_evaluation = "reevaluated_after_integrity_gate",
    spectrum_hs_calls = 1L,
    benchmark_callback_counted = FALSE,
    hessian_callback_ok = hessian_result$ok,
    hessian_callback_error = hessian_result$error,
    spectrum,
    stringsAsFactors = FALSE
  )
}

hessian_raw_newton_control <- function() {
  list(
    residual_absolute_tolerance = 0,
    residual_relative_tolerance = 100 * sqrt(.Machine$double.eps)
  )
}

hessian_raw_newton_sign <- function(value) {
  if (length(value) != 1L || is.na(value)) {
    return("unavailable")
  }
  if (!is.finite(value)) {
    return("nonfinite")
  }
  if (value < 0) {
    return("negative")
  }
  if (value > 0) {
    return("positive")
  }
  "zero"
}

hessian_raw_newton_eligibility <- function(spectrum) {
  required <- c(
    "spectral_classification",
    "singularity",
    "spectral_condition_estimate",
    "spectral_relative_tolerance",
    "inertia_zero",
    "inertia_unresolved"
  )
  if (
    !is.data.frame(spectrum) ||
      nrow(spectrum) != 1L ||
      !all(required %in% names(spectrum))
  ) {
    stop(
      "spectrum must be one Hessian spectrum row with solve-eligibility fields",
      call. = FALSE
    )
  }

  relative_tolerance <- spectrum$spectral_relative_tolerance[[1L]]
  condition_limit <- if (
    is.numeric(relative_tolerance) &&
      length(relative_tolerance) == 1L &&
      is.finite(relative_tolerance) &&
      relative_tolerance > 0
  ) {
    1 / relative_tolerance
  } else {
    NA_real_
  }
  classification <- spectrum$spectral_classification[[1L]]
  singularity <- spectrum$singularity[[1L]]
  condition_estimate <- spectrum$spectral_condition_estimate[[1L]]
  inertia_zero <- spectrum$inertia_zero[[1L]]
  inertia_unresolved <- spectrum$inertia_unresolved[[1L]]
  resolved_nonsingular_classes <- c(
    "positive_definite",
    "negative_definite",
    "indefinite"
  )

  basis <- if (identical(classification, "calculation_failed")) {
    "spectrum_calculation_failed"
  } else if (
    identical(singularity, "exactly_singular") ||
      (!is.na(inertia_zero) && inertia_zero > 0L)
  ) {
    "exactly_singular_spectrum"
  } else if (
    identical(classification, "numerically_unresolved") ||
      identical(singularity, "numerically_unresolved") ||
      (!is.na(inertia_unresolved) && inertia_unresolved > 0L)
  ) {
    "numerically_unresolved_spectrum"
  } else if (!identical(singularity, "resolved_nonsingular")) {
    "spectrum_not_resolved_nonsingular"
  } else if (!classification %in% resolved_nonsingular_classes) {
    "spectrum_class_not_resolved_nonsingular"
  } else if (!is.finite(condition_limit)) {
    "condition_limit_unavailable"
  } else if (
    !is.numeric(condition_estimate) ||
      length(condition_estimate) != 1L ||
      !is.finite(condition_estimate) ||
      condition_estimate <= 0
  ) {
    "condition_estimate_unavailable"
  } else if (condition_estimate > condition_limit) {
    "condition_limit_exceeded"
  } else {
    "resolved_nonsingular_spectrum_within_condition_limit"
  }

  list(
    eligible = identical(
      basis,
      "resolved_nonsingular_spectrum_within_condition_limit"
    ),
    basis = basis,
    condition_limit = condition_limit
  )
}

probe_raw_newton_direction <- function(
  fg,
  point,
  spectrum,
  control = hessian_raw_newton_control(),
  linear_solver = base::solve,
  linear_solver_name = "base_solve"
) {
  if (!is.list(fg) || !is.function(fg$gr) || !is.function(fg$hs)) {
    stop("fg must provide function callbacks fg$gr and fg$hs", call. = FALSE)
  }
  if (
    !is.numeric(point) ||
      !is.null(dim(point)) ||
      length(point) == 0L ||
      !all(is.finite(point))
  ) {
    stop("point must be a non-empty finite numeric vector", call. = FALSE)
  }

  required_controls <- c(
    "residual_absolute_tolerance",
    "residual_relative_tolerance"
  )
  if (
    !is.list(control) ||
      !all(required_controls %in% names(control)) ||
      any(!vapply(control[required_controls], is.numeric, logical(1))) ||
      any(vapply(control[required_controls], length, integer(1)) != 1L) ||
      any(!is.finite(unlist(control[required_controls]))) ||
      control$residual_absolute_tolerance < 0 ||
      control$residual_relative_tolerance <= 0
  ) {
    stop(
      paste(
        "raw-Newton residual tolerances must be finite numeric scalars,",
        "with absolute tolerance nonnegative and relative tolerance positive"
      ),
      call. = FALSE
    )
  }
  if (!is.function(linear_solver)) {
    stop("linear_solver must be a function", call. = FALSE)
  }
  if (
    !is.character(linear_solver_name) ||
      length(linear_solver_name) != 1L ||
      is.na(linear_solver_name) ||
      !nzchar(linear_solver_name)
  ) {
    stop("linear_solver_name must be a non-empty string", call. = FALSE)
  }

  eligibility <- hessian_raw_newton_eligibility(spectrum)
  result <- list(
    direction_evaluation = if (eligibility$eligible) {
      "reevaluated_after_integrity_and_spectrum_gates"
    } else {
      "not_reevaluated_spectrum_ineligible"
    },
    direction_gr_calls = 0L,
    direction_hs_calls = 0L,
    benchmark_callback_counted = FALSE,
    stable_solve_eligible = eligibility$eligible,
    stable_solve_basis = eligibility$basis,
    stable_condition_limit = eligibility$condition_limit,
    stable_condition_limit_definition = "1 / spectral_relative_tolerance",
    linear_solver = linear_solver_name,
    solve_attempted = FALSE,
    solve_success = FALSE,
    solve_status = "spectrum_ineligible",
    solve_error = "",
    gradient_callback_ok = NA,
    gradient_callback_error = "",
    gradient_shape_ok = NA,
    gradient_finite = NA,
    gradient_norm = NA_real_,
    hessian_callback_ok = NA,
    hessian_callback_error = "",
    hessian_shape_ok = NA,
    hessian_finite = NA,
    symmetrization = NA_character_,
    symmetrized_hessian_finite = NA,
    raw_direction_values = "",
    raw_direction_shape_ok = NA,
    raw_direction_finite = NA,
    raw_direction_norm = NA_real_,
    system_product_norm = NA_real_,
    solve_residual_norm = NA_real_,
    solve_residual_scale = NA_real_,
    solve_residual_absolute_tolerance = control$residual_absolute_tolerance,
    solve_residual_relative_tolerance = control$residual_relative_tolerance,
    solve_residual_absolute_tolerance_term = control$residual_absolute_tolerance,
    solve_residual_relative_tolerance_term = NA_real_,
    solve_residual_threshold = NA_real_,
    solve_residual_values_finite = NA,
    solve_residual_pass = NA,
    directional_slope = NA_real_,
    directional_slope_sign = "unavailable",
    raw_direction_is_descent = NA,
    quadratic_curvature = NA_real_,
    quadratic_curvature_sign = "unavailable",
    predicted_decrease = NA_real_,
    predicted_decrease_sign = "unavailable",
    predicted_decrease_positive = NA,
    predicted_decrease_definition = "-(g' p + 0.5 p' H p)",
    direction_diagnostics_complete = FALSE
  )

  if (!eligibility$eligible) {
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  gradient_result <- hessian_probe_callback(fg$gr, point)
  hessian_result <- hessian_probe_callback(fg$hs, point)
  result$direction_gr_calls <- 1L
  result$direction_hs_calls <- 1L
  result$gradient_callback_ok <- gradient_result$ok
  result$gradient_callback_error <- gradient_result$error
  result$hessian_callback_ok <- hessian_result$ok
  result$hessian_callback_error <- hessian_result$error

  dimension <- length(point)
  gradient <- gradient_result$value
  result$gradient_shape_ok <- gradient_result$ok &&
    is.numeric(gradient) &&
    is.null(dim(gradient)) &&
    length(gradient) == dimension
  result$gradient_finite <- result$gradient_shape_ok &&
    all(is.finite(gradient))
  hessian <- hessian_result$value
  result$hessian_shape_ok <- hessian_result$ok &&
    is.matrix(hessian) &&
    is.numeric(hessian) &&
    identical(dim(hessian), c(dimension, dimension))
  result$hessian_finite <- result$hessian_shape_ok &&
    all(is.finite(hessian))

  if (!gradient_result$ok) {
    result$solve_status <- "gradient_callback_failed"
    result$solve_error <- gradient_result$error
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  if (!result$gradient_shape_ok || !result$gradient_finite) {
    result$solve_status <- "gradient_value_invalid"
    result$solve_error <- paste0(
      "gradient must be a finite numeric vector of length ",
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  result$gradient_norm <- hessian_probe_norm(gradient)
  if (!hessian_result$ok) {
    result$solve_status <- "hessian_callback_failed"
    result$solve_error <- hessian_result$error
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  if (!result$hessian_shape_ok || !result$hessian_finite) {
    result$solve_status <- "hessian_value_invalid"
    result$solve_error <- paste0(
      "Hessian must be a finite ",
      dimension,
      " x ",
      dimension,
      " numeric matrix"
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  symmetrized <- tryCatch(
    hessian_spectrum_symmetrize(hessian),
    error = identity
  )
  if (inherits(symmetrized, "error")) {
    result$solve_status <- "hessian_symmetrization_failed"
    result$solve_error <- conditionMessage(symmetrized)
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  system_hessian <- symmetrized$value
  result$symmetrization <- symmetrized$method
  result$symmetrized_hessian_finite <- all(is.finite(system_hessian))
  if (!result$symmetrized_hessian_finite) {
    result$solve_status <- "hessian_symmetrization_failed"
    result$solve_error <- "symmetrized Hessian is nonfinite"
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  result$solve_attempted <- TRUE
  solve_result <- tryCatch(
    linear_solver(system_hessian, -gradient),
    error = identity
  )
  if (inherits(solve_result, "error")) {
    result$solve_status <- "solver_error"
    result$solve_error <- conditionMessage(solve_result)
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  direction <- solve_result
  result$raw_direction_shape_ok <- is.numeric(direction) &&
    is.null(dim(direction)) &&
    length(direction) == dimension
  result$raw_direction_finite <- result$raw_direction_shape_ok &&
    all(is.finite(direction))
  if (!result$raw_direction_shape_ok || !result$raw_direction_finite) {
    result$solve_status <- "solver_return_invalid"
    result$solve_error <- paste0(
      "linear solver must return a finite numeric vector of length ",
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  direction <- as.numeric(direction)
  result$raw_direction_values <- hessian_probe_point_values(direction)
  result$raw_direction_norm <- hessian_probe_norm(direction)

  system_product <- as.vector(system_hessian %*% direction)
  residual <- system_product + gradient
  result$system_product_norm <- hessian_probe_norm(system_product)
  result$solve_residual_norm <- hessian_probe_norm(residual)
  result$solve_residual_scale <- max(
    1,
    result$system_product_norm,
    result$gradient_norm
  )
  result$solve_residual_relative_tolerance_term <-
    control$residual_relative_tolerance * result$solve_residual_scale
  result$solve_residual_threshold <-
    control$residual_absolute_tolerance +
    result$solve_residual_relative_tolerance_term
  result$solve_residual_values_finite <- all(is.finite(c(
    system_product,
    residual,
    result$system_product_norm,
    result$solve_residual_norm,
    result$solve_residual_scale,
    result$solve_residual_relative_tolerance_term,
    result$solve_residual_threshold
  )))
  result$solve_residual_pass <- result$solve_residual_values_finite &&
    result$solve_residual_norm <= result$solve_residual_threshold

  if (!result$solve_residual_values_finite) {
    result$solve_status <- "solve_residual_nonfinite"
    result$solve_error <- "solve residual evidence is nonfinite"
  } else if (!result$solve_residual_pass) {
    result$solve_status <- "solve_residual_tolerance_failed"
    result$solve_error <- "solve residual exceeds the scale-aware threshold"
  } else {
    result$solve_success <- TRUE
  }

  result$directional_slope <- sum(gradient * direction)
  result$quadratic_curvature <- sum(direction * system_product)
  result$predicted_decrease <- -(result$directional_slope +
    result$quadratic_curvature / 2)
  result$directional_slope_sign <- hessian_raw_newton_sign(
    result$directional_slope
  )
  result$raw_direction_is_descent <- if (is.finite(result$directional_slope)) {
    result$directional_slope < 0
  } else {
    NA
  }
  result$quadratic_curvature_sign <- hessian_raw_newton_sign(
    result$quadratic_curvature
  )
  result$predicted_decrease_sign <- hessian_raw_newton_sign(
    result$predicted_decrease
  )
  result$predicted_decrease_positive <- if (
    is.finite(result$predicted_decrease)
  ) {
    result$predicted_decrease > 0
  } else {
    NA
  }
  result$direction_diagnostics_complete <- all(is.finite(c(
    result$raw_direction_norm,
    result$directional_slope,
    result$quadratic_curvature,
    result$predicted_decrease
  )))
  if (result$solve_success) {
    result$solve_status <- if (result$direction_diagnostics_complete) {
      "solved"
    } else {
      "solved_direction_metrics_nonfinite"
    }
  }

  as.data.frame(result, stringsAsFactors = FALSE)
}

hessian_public_newton_control <- function() {
  list(
    comparison_absolute_tolerance = 0,
    comparison_relative_tolerance = 100 * sqrt(.Machine$double.eps)
  )
}

hessian_public_newton_provenance <- function(direction_reason) {
  provenance <- c(
    hessian_solve = "ordinary_cholesky_solve",
    cholesky_fallback = "cholesky_failure_steepest_descent",
    direction_check_fallback = "direction_check_failure_steepest_descent"
  )
  if (
    !is.character(direction_reason) ||
      length(direction_reason) != 1L ||
      is.na(direction_reason) ||
      !direction_reason %in% names(provenance)
  ) {
    return("evaluation_failure")
  }

  unname(provenance[[direction_reason]])
}

hessian_safeguarded_newton_branch_provenance <- function(direction_reason) {
  provenance <- c(
    hessian_solve = "cholesky_factor_solve",
    cholesky_fallback = "cholesky_failure_steepest_descent",
    direction_check_fallback = "direction_check_failure_steepest_descent"
  )
  if (
    !is.character(direction_reason) ||
      length(direction_reason) != 1L ||
      is.na(direction_reason) ||
      !direction_reason %in% names(provenance)
  ) {
    return("evaluation_failure")
  }

  unname(provenance[[direction_reason]])
}

hessian_safeguarded_newton_provenance <- function(
  direction_reason,
  ordinary_cholesky_success
) {
  if (
    !is.logical(ordinary_cholesky_success) ||
      length(ordinary_cholesky_success) != 1L ||
      is.na(ordinary_cholesky_success)
  ) {
    return("evaluation_failure")
  }

  if (identical(direction_reason, "hessian_solve")) {
    if (ordinary_cholesky_success) {
      return("ordinary_cholesky_solve")
    }
    return("eigenvalue_floor_repair_solve")
  }
  if (identical(direction_reason, "direction_check_fallback")) {
    if (ordinary_cholesky_success) {
      return("ordinary_cholesky_direction_check_fallback")
    }
    return("ambiguous_repair_or_failed_repair_direction_check_fallback")
  }
  if (identical(direction_reason, "cholesky_fallback")) {
    if (ordinary_cholesky_success) {
      return("unexpected_ordinary_cholesky_fallback")
    }
    return("eigenvalue_floor_repair_failed_steepest_descent")
  }

  "evaluation_failure"
}

hessian_public_newton_model_hessian <- function(hessian, dimension) {
  block_shape_ok <- is.matrix(hessian) &&
    is.numeric(hessian) &&
    nrow(hessian) > 0L &&
    ncol(hessian) == nrow(hessian) &&
    dimension %% nrow(hessian) == 0L
  finite <- block_shape_ok && all(is.finite(hessian))
  symmetric <- finite &&
    isTRUE(isSymmetric(
      hessian,
      tol = sqrt(.Machine$double.eps),
      check.attributes = FALSE
    ))
  model_hessian <- NULL
  if (symmetric) {
    block_dimension <- nrow(hessian)
    block_count <- dimension %/% block_dimension
    if (block_count == 1L) {
      model_hessian <- hessian
    } else {
      model_hessian <- kronecker(diag(block_count), hessian)
    }
  }

  list(
    shape_ok = block_shape_ok,
    finite = finite,
    symmetric = symmetric,
    value = model_hessian
  )
}

hessian_public_newton_raw_direction <- function(values, dimension) {
  if (
    !is.character(values) ||
      length(values) != 1L ||
      is.na(values) ||
      !nzchar(values)
  ) {
    return(NULL)
  }

  parsed <- suppressWarnings(as.numeric(strsplit(values, ";", fixed = TRUE)[[
    1
  ]]))
  if (length(parsed) != dimension || any(!is.finite(parsed))) {
    return(NULL)
  }
  parsed
}

probe_public_newton_direction <- function(
  fg,
  point,
  raw_newton,
  direction_constructor,
  control = hessian_public_newton_control(),
  try_safe_chol = FALSE
) {
  if (!is.list(fg) || !is.function(fg$gr) || !is.function(fg$hs)) {
    stop("fg must provide function callbacks fg$gr and fg$hs", call. = FALSE)
  }
  if (
    !is.numeric(point) ||
      !is.null(dim(point)) ||
      length(point) == 0L ||
      !all(is.finite(point))
  ) {
    stop("point must be a non-empty finite numeric vector", call. = FALSE)
  }
  if (!is.function(direction_constructor)) {
    stop("direction_constructor must be a function", call. = FALSE)
  }
  if (
    !is.logical(try_safe_chol) ||
      length(try_safe_chol) != 1L ||
      is.na(try_safe_chol)
  ) {
    stop("try_safe_chol must be TRUE or FALSE", call. = FALSE)
  }
  required_raw_fields <- c(
    "solve_success",
    "solve_status",
    "raw_direction_values"
  )
  if (
    !is.data.frame(raw_newton) ||
      nrow(raw_newton) != 1L ||
      !all(required_raw_fields %in% names(raw_newton))
  ) {
    stop("raw_newton must be one raw exact-Newton direction row", call. = FALSE)
  }

  required_controls <- c(
    "comparison_absolute_tolerance",
    "comparison_relative_tolerance"
  )
  if (
    !is.list(control) ||
      !all(required_controls %in% names(control)) ||
      any(!vapply(control[required_controls], is.numeric, logical(1))) ||
      any(vapply(control[required_controls], length, integer(1)) != 1L) ||
      any(!is.finite(unlist(control[required_controls]))) ||
      control$comparison_absolute_tolerance < 0 ||
      control$comparison_relative_tolerance <= 0
  ) {
    stop(
      paste(
        "public/raw comparison tolerances must be finite numeric scalars,",
        "with absolute tolerance nonnegative and relative tolerance positive"
      ),
      call. = FALSE
    )
  }

  raw_solve_success <- isTRUE(raw_newton$solve_success[[1L]])
  raw_solve_status <- as.character(raw_newton$solve_status[[1L]])
  raw_direction_values <- as.character(
    raw_newton$raw_direction_values[[1L]]
  )
  result <- list(
    public_direction_evaluation = paste(
      "gradient_seed_then",
      if (try_safe_chol) {
        "mize_newton_direction_try_safe_chol_true_calculate"
      } else {
        "mize_newton_direction_try_safe_chol_false_calculate"
      },
      sep = "_"
    ),
    public_direction_implementation = paste0(
      "mize:::newton_direction(try_safe_chol = ",
      if (try_safe_chol) "TRUE" else "FALSE",
      ")$calculate"
    ),
    public_direction_gr_calls = 0L,
    public_direction_hs_calls = 0L,
    benchmark_callback_counted = FALSE,
    public_direction_invoked = FALSE,
    public_direction_success = FALSE,
    public_direction_status = "not_evaluated",
    public_direction_error = "",
    gradient_callback_ok = NA,
    gradient_callback_error = "",
    gradient_shape_ok = NA,
    gradient_finite = NA,
    gradient_norm = NA_real_,
    hessian_callback_ok = NA,
    hessian_callback_error = "",
    hessian_shape_ok = NA,
    hessian_finite = NA,
    hessian_symmetric = NA,
    symmetrization = NA_character_,
    symmetrized_hessian_finite = NA,
    direction_reason = NA_character_,
    public_direction_provenance = "evaluation_failure",
    public_direction_values = "",
    public_direction_shape_ok = NA,
    public_direction_finite = NA,
    public_direction_norm = NA_real_,
    public_directional_slope = NA_real_,
    public_directional_slope_sign = "unavailable",
    public_direction_is_descent = NA,
    public_quadratic_curvature = NA_real_,
    public_quadratic_curvature_sign = "unavailable",
    public_predicted_decrease = NA_real_,
    public_predicted_decrease_sign = "unavailable",
    public_predicted_decrease_positive = NA,
    public_predicted_decrease_definition = "-(g' p + 0.5 p' H p)",
    public_direction_diagnostics_complete = FALSE,
    raw_solve_success = raw_solve_success,
    raw_solve_status = raw_solve_status,
    raw_direction_values = raw_direction_values,
    raw_comparison_available = FALSE,
    raw_comparison_basis = if (raw_solve_success) {
      "public_direction_unavailable"
    } else {
      paste0("raw_solve_status_", raw_solve_status)
    },
    raw_public_difference_norm = NA_real_,
    raw_public_comparison_scale = NA_real_,
    raw_public_comparison_absolute_tolerance = control$comparison_absolute_tolerance,
    raw_public_comparison_relative_tolerance = control$comparison_relative_tolerance,
    raw_public_comparison_absolute_tolerance_term = control$comparison_absolute_tolerance,
    raw_public_comparison_relative_tolerance_term = NA_real_,
    raw_public_comparison_threshold = NA_real_,
    raw_public_scaled_difference = NA_real_,
    raw_public_match = NA,
    raw_public_match_status = if (raw_solve_success) {
      "unavailable_public_direction"
    } else {
      "unavailable_raw_solve"
    }
  )
  if (try_safe_chol) {
    result <- append(
      result,
      list(
        ordinary_cholesky_attempted = FALSE,
        ordinary_cholesky_success = NA,
        ordinary_cholesky_error = "",
        safeguard_repair_floor = 1e-10,
        safeguard_repair_required = NA,
        safeguard_repair_attempted = FALSE,
        safeguard_repair_succeeded = NA,
        safeguard_provenance = "evaluation_failure"
      ),
      after = match("public_direction_provenance", names(result))
    )
  }

  gradient_result <- hessian_probe_callback(fg$gr, point)
  result$public_direction_gr_calls <- 1L
  result$gradient_callback_ok <- gradient_result$ok
  result$gradient_callback_error <- gradient_result$error
  if (!gradient_result$ok) {
    result$public_direction_status <- "gradient_callback_failed"
    result$public_direction_error <- gradient_result$error
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  dimension <- length(point)
  gradient <- gradient_result$value
  result$gradient_shape_ok <- is.numeric(gradient) &&
    is.null(dim(gradient)) &&
    length(gradient) == dimension
  result$gradient_finite <- result$gradient_shape_ok &&
    all(is.finite(gradient))
  if (!result$gradient_shape_ok || !result$gradient_finite) {
    result$public_direction_status <- "gradient_value_invalid"
    result$public_direction_error <- paste0(
      "gradient must be a finite numeric vector of length ",
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  gradient <- as.numeric(gradient)
  result$gradient_norm <- hessian_probe_norm(gradient)

  hessian_capture <- new.env(parent = emptyenv())
  hessian_capture$calls <- 0L
  hessian_capture$ok <- NA
  hessian_capture$value <- NULL
  hessian_capture$error <- ""
  captured_hessian <- function(x) {
    hessian_capture$calls <- hessian_capture$calls + 1L
    callback_result <- hessian_probe_callback(fg$hs, x)
    hessian_capture$ok <- callback_result$ok
    hessian_capture$value <- callback_result$value
    hessian_capture$error <- callback_result$error
    if (!callback_result$ok) {
      stop(callback_result$error, call. = FALSE)
    }
    callback_result$value
  }

  result$public_direction_invoked <- TRUE
  direction_result <- tryCatch(
    {
      direction <- direction_constructor(try_safe_chol = try_safe_chol)
      if (!is.list(direction) || !is.function(direction$calculate)) {
        stop(
          "newton_direction() did not return a calculable direction object",
          call. = FALSE
        )
      }
      direction$calculate(
        opt = list(cache = list(gr_curr = gradient, gr_curr_iter = 1L)),
        stage = list(),
        sub_stage = direction,
        par = point,
        fg = list(hs = captured_hessian),
        iter = 1L
      )
    },
    error = identity
  )
  result$public_direction_hs_calls <- hessian_capture$calls
  result$hessian_callback_ok <- hessian_capture$ok
  result$hessian_callback_error <- hessian_capture$error

  hessian_model <- hessian_public_newton_model_hessian(
    hessian_capture$value,
    dimension
  )
  result$hessian_shape_ok <- hessian_model$shape_ok
  result$hessian_finite <- hessian_model$finite
  result$hessian_symmetric <- hessian_model$symmetric
  if (try_safe_chol && hessian_model$symmetric) {
    result$ordinary_cholesky_attempted <- TRUE
    ordinary_cholesky <- tryCatch(
      chol(hessian_capture$value),
      error = identity
    )
    result$ordinary_cholesky_success <- !inherits(
      ordinary_cholesky,
      "error"
    )
    result$ordinary_cholesky_error <- if (result$ordinary_cholesky_success) {
      ""
    } else {
      conditionMessage(ordinary_cholesky)
    }
    result$safeguard_repair_required <-
      !result$ordinary_cholesky_success
    result$safeguard_repair_attempted <-
      result$safeguard_repair_required
  }

  if (inherits(direction_result, "error")) {
    result$public_direction_error <- conditionMessage(direction_result)
    result$public_direction_status <- if (
      identical(hessian_capture$ok, FALSE)
    ) {
      "hessian_callback_failed"
    } else if (
      identical(hessian_capture$ok, TRUE) &&
        (!hessian_model$shape_ok ||
          !hessian_model$finite ||
          !hessian_model$symmetric)
    ) {
      "hessian_value_invalid"
    } else {
      "direction_invocation_failed"
    }
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  valid_result <- is.list(direction_result) &&
    is.list(direction_result$sub_stage)
  if (!valid_result) {
    result$public_direction_status <- "direction_return_invalid"
    result$public_direction_error <- paste(
      "direction calculate must return a list containing sub_stage"
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  public_direction <- direction_result$sub_stage$value
  direction_reason <- direction_result$sub_stage$direction_reason
  result$direction_reason <- if (
    is.character(direction_reason) &&
      length(direction_reason) == 1L &&
      !is.na(direction_reason)
  ) {
    direction_reason
  } else {
    NA_character_
  }
  result$public_direction_provenance <- if (try_safe_chol) {
    hessian_safeguarded_newton_branch_provenance(result$direction_reason)
  } else {
    hessian_public_newton_provenance(result$direction_reason)
  }
  if (try_safe_chol) {
    result$safeguard_provenance <- hessian_safeguarded_newton_provenance(
      result$direction_reason,
      result$ordinary_cholesky_success
    )
    if (isTRUE(result$safeguard_repair_required)) {
      result$safeguard_repair_succeeded <- if (
        identical(result$direction_reason, "hessian_solve")
      ) {
        TRUE
      } else if (identical(result$direction_reason, "cholesky_fallback")) {
        FALSE
      } else {
        NA
      }
    }
  }
  result$public_direction_shape_ok <- is.numeric(public_direction) &&
    is.null(dim(public_direction)) &&
    length(public_direction) == dimension
  result$public_direction_finite <- result$public_direction_shape_ok &&
    all(is.finite(public_direction))
  if (!result$public_direction_shape_ok || !result$public_direction_finite) {
    result$public_direction_status <- "direction_return_invalid"
    result$public_direction_error <- paste0(
      if (try_safe_chol) {
        "safeguarded direction must be a finite numeric vector of length "
      } else {
        "public direction must be a finite numeric vector of length "
      },
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  if (identical(result$public_direction_provenance, "evaluation_failure")) {
    result$public_direction_status <- "direction_reason_invalid"
    result$public_direction_error <- paste(
      if (try_safe_chol) {
        paste(
          "safeguarded matrix-Hessian direction returned an unexpected",
          "direction_reason"
        )
      } else {
        paste(
          "public matrix-Hessian direction returned an unexpected",
          "direction_reason"
        )
      }
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  public_direction <- as.numeric(public_direction)
  result$public_direction_success <- TRUE
  result$public_direction_values <- hessian_probe_point_values(public_direction)
  result$public_direction_norm <- hessian_probe_norm(public_direction)
  result$public_directional_slope <- sum(gradient * public_direction)
  result$public_directional_slope_sign <- hessian_raw_newton_sign(
    result$public_directional_slope
  )
  result$public_direction_is_descent <- if (
    is.finite(result$public_directional_slope)
  ) {
    result$public_directional_slope < 0
  } else {
    NA
  }

  if (!is.null(hessian_model$value)) {
    symmetrized <- tryCatch(
      hessian_spectrum_symmetrize(hessian_model$value),
      error = identity
    )
    if (!inherits(symmetrized, "error")) {
      result$symmetrization <- symmetrized$method
      result$symmetrized_hessian_finite <- all(is.finite(symmetrized$value))
      if (result$symmetrized_hessian_finite) {
        system_product <- as.vector(symmetrized$value %*% public_direction)
        result$public_quadratic_curvature <- sum(
          public_direction * system_product
        )
        result$public_predicted_decrease <- -(result$public_directional_slope +
          result$public_quadratic_curvature / 2)
      }
    }
  }
  result$public_quadratic_curvature_sign <- hessian_raw_newton_sign(
    result$public_quadratic_curvature
  )
  result$public_predicted_decrease_sign <- hessian_raw_newton_sign(
    result$public_predicted_decrease
  )
  result$public_predicted_decrease_positive <- if (
    is.finite(result$public_predicted_decrease)
  ) {
    result$public_predicted_decrease > 0
  } else {
    NA
  }
  result$public_direction_diagnostics_complete <- all(is.finite(c(
    result$public_direction_norm,
    result$public_directional_slope,
    result$public_quadratic_curvature,
    result$public_predicted_decrease
  )))
  safeguard_provenance_ambiguous <- try_safe_chol &&
    identical(
      result$safeguard_provenance,
      "ambiguous_repair_or_failed_repair_direction_check_fallback"
    )
  result$public_direction_status <- if (safeguard_provenance_ambiguous) {
    if (result$public_direction_diagnostics_complete) {
      "evaluated_safeguard_provenance_ambiguous"
    } else {
      paste(
        "evaluated_safeguard_provenance_ambiguous",
        "direction_metrics_nonfinite",
        sep = "_"
      )
    }
  } else if (result$public_direction_diagnostics_complete) {
    "evaluated"
  } else {
    "evaluated_direction_metrics_nonfinite"
  }

  if (raw_solve_success) {
    raw_direction <- hessian_public_newton_raw_direction(
      raw_direction_values,
      dimension
    )
    if (is.null(raw_direction)) {
      result$raw_comparison_basis <- "raw_direction_values_invalid"
      result$raw_public_match_status <- "unavailable_raw_direction"
    } else {
      difference <- public_direction - raw_direction
      result$raw_public_difference_norm <- hessian_probe_norm(difference)
      raw_norm <- hessian_probe_norm(raw_direction)
      result$raw_public_comparison_scale <- max(
        1,
        result$public_direction_norm,
        raw_norm
      )
      result$raw_public_comparison_relative_tolerance_term <-
        control$comparison_relative_tolerance *
        result$raw_public_comparison_scale
      result$raw_public_comparison_threshold <-
        control$comparison_absolute_tolerance +
        result$raw_public_comparison_relative_tolerance_term
      result$raw_public_scaled_difference <-
        result$raw_public_difference_norm /
        result$raw_public_comparison_scale
      comparison_finite <- all(is.finite(c(
        difference,
        result$raw_public_difference_norm,
        result$raw_public_comparison_scale,
        result$raw_public_comparison_relative_tolerance_term,
        result$raw_public_comparison_threshold,
        result$raw_public_scaled_difference
      )))
      if (comparison_finite) {
        result$raw_comparison_available <- TRUE
        result$raw_comparison_basis <- "raw_and_public_directions_finite"
        result$raw_public_match <-
          result$raw_public_difference_norm <=
            result$raw_public_comparison_threshold
        result$raw_public_match_status <- if (result$raw_public_match) {
          "match"
        } else {
          "different"
        }
      } else {
        result$raw_comparison_basis <- "comparison_values_nonfinite"
        result$raw_public_match_status <- "unavailable_nonfinite_comparison"
      }
    }
  }

  as.data.frame(result, stringsAsFactors = FALSE)
}

probe_safeguarded_newton_direction <- function(
  fg,
  point,
  raw_newton,
  current_public_newton,
  direction_constructor,
  control = hessian_public_newton_control()
) {
  required_public_fields <- c(
    "public_direction_success",
    "public_direction_status",
    "public_direction_values"
  )
  if (
    !is.data.frame(current_public_newton) ||
      nrow(current_public_newton) != 1L ||
      !all(required_public_fields %in% names(current_public_newton))
  ) {
    stop(
      "current_public_newton must be one current-public direction row",
      call. = FALSE
    )
  }

  safeguarded <- probe_public_newton_direction(
    fg = fg,
    point = point,
    raw_newton = raw_newton,
    direction_constructor = direction_constructor,
    control = control,
    try_safe_chol = TRUE
  )
  safeguarded_names <- names(safeguarded)
  safeguarded_names <- sub(
    "^raw_public_",
    "raw_safeguarded_",
    safeguarded_names
  )
  safeguarded_names <- sub(
    "^public_",
    "safeguarded_",
    safeguarded_names
  )
  safeguarded_names[safeguarded_names == "raw_comparison_available"] <-
    "raw_safeguarded_comparison_available"
  safeguarded_names[safeguarded_names == "raw_comparison_basis"] <-
    "raw_safeguarded_comparison_basis"
  names(safeguarded) <- safeguarded_names
  safeguarded$raw_safeguarded_comparison_basis <- sub(
    "public",
    "safeguarded",
    safeguarded$raw_safeguarded_comparison_basis,
    fixed = TRUE
  )
  safeguarded$raw_safeguarded_match_status <- sub(
    "public",
    "safeguarded",
    safeguarded$raw_safeguarded_match_status,
    fixed = TRUE
  )

  current_public_success <- isTRUE(
    current_public_newton$public_direction_success[[1L]]
  )
  current_public_status <- as.character(
    current_public_newton$public_direction_status[[1L]]
  )
  current_public_values <- as.character(
    current_public_newton$public_direction_values[[1L]]
  )
  safeguarded_success <- isTRUE(
    safeguarded$safeguarded_direction_success[[1L]]
  )
  comparison <- list(
    current_public_direction_success = current_public_success,
    current_public_direction_status = current_public_status,
    current_public_direction_values = current_public_values,
    public_safeguarded_comparison_available = FALSE,
    public_safeguarded_comparison_basis = if (!current_public_success) {
      paste0("current_public_status_", current_public_status)
    } else if (!safeguarded_success) {
      "safeguarded_direction_unavailable"
    } else {
      "not_evaluated"
    },
    public_safeguarded_difference_norm = NA_real_,
    public_safeguarded_comparison_scale = NA_real_,
    public_safeguarded_comparison_absolute_tolerance = control$comparison_absolute_tolerance,
    public_safeguarded_comparison_relative_tolerance = control$comparison_relative_tolerance,
    public_safeguarded_comparison_absolute_tolerance_term = control$comparison_absolute_tolerance,
    public_safeguarded_comparison_relative_tolerance_term = NA_real_,
    public_safeguarded_comparison_threshold = NA_real_,
    public_safeguarded_scaled_difference = NA_real_,
    public_safeguarded_match = NA,
    public_safeguarded_match_status = if (!current_public_success) {
      "unavailable_current_public_direction"
    } else if (!safeguarded_success) {
      "unavailable_safeguarded_direction"
    } else {
      "not_evaluated"
    }
  )

  if (current_public_success && safeguarded_success) {
    dimension <- length(point)
    current_public_direction <- hessian_public_newton_raw_direction(
      current_public_values,
      dimension
    )
    safeguarded_direction <- hessian_public_newton_raw_direction(
      safeguarded$safeguarded_direction_values[[1L]],
      dimension
    )
    if (is.null(current_public_direction)) {
      comparison$public_safeguarded_comparison_basis <-
        "current_public_direction_values_invalid"
      comparison$public_safeguarded_match_status <-
        "unavailable_current_public_direction"
    } else if (is.null(safeguarded_direction)) {
      comparison$public_safeguarded_comparison_basis <-
        "safeguarded_direction_values_invalid"
      comparison$public_safeguarded_match_status <-
        "unavailable_safeguarded_direction"
    } else {
      difference <- safeguarded_direction - current_public_direction
      comparison$public_safeguarded_difference_norm <-
        hessian_probe_norm(difference)
      current_public_norm <- hessian_probe_norm(current_public_direction)
      safeguarded_norm <- hessian_probe_norm(safeguarded_direction)
      comparison$public_safeguarded_comparison_scale <- max(
        1,
        current_public_norm,
        safeguarded_norm
      )
      comparison$public_safeguarded_comparison_relative_tolerance_term <-
        control$comparison_relative_tolerance *
        comparison$public_safeguarded_comparison_scale
      comparison$public_safeguarded_comparison_threshold <-
        control$comparison_absolute_tolerance +
        comparison$public_safeguarded_comparison_relative_tolerance_term
      comparison$public_safeguarded_scaled_difference <-
        comparison$public_safeguarded_difference_norm /
        comparison$public_safeguarded_comparison_scale
      comparison_finite <- all(is.finite(c(
        difference,
        comparison$public_safeguarded_difference_norm,
        comparison$public_safeguarded_comparison_scale,
        comparison$public_safeguarded_comparison_relative_tolerance_term,
        comparison$public_safeguarded_comparison_threshold,
        comparison$public_safeguarded_scaled_difference
      )))
      if (comparison_finite) {
        comparison$public_safeguarded_comparison_available <- TRUE
        comparison$public_safeguarded_comparison_basis <-
          "current_public_and_safeguarded_directions_finite"
        comparison$public_safeguarded_match <-
          comparison$public_safeguarded_difference_norm <=
            comparison$public_safeguarded_comparison_threshold
        comparison$public_safeguarded_match_status <- if (
          comparison$public_safeguarded_match
        ) {
          "match"
        } else {
          "different"
        }
      } else {
        comparison$public_safeguarded_comparison_basis <-
          "comparison_values_nonfinite"
        comparison$public_safeguarded_match_status <-
          "unavailable_nonfinite_comparison"
      }
    }
  }

  data.frame(
    safeguarded,
    as.data.frame(comparison, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )
}

hessian_tn_direction_control <- function() {
  list(
    init = 0,
    exit_criterion = "curvature",
    preconditioner = "",
    max_inner_iterations = 40L,
    comparison_absolute_tolerance = 0,
    comparison_relative_tolerance = 100 * sqrt(.Machine$double.eps)
  )
}

hessian_tn_returned_inner_count <- function(value) {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      !is.null(dim(value)) ||
      !is.finite(value) ||
      value < 0 ||
      value != floor(value) ||
      value > .Machine$integer.max
  ) {
    return(NA_integer_)
  }
  as.integer(value)
}

hessian_tn_capture_inner_state <- function() {
  calls <- sys.calls()
  bd_indices <- which(vapply(
    calls,
    function(call) {
      head <- call[[1L]]
      is.symbol(head) && identical(as.character(head), "bd_approx")
    },
    logical(1)
  ))
  if (length(bd_indices) == 0L) {
    return(list(ok = FALSE, error = "bd_approx call frame was unavailable"))
  }

  bd_index <- tail(bd_indices, 1L)
  if (bd_index <= 1L) {
    return(list(ok = FALSE, error = "tn_inner_cg call frame was unavailable"))
  }
  bd_frame <- sys.frame(bd_index)
  inner_frame <- sys.frame(bd_index - 1L)
  bd_fields <- c("par", "dm", "gm", "h", "step")
  inner_fields <- c(
    "j",
    "zm",
    "rm",
    "dot_rm_ym",
    "eps",
    "exit_criterion"
  )
  if (
    !all(vapply(
      bd_fields,
      exists,
      logical(1),
      envir = bd_frame,
      inherits = FALSE
    )) ||
      !all(vapply(
        inner_fields,
        exists,
        logical(1),
        envir = inner_frame,
        inherits = FALSE
      ))
  ) {
    return(list(
      ok = FALSE,
      error = "required TN callback-frame state was unavailable"
    ))
  }

  state <- c(
    mget(bd_fields, envir = bd_frame, inherits = FALSE),
    mget(inner_fields, envir = inner_frame, inherits = FALSE)
  )
  state$ok <- TRUE
  state$error <- ""
  state
}

hessian_tn_inner_diagnostics <- function(
  records,
  gradient,
  returned_direction,
  direction_returned,
  returned_inner_count,
  control = hessian_tn_direction_control()
) {
  returned_inner_count <- hessian_tn_returned_inner_count(
    returned_inner_count
  )
  result <- list(
    tn_runtime_gradient_norm = sqrt(sum(gradient * gradient)),
    tn_inner_stop_reason = "observer_unavailable_or_inconsistent",
    tn_inner_stop_iteration = NA_integer_,
    tn_inner_hvp_attempts = length(records),
    tn_inner_accepted_updates = 0L,
    tn_inner_result_kind = "unavailable",
    tn_inner_result_values = "",
    tn_inner_directional_slope = NA_real_,
    tn_final_descent_fallback = NA,
    tn_inner_provenance_complete = FALSE,
    tn_inner_observer_error = "",
    tn_hvp_iteration_indices = "",
    tn_hvp_step_sizes = "",
    tn_hvp_curvatures = "",
    tn_hvp_step_coefficients = "",
    tn_hvp_residual_norms = "",
    tn_hvp_nonfinite_products_zeroed = 0L,
    tn_inner_count_matches_callbacks = NA
  )
  if (length(records) == 0L) {
    runtime_norm <- result$tn_runtime_gradient_norm
    if (!is.finite(runtime_norm)) {
      inner_direction <- rep(0, length(gradient))
      result$tn_inner_stop_reason <- "gradient_norm_unusable"
    } else if (runtime_norm == 0) {
      inner_direction <- rep(0, length(gradient))
      result$tn_inner_stop_reason <- "stationary_gradient"
    } else {
      result$tn_inner_observer_error <- paste(
        "nonstationary finite gradient returned without a TN HVP callback"
      )
      return(result)
    }
    result$tn_inner_result_kind <- "zero_vector"
    result$tn_inner_result_values <- hessian_probe_point_values(inner_direction)
    result$tn_inner_directional_slope <- sum(gradient * inner_direction)
    result$tn_final_descent_fallback <-
      result$tn_inner_directional_slope >= 0
    result$tn_inner_count_matches_callbacks <- if (
      is.na(returned_inner_count)
    ) {
      NA
    } else {
      returned_inner_count == 0L
    }
    expected_direction <- if (result$tn_final_descent_fallback) {
      -gradient
    } else {
      inner_direction
    }
    direction_matches <- direction_returned &&
      identical(as.numeric(returned_direction), as.numeric(expected_direction))
    result$tn_inner_provenance_complete <- isTRUE(direction_matches) &&
      isTRUE(result$tn_inner_count_matches_callbacks)
    if (!result$tn_inner_provenance_complete) {
      result$tn_inner_observer_error <- paste(
        "returned TN direction or inner count did not match the observed",
        "zero-HVP branch"
      )
    }
    return(result)
  }

  observer_ok <- vapply(
    records,
    function(record) {
      isTRUE(record$state$ok)
    },
    logical(1)
  )
  if (!all(observer_ok)) {
    errors <- vapply(
      records[!observer_ok],
      function(record) {
        record$state$error
      },
      character(1)
    )
    result$tn_inner_observer_error <- paste(unique(errors), collapse = "; ")
    return(result)
  }

  iterations <- vapply(
    records,
    function(record) {
      as.integer(record$state$j)
    },
    integer(1)
  )
  steps <- vapply(
    records,
    function(record) {
      as.numeric(record$state$h)
    },
    numeric(1)
  )
  result$tn_hvp_iteration_indices <- hessian_probe_point_values(iterations)
  result$tn_hvp_step_sizes <- hessian_probe_point_values(steps)
  curvatures <- rep(NA_real_, length(records))
  coefficients <- rep(NA_real_, length(records))
  residual_norms <- rep(NA_real_, length(records))
  expected_inner <- NULL
  stop_reason <- NULL
  stop_iteration <- NA_integer_
  accepted_updates <- 0L
  nonfinite_products_zeroed <- 0L

  for (index in seq_along(records)) {
    record <- records[[index]]
    state <- record$state
    shape_ok <- isTRUE(record$callback_ok) &&
      is.numeric(record$value) &&
      is.null(dim(record$value)) &&
      length(record$value) == length(gradient)
    if (!shape_ok) {
      stop_reason <- "callback_or_validation_failure"
      stop_iteration <- as.integer(state$j)
      break
    }
    if (
      !is.numeric(state$h) ||
        length(state$h) != 1L ||
        !is.finite(state$h) ||
        state$h <= 0 ||
        !identical(state$exit_criterion, control$exit_criterion)
    ) {
      result$tn_inner_observer_error <- paste(
        "captured TN finite-difference state was inconsistent"
      )
      break
    }

    product <- (as.numeric(record$value) - state$gm) / state$h
    if (any(!is.finite(product))) {
      product <- rep(0, length(state$dm))
      nonfinite_products_zeroed <- nonfinite_products_zeroed + 1L
    }
    curvature <- sum(state$dm * product)
    curvatures[[index]] <- curvature
    if (!is.finite(curvature) || curvature <= 0) {
      stop_reason <- if (is.finite(curvature)) {
        "nonpositive_curvature"
      } else {
        "curvature_unusable"
      }
      stop_iteration <- as.integer(state$j)
      expected_inner <- if (state$j == 0) -gradient else state$zm
      break
    }

    coefficient <- state$dot_rm_ym / curvature
    coefficients[[index]] <- coefficient
    if (!is.finite(coefficient)) {
      stop_reason <- "step_coefficient_unusable"
      stop_iteration <- as.integer(state$j)
      expected_inner <- if (state$j == 0) -gradient else state$zm
      break
    }
    candidate <- state$zm + coefficient * state$dm
    accepted_updates <- accepted_updates + 1L
    residual <- state$rm + coefficient * product
    residual_norm <- sqrt(sum(residual * residual))
    residual_norms[[index]] <- residual_norm
    if (!is.finite(residual_norm) || residual_norm < state$eps) {
      stop_reason <- if (is.finite(residual_norm)) {
        "residual_converged"
      } else {
        "residual_unusable"
      }
      stop_iteration <- as.integer(state$j)
      expected_inner <- candidate
      break
    }
    if (index == control$max_inner_iterations) {
      stop_reason <- "inner_iteration_limit"
      stop_iteration <- as.integer(state$j)
      expected_inner <- candidate
      break
    }
  }

  result$tn_hvp_curvatures <- hessian_probe_point_values(curvatures)
  result$tn_hvp_step_coefficients <- hessian_probe_point_values(coefficients)
  result$tn_hvp_residual_norms <- hessian_probe_point_values(residual_norms)
  result$tn_hvp_nonfinite_products_zeroed <- nonfinite_products_zeroed
  result$tn_inner_accepted_updates <- accepted_updates
  result$tn_inner_stop_iteration <- stop_iteration
  if (is.null(stop_reason)) {
    if (!nzchar(result$tn_inner_observer_error)) {
      result$tn_inner_observer_error <- paste(
        "captured TN predicates did not identify the realized stop branch"
      )
    }
    return(result)
  }
  result$tn_inner_stop_reason <- stop_reason
  if (identical(stop_reason, "callback_or_validation_failure")) {
    result$tn_inner_result_kind <- "unavailable"
    return(result)
  }
  if (is.null(expected_inner)) {
    result$tn_inner_observer_error <- "observed TN stop had no inner result"
    return(result)
  }

  initial_fallback_reasons <- c(
    "nonpositive_curvature",
    "curvature_unusable",
    "step_coefficient_unusable"
  )
  result$tn_inner_result_kind <- if (
    stop_reason %in% initial_fallback_reasons && stop_iteration == 0L
  ) {
    "steepest_descent"
  } else {
    "retained_cg_iterate"
  }
  result$tn_inner_result_values <- hessian_probe_point_values(expected_inner)
  result$tn_inner_directional_slope <- sum(gradient * expected_inner)
  descent_check <- result$tn_inner_directional_slope >= 0
  if (length(descent_check) != 1L || is.na(descent_check)) {
    result$tn_inner_observer_error <- paste(
      "the TN final descent predicate was not a logical scalar"
    )
    return(result)
  }
  result$tn_final_descent_fallback <- descent_check
  expected_direction <- if (descent_check) -gradient else expected_inner
  direction_matches <- direction_returned &&
    identical(as.numeric(returned_direction), as.numeric(expected_direction))
  callback_count <- sum(vapply(
    records,
    function(record) {
      isTRUE(record$callback_ok) &&
        is.numeric(record$value) &&
        is.null(dim(record$value)) &&
        length(record$value) == length(gradient)
    },
    logical(1)
  ))
  result$tn_inner_count_matches_callbacks <- if (is.na(returned_inner_count)) {
    NA
  } else {
    returned_inner_count == callback_count
  }
  result$tn_inner_provenance_complete <- isTRUE(direction_matches) &&
    isTRUE(result$tn_inner_count_matches_callbacks)
  if (!result$tn_inner_provenance_complete) {
    result$tn_inner_observer_error <- paste(
      "returned TN direction or inner count differed from observed predicates"
    )
  }
  result
}

hessian_tn_direction_comparison <- function(
  reference,
  reference_success,
  reference_status,
  reference_values,
  tn_success,
  tn_values,
  dimension,
  control = hessian_tn_direction_control()
) {
  comparison_prefix <- paste0(reference, "_tn_")
  unavailable_reference <- if (identical(reference, "raw")) {
    "unavailable_raw_solve"
  } else {
    paste0("unavailable_", reference, "_direction")
  }
  result <- list(
    reference_direction_success = reference_success,
    reference_direction_status = reference_status,
    reference_direction_values = reference_values,
    comparison_available = FALSE,
    comparison_basis = if (!reference_success) {
      if (identical(reference, "raw")) {
        paste0("raw_solve_status_", reference_status)
      } else {
        unavailable_reference
      }
    } else if (!tn_success) {
      "tn_direction_unavailable"
    } else {
      "comparison_not_evaluated"
    },
    difference_norm = NA_real_,
    comparison_scale = NA_real_,
    comparison_absolute_tolerance = control$comparison_absolute_tolerance,
    comparison_relative_tolerance = control$comparison_relative_tolerance,
    comparison_absolute_tolerance_term = control$comparison_absolute_tolerance,
    comparison_relative_tolerance_term = NA_real_,
    comparison_threshold = NA_real_,
    scaled_difference = NA_real_,
    match = NA,
    match_status = if (!reference_success) {
      unavailable_reference
    } else if (!tn_success) {
      "unavailable_tn_direction"
    } else {
      "unavailable_comparison"
    }
  )
  names(result)[1:3] <- if (identical(reference, "raw")) {
    c("raw_solve_success", "raw_solve_status", "raw_direction_values")
  } else {
    paste0(
      reference,
      c(
        "_direction_success",
        "_direction_status",
        "_direction_values"
      )
    )
  }
  names(result)[4:length(result)] <- paste0(
    comparison_prefix,
    names(result)[4:length(result)]
  )
  if (!reference_success || !tn_success) {
    return(result)
  }

  reference_direction <- hessian_public_newton_raw_direction(
    reference_values,
    dimension
  )
  tn_direction <- hessian_public_newton_raw_direction(tn_values, dimension)
  if (is.null(reference_direction)) {
    result[[paste0(comparison_prefix, "comparison_basis")]] <-
      paste0(reference, "_direction_values_invalid")
    result[[paste0(comparison_prefix, "match_status")]] <-
      unavailable_reference
    return(result)
  }
  if (is.null(tn_direction)) {
    result[[paste0(comparison_prefix, "comparison_basis")]] <-
      "tn_direction_values_invalid"
    result[[paste0(comparison_prefix, "match_status")]] <-
      "unavailable_tn_direction"
    return(result)
  }

  difference_norm <- hessian_probe_norm(tn_direction - reference_direction)
  comparison_scale <- max(
    1,
    hessian_probe_norm(reference_direction),
    hessian_probe_norm(tn_direction)
  )
  relative_term <- control$comparison_relative_tolerance * comparison_scale
  threshold <- control$comparison_absolute_tolerance + relative_term
  scaled_difference <- difference_norm / comparison_scale
  finite <- all(is.finite(c(
    difference_norm,
    comparison_scale,
    relative_term,
    threshold,
    scaled_difference
  )))
  if (!finite) {
    result[[paste0(comparison_prefix, "comparison_basis")]] <-
      "comparison_values_nonfinite"
    result[[paste0(comparison_prefix, "match_status")]] <-
      "unavailable_nonfinite_comparison"
    return(result)
  }

  match <- difference_norm <= threshold
  result[[paste0(comparison_prefix, "comparison_available")]] <- TRUE
  result[[paste0(comparison_prefix, "comparison_basis")]] <- paste0(
    reference,
    "_and_tn_directions_finite"
  )
  result[[paste0(comparison_prefix, "difference_norm")]] <- difference_norm
  result[[paste0(comparison_prefix, "comparison_scale")]] <- comparison_scale
  result[[paste0(
    comparison_prefix,
    "comparison_relative_tolerance_term"
  )]] <- relative_term
  result[[paste0(comparison_prefix, "comparison_threshold")]] <- threshold
  result[[paste0(comparison_prefix, "scaled_difference")]] <-
    scaled_difference
  result[[paste0(comparison_prefix, "match")]] <- match
  result[[paste0(comparison_prefix, "match_status")]] <- if (match) {
    "match"
  } else {
    "different"
  }
  result
}

probe_tn_direction <- function(
  fg,
  point,
  raw_newton,
  current_public_newton,
  safeguarded_newton,
  direction_constructor,
  control = hessian_tn_direction_control()
) {
  if (!is.list(fg) || !is.function(fg$gr) || !is.function(fg$hs)) {
    stop("fg must provide function callbacks fg$gr and fg$hs", call. = FALSE)
  }
  if (
    !is.numeric(point) ||
      !is.null(dim(point)) ||
      length(point) == 0L ||
      !all(is.finite(point))
  ) {
    stop("point must be a non-empty finite numeric vector", call. = FALSE)
  }
  if (!is.function(direction_constructor)) {
    stop("direction_constructor must be a function", call. = FALSE)
  }
  input_specs <- list(
    raw_newton = list(
      value = raw_newton,
      fields = c("solve_success", "solve_status", "raw_direction_values")
    ),
    current_public_newton = list(
      value = current_public_newton,
      fields = c(
        "public_direction_success",
        "public_direction_status",
        "public_direction_values"
      )
    ),
    safeguarded_newton = list(
      value = safeguarded_newton,
      fields = c(
        "safeguarded_direction_success",
        "safeguarded_direction_status",
        "safeguarded_direction_values"
      )
    )
  )
  for (input_name in names(input_specs)) {
    spec <- input_specs[[input_name]]
    if (
      !is.data.frame(spec$value) ||
        nrow(spec$value) != 1L ||
        !all(spec$fields %in% names(spec$value))
    ) {
      stop(input_name, " must be one compatible direction row", call. = FALSE)
    }
  }

  dimension <- length(point)
  result <- list(
    tn_direction_evaluation = paste(
      "gradient_seed_then_mize_tn_direction_default_calculate"
    ),
    tn_direction_implementation = paste(
      "mize:::tn_direction(init = 0, exit_criterion = 'curvature',",
      "preconditioner = '')$calculate"
    ),
    tn_init = "zero",
    tn_exit_criterion = control$exit_criterion,
    tn_preconditioner = "none",
    tn_max_inner_iterations = control$max_inner_iterations,
    tn_gradient_seed_calls = 0L,
    tn_hvp_gradient_calls = 0L,
    tn_total_gradient_calls = 0L,
    tn_algorithm_hs_calls = 0L,
    tn_reporting_hs_calls = 0L,
    benchmark_callback_counted = FALSE,
    tn_direction_invoked = FALSE,
    tn_direction_success = FALSE,
    tn_direction_status = "not_evaluated",
    tn_direction_error = "",
    gradient_seed_callback_ok = NA,
    gradient_seed_callback_error = "",
    gradient_seed_shape_ok = NA,
    gradient_seed_finite = NA,
    gradient_norm = NA_real_,
    tn_hvp_callbacks_returned = 0L,
    tn_hvp_shape_valid_results = 0L,
    tn_hvp_finite_results = 0L,
    tn_hvp_callback_errors = "",
    tn_returned_inner_gradient_count = NA_integer_,
    reporting_hessian_callback_ok = NA,
    reporting_hessian_callback_error = "",
    reporting_hessian_shape_ok = NA,
    reporting_hessian_finite = NA,
    reporting_hessian_symmetric = NA,
    symmetrization = NA_character_,
    symmetrized_hessian_finite = NA,
    tn_direction_values = "",
    tn_direction_shape_ok = NA,
    tn_direction_finite = NA,
    tn_direction_norm = NA_real_,
    tn_directional_slope = NA_real_,
    tn_directional_slope_sign = "unavailable",
    tn_direction_is_descent = NA,
    tn_quadratic_curvature = NA_real_,
    tn_quadratic_curvature_sign = "unavailable",
    tn_predicted_decrease = NA_real_,
    tn_predicted_decrease_sign = "unavailable",
    tn_predicted_decrease_positive = NA,
    tn_predicted_decrease_definition = "-(g' p + 0.5 p' H p)",
    tn_direction_metrics_complete = FALSE,
    tn_callback_accounting_complete = FALSE,
    tn_diagnostics_complete = FALSE,
    tn_runtime_gradient_norm = NA_real_,
    tn_inner_stop_reason = "not_evaluated",
    tn_inner_stop_iteration = NA_integer_,
    tn_inner_hvp_attempts = 0L,
    tn_inner_accepted_updates = 0L,
    tn_inner_result_kind = "unavailable",
    tn_inner_result_values = "",
    tn_inner_directional_slope = NA_real_,
    tn_final_descent_fallback = NA,
    tn_inner_provenance_complete = FALSE,
    tn_inner_observer_error = "",
    tn_hvp_iteration_indices = "",
    tn_hvp_step_sizes = "",
    tn_hvp_curvatures = "",
    tn_hvp_step_coefficients = "",
    tn_hvp_residual_norms = "",
    tn_hvp_nonfinite_products_zeroed = 0L,
    tn_inner_count_matches_callbacks = NA
  )
  comparisons <- c(
    hessian_tn_direction_comparison(
      reference = "raw",
      reference_success = isTRUE(raw_newton$solve_success[[1L]]),
      reference_status = as.character(raw_newton$solve_status[[1L]]),
      reference_values = as.character(raw_newton$raw_direction_values[[1L]]),
      tn_success = FALSE,
      tn_values = "",
      dimension = dimension,
      control = control
    ),
    hessian_tn_direction_comparison(
      reference = "current_public",
      reference_success = isTRUE(
        current_public_newton$public_direction_success[[1L]]
      ),
      reference_status = as.character(
        current_public_newton$public_direction_status[[1L]]
      ),
      reference_values = as.character(
        current_public_newton$public_direction_values[[1L]]
      ),
      tn_success = FALSE,
      tn_values = "",
      dimension = dimension,
      control = control
    ),
    hessian_tn_direction_comparison(
      reference = "safeguarded",
      reference_success = isTRUE(
        safeguarded_newton$safeguarded_direction_success[[1L]]
      ),
      reference_status = as.character(
        safeguarded_newton$safeguarded_direction_status[[1L]]
      ),
      reference_values = as.character(
        safeguarded_newton$safeguarded_direction_values[[1L]]
      ),
      tn_success = FALSE,
      tn_values = "",
      dimension = dimension,
      control = control
    )
  )
  result <- c(result, comparisons)

  seed <- hessian_probe_callback(fg$gr, point)
  result$tn_gradient_seed_calls <- 1L
  result$tn_total_gradient_calls <- 1L
  result$gradient_seed_callback_ok <- seed$ok
  result$gradient_seed_callback_error <- seed$error
  if (!seed$ok) {
    result$tn_direction_status <- "gradient_seed_callback_failed"
    result$tn_direction_error <- seed$error
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  gradient <- seed$value
  result$gradient_seed_shape_ok <- is.numeric(gradient) &&
    is.null(dim(gradient)) &&
    length(gradient) == dimension
  result$gradient_seed_finite <- result$gradient_seed_shape_ok &&
    all(is.finite(gradient))
  if (!result$gradient_seed_shape_ok || !result$gradient_seed_finite) {
    result$tn_direction_status <- "gradient_seed_value_invalid"
    result$tn_direction_error <- paste0(
      "gradient seed must be a finite numeric vector of length ",
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  gradient <- as.numeric(gradient)
  result$gradient_norm <- hessian_probe_norm(gradient)

  hvp_capture <- new.env(parent = emptyenv())
  hvp_capture$records <- list()
  captured_gradient <- function(x) {
    record <- list(
      state = hessian_tn_capture_inner_state(),
      point = x,
      callback_ok = FALSE,
      value = NULL,
      error = ""
    )
    callback_result <- hessian_probe_callback(fg$gr, x)
    record$callback_ok <- callback_result$ok
    record$value <- callback_result$value
    record$error <- callback_result$error
    hvp_capture$records[[length(hvp_capture$records) + 1L]] <- record
    if (!callback_result$ok) {
      stop(callback_result$error, call. = FALSE)
    }
    callback_result$value
  }

  result$tn_direction_invoked <- TRUE
  direction_result <- tryCatch(
    {
      direction <- direction_constructor(
        init = control$init,
        exit_criterion = control$exit_criterion,
        preconditioner = control$preconditioner
      )
      if (!is.list(direction) || !is.function(direction$calculate)) {
        stop(
          "tn_direction() did not return a calculable direction object",
          call. = FALSE
        )
      }
      direction$calculate(
        opt = list(
          cache = list(gr_curr = gradient, gr_curr_iter = 1L),
          counts = list(fn = 0L, gr = 0L),
          convergence = list(max_fn = Inf, max_gr = Inf, max_fg = Inf)
        ),
        stage = list(),
        sub_stage = direction,
        par = point,
        fg = list(gr = captured_gradient),
        iter = 1L
      )
    },
    error = identity
  )
  records <- hvp_capture$records
  result$tn_hvp_gradient_calls <- length(records)
  result$tn_total_gradient_calls <- 1L + length(records)
  result$tn_hvp_callbacks_returned <- sum(vapply(
    records,
    function(record) {
      isTRUE(record$callback_ok)
    },
    logical(1)
  ))
  result$tn_hvp_shape_valid_results <- sum(vapply(
    records,
    function(record) {
      isTRUE(record$callback_ok) &&
        is.numeric(record$value) &&
        is.null(dim(record$value)) &&
        length(record$value) == dimension
    },
    logical(1)
  ))
  result$tn_hvp_finite_results <- sum(vapply(
    records,
    function(record) {
      isTRUE(record$callback_ok) &&
        is.numeric(record$value) &&
        is.null(dim(record$value)) &&
        length(record$value) == dimension &&
        all(is.finite(record$value))
    },
    logical(1)
  ))
  callback_errors <- vapply(
    records,
    function(record) {
      if (isTRUE(record$callback_ok)) "" else record$error
    },
    character(1)
  )
  result$tn_hvp_callback_errors <- paste(
    unique(callback_errors[nzchar(callback_errors)]),
    collapse = "; "
  )

  valid_direction_result <- !inherits(direction_result, "error") &&
    is.list(direction_result) &&
    is.list(direction_result$sub_stage)
  returned_direction <- if (valid_direction_result) {
    direction_result$sub_stage$value
  } else {
    NULL
  }
  if (
    valid_direction_result &&
      is.list(direction_result$opt) &&
      is.list(direction_result$opt$counts)
  ) {
    result$tn_returned_inner_gradient_count <-
      hessian_tn_returned_inner_count(direction_result$opt$counts$gr)
  }
  inner <- hessian_tn_inner_diagnostics(
    records = records,
    gradient = gradient,
    returned_direction = returned_direction,
    direction_returned = valid_direction_result,
    returned_inner_count = result$tn_returned_inner_gradient_count,
    control = control
  )
  result[names(inner)] <- inner
  result$tn_callback_accounting_complete <-
    result$tn_gradient_seed_calls == 1L &&
    result$tn_total_gradient_calls ==
      result$tn_gradient_seed_calls + result$tn_hvp_gradient_calls &&
    result$tn_hvp_callbacks_returned <= result$tn_hvp_gradient_calls &&
    isTRUE(result$tn_inner_count_matches_callbacks)

  if (inherits(direction_result, "error")) {
    result$tn_direction_error <- conditionMessage(direction_result)
    result$tn_direction_status <- if (
      result$tn_hvp_callbacks_returned < result$tn_hvp_gradient_calls
    ) {
      "hvp_gradient_callback_failed"
    } else if (
      result$tn_hvp_shape_valid_results < result$tn_hvp_callbacks_returned
    ) {
      "hvp_gradient_value_invalid"
    } else {
      "direction_invocation_failed"
    }
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }
  if (!valid_direction_result) {
    result$tn_direction_status <- "direction_return_invalid"
    result$tn_direction_error <- paste(
      "TN calculate must return a list containing sub_stage"
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  result$tn_direction_shape_ok <- is.numeric(returned_direction) &&
    is.null(dim(returned_direction)) &&
    length(returned_direction) == dimension
  result$tn_direction_finite <- result$tn_direction_shape_ok &&
    all(is.finite(returned_direction))
  if (!result$tn_direction_shape_ok || !result$tn_direction_finite) {
    result$tn_direction_status <- "direction_return_invalid"
    result$tn_direction_error <- paste0(
      "TN direction must be a finite numeric vector of length ",
      dimension
    )
    return(as.data.frame(result, stringsAsFactors = FALSE))
  }

  tn_direction <- as.numeric(returned_direction)
  result$tn_direction_success <- TRUE
  result$tn_direction_values <- hessian_probe_point_values(tn_direction)
  result$tn_direction_norm <- hessian_probe_norm(tn_direction)
  result$tn_directional_slope <- sum(gradient * tn_direction)
  result$tn_directional_slope_sign <- hessian_raw_newton_sign(
    result$tn_directional_slope
  )
  result$tn_direction_is_descent <- if (
    is.finite(result$tn_directional_slope)
  ) {
    result$tn_directional_slope < 0
  } else {
    NA
  }

  hessian_result <- hessian_probe_callback(fg$hs, point)
  result$tn_reporting_hs_calls <- 1L
  result$reporting_hessian_callback_ok <- hessian_result$ok
  result$reporting_hessian_callback_error <- hessian_result$error
  hessian_model <- hessian_public_newton_model_hessian(
    hessian_result$value,
    dimension
  )
  result$reporting_hessian_shape_ok <- hessian_model$shape_ok
  result$reporting_hessian_finite <- hessian_model$finite
  result$reporting_hessian_symmetric <- hessian_model$symmetric
  if (!hessian_result$ok) {
    result$tn_direction_status <- "evaluated_reporting_hessian_callback_failed"
    result$tn_direction_error <- hessian_result$error
  } else if (is.null(hessian_model$value)) {
    result$tn_direction_status <- "evaluated_reporting_hessian_value_invalid"
    result$tn_direction_error <- paste(
      "reporting Hessian must be finite, symmetric, and dimension-compatible"
    )
  } else {
    symmetrized <- tryCatch(
      hessian_spectrum_symmetrize(hessian_model$value),
      error = identity
    )
    if (!inherits(symmetrized, "error")) {
      result$symmetrization <- symmetrized$method
      result$symmetrized_hessian_finite <- all(is.finite(symmetrized$value))
      if (result$symmetrized_hessian_finite) {
        system_product <- as.vector(symmetrized$value %*% tn_direction)
        result$tn_quadratic_curvature <- sum(tn_direction * system_product)
        result$tn_predicted_decrease <- -(result$tn_directional_slope +
          result$tn_quadratic_curvature / 2)
      }
    }
  }
  result$tn_quadratic_curvature_sign <- hessian_raw_newton_sign(
    result$tn_quadratic_curvature
  )
  result$tn_predicted_decrease_sign <- hessian_raw_newton_sign(
    result$tn_predicted_decrease
  )
  result$tn_predicted_decrease_positive <- if (
    is.finite(result$tn_predicted_decrease)
  ) {
    result$tn_predicted_decrease > 0
  } else {
    NA
  }
  result$tn_direction_metrics_complete <- all(is.finite(c(
    result$tn_direction_norm,
    result$tn_directional_slope,
    result$tn_quadratic_curvature,
    result$tn_predicted_decrease
  )))
  result$tn_diagnostics_complete <- result$tn_direction_metrics_complete &&
    result$tn_callback_accounting_complete &&
    result$tn_inner_provenance_complete
  if (identical(result$tn_direction_status, "not_evaluated")) {
    result$tn_direction_status <- if (!result$tn_inner_provenance_complete) {
      "evaluated_inner_provenance_incomplete"
    } else if (!result$tn_direction_metrics_complete) {
      "evaluated_direction_metrics_incomplete"
    } else {
      "evaluated"
    }
  }

  raw_comparison <- hessian_tn_direction_comparison(
    reference = "raw",
    reference_success = isTRUE(raw_newton$solve_success[[1L]]),
    reference_status = as.character(raw_newton$solve_status[[1L]]),
    reference_values = as.character(raw_newton$raw_direction_values[[1L]]),
    tn_success = result$tn_direction_success,
    tn_values = result$tn_direction_values,
    dimension = dimension,
    control = control
  )
  public_comparison <- hessian_tn_direction_comparison(
    reference = "current_public",
    reference_success = isTRUE(
      current_public_newton$public_direction_success[[1L]]
    ),
    reference_status = as.character(
      current_public_newton$public_direction_status[[1L]]
    ),
    reference_values = as.character(
      current_public_newton$public_direction_values[[1L]]
    ),
    tn_success = result$tn_direction_success,
    tn_values = result$tn_direction_values,
    dimension = dimension,
    control = control
  )
  safeguarded_comparison <- hessian_tn_direction_comparison(
    reference = "safeguarded",
    reference_success = isTRUE(
      safeguarded_newton$safeguarded_direction_success[[1L]]
    ),
    reference_status = as.character(
      safeguarded_newton$safeguarded_direction_status[[1L]]
    ),
    reference_values = as.character(
      safeguarded_newton$safeguarded_direction_values[[1L]]
    ),
    tn_success = result$tn_direction_success,
    tn_values = result$tn_direction_values,
    dimension = dimension,
    control = control
  )
  comparisons <- c(raw_comparison, public_comparison, safeguarded_comparison)
  result[names(comparisons)] <- comparisons
  as.data.frame(result, stringsAsFactors = FALSE)
}

hessian_probe_centered_difference <- function(plus, minus, step_size) {
  if (step_size > .Machine$double.xmax / 2) {
    return(plus / step_size / 2 - minus / step_size / 2)
  }

  denominator <- 2 * step_size
  if (denominator >= 1) {
    return(plus / denominator - minus / denominator)
  }

  (plus - minus) / denominator
}

probe_hessian_integrity <- function(
  fg,
  point,
  point_label = "probe",
  control = hessian_probe_control()
) {
  if (!is.list(fg) || !is.function(fg$gr) || !is.function(fg$hs)) {
    stop("fg must provide function callbacks fg$gr and fg$hs", call. = FALSE)
  }
  if (
    !is.numeric(point) ||
      !is.null(dim(point)) ||
      length(point) == 0L ||
      !all(is.finite(point))
  ) {
    stop("point must be a non-empty finite numeric vector", call. = FALSE)
  }

  required_controls <- c(
    "step_factor",
    "symmetry_relative_tolerance",
    "directional_absolute_tolerance",
    "directional_relative_tolerance",
    "directional_roundoff_multiplier"
  )
  if (
    !is.list(control) ||
      !all(required_controls %in% names(control)) ||
      any(!vapply(control[required_controls], is.numeric, logical(1))) ||
      any(vapply(control[required_controls], length, integer(1)) != 1L) ||
      any(!is.finite(unlist(control[required_controls]))) ||
      any(unlist(control[required_controls]) <= 0)
  ) {
    stop(
      "control values must be finite positive numeric scalars",
      call. = FALSE
    )
  }

  dimension <- length(point)
  directions <- hessian_probe_directions(point)
  hessian_result <- hessian_probe_callback(fg$hs, point)
  hessian <- hessian_result$value
  hessian_shape_ok <- hessian_result$ok &&
    is.matrix(hessian) &&
    is.numeric(hessian) &&
    identical(dim(hessian), c(dimension, dimension))
  hessian_finite <- hessian_shape_ok && all(is.finite(hessian))

  symmetry_absolute_residual <- NA_real_
  symmetry_scale <- NA_real_
  symmetry_scaled_residual <- NA_real_
  symmetry_threshold <- NA_real_
  symmetry_values_finite <- FALSE
  symmetry_pass <- FALSE
  if (hessian_finite) {
    symmetry_absolute_residual <- max(abs(hessian - t(hessian)))
    symmetry_scale <- max(1, max(abs(hessian)))
    symmetry_scaled_residual <- symmetry_absolute_residual / symmetry_scale
    symmetry_threshold <-
      control$symmetry_relative_tolerance * symmetry_scale
    symmetry_values_finite <- all(is.finite(c(
      symmetry_absolute_residual,
      symmetry_scale,
      symmetry_scaled_residual,
      symmetry_threshold
    )))
    symmetry_pass <- symmetry_values_finite &&
      symmetry_absolute_residual <= symmetry_threshold
  }

  rows <- lapply(names(directions), function(direction_name) {
    direction_spec <- directions[[direction_name]]
    direction <- direction_spec$direction
    direction_norm <- hessian_probe_norm(direction)
    step_size <- control$step_factor * direction_spec$direction_scale
    direction_geometry_finite <- all(is.finite(c(
      direction,
      direction_norm,
      direction_spec$direction_scale,
      step_size
    ))) &&
      step_size > 0
    point_plus <- point + step_size * direction
    point_minus <- point - step_size * direction
    probe_points_finite <- direction_geometry_finite &&
      all(is.finite(point_plus)) &&
      all(is.finite(point_minus))
    if (probe_points_finite) {
      gradient_plus_result <- hessian_probe_callback(fg$gr, point_plus)
      gradient_minus_result <- hessian_probe_callback(fg$gr, point_minus)
    } else {
      nonfinite_probe <- list(
        ok = FALSE,
        value = NULL,
        error = "centered probe point is nonfinite"
      )
      gradient_plus_result <- nonfinite_probe
      gradient_minus_result <- nonfinite_probe
    }
    gradient_plus <- gradient_plus_result$value
    gradient_minus <- gradient_minus_result$value
    gradient_shape_ok <- gradient_plus_result$ok &&
      gradient_minus_result$ok &&
      is.numeric(gradient_plus) &&
      is.null(dim(gradient_plus)) &&
      length(gradient_plus) == dimension &&
      is.numeric(gradient_minus) &&
      is.null(dim(gradient_minus)) &&
      length(gradient_minus) == dimension
    gradient_finite <- gradient_shape_ok &&
      all(is.finite(gradient_plus)) &&
      all(is.finite(gradient_minus))

    directional_absolute_residual <- NA_real_
    directional_scale <- NA_real_
    directional_scaled_residual <- NA_real_
    gradient_scale <- NA_real_
    directional_scaled_threshold <- NA_real_
    directional_roundoff_allowance <- NA_real_
    directional_threshold <- NA_real_
    directional_products_finite <- FALSE
    directional_values_finite <- FALSE
    directional_pass <- FALSE
    if (hessian_finite && gradient_finite && direction_geometry_finite) {
      analytic_product <- as.vector(hessian %*% direction)
      finite_difference_product <- hessian_probe_centered_difference(
        gradient_plus,
        gradient_minus,
        step_size
      )
      directional_products_finite <-
        all(is.finite(analytic_product)) &&
        all(is.finite(finite_difference_product))
      if (directional_products_finite) {
        directional_absolute_residual <- hessian_probe_norm(
          analytic_product - finite_difference_product
        )
        directional_scale <- max(
          1,
          hessian_probe_norm(analytic_product),
          hessian_probe_norm(finite_difference_product)
        )
        directional_scaled_residual <-
          directional_absolute_residual / directional_scale
        gradient_scale <- max(
          1,
          hessian_probe_norm(gradient_plus),
          hessian_probe_norm(gradient_minus)
        )
        directional_roundoff_allowance <-
          control$directional_roundoff_multiplier *
          .Machine$double.eps *
          gradient_scale /
          step_size
        directional_threshold <-
          control$directional_absolute_tolerance +
          control$directional_relative_tolerance * directional_scale +
          directional_roundoff_allowance
        directional_scaled_threshold <-
          directional_threshold / directional_scale
        directional_values_finite <- all(is.finite(c(
          directional_absolute_residual,
          directional_scale,
          directional_scaled_residual,
          gradient_scale,
          directional_roundoff_allowance,
          directional_threshold,
          directional_scaled_threshold
        )))
        directional_pass <- directional_values_finite &&
          directional_absolute_residual <= directional_threshold
      }
    }

    callback_errors <- c(
      if (!hessian_result$ok) paste0("hs: ", hessian_result$error),
      if (!gradient_plus_result$ok) {
        paste0("gr(+h): ", gradient_plus_result$error)
      },
      if (!gradient_minus_result$ok) {
        paste0("gr(-h): ", gradient_minus_result$error)
      }
    )

    data.frame(
      point = point_label,
      dimension = dimension,
      direction = direction_name,
      direction_norm = direction_norm,
      coordinate_scale_min = min(pmax(1, abs(point))),
      coordinate_scale_max = max(pmax(1, abs(point))),
      direction_scale = direction_spec$direction_scale,
      step_factor = control$step_factor,
      step_size = step_size,
      direction_geometry_finite = direction_geometry_finite,
      probe_points_finite = probe_points_finite,
      hessian_shape_ok = hessian_shape_ok,
      hessian_finite = hessian_finite,
      gradient_shape_ok = gradient_shape_ok,
      gradient_finite = gradient_finite,
      symmetry_absolute_residual = symmetry_absolute_residual,
      symmetry_scale = symmetry_scale,
      symmetry_scaled_residual = symmetry_scaled_residual,
      symmetry_relative_tolerance = control$symmetry_relative_tolerance,
      symmetry_threshold = symmetry_threshold,
      symmetry_values_finite = symmetry_values_finite,
      symmetry_pass = symmetry_pass,
      directional_absolute_residual = directional_absolute_residual,
      directional_scale = directional_scale,
      directional_scaled_residual = directional_scaled_residual,
      gradient_scale = gradient_scale,
      directional_absolute_tolerance = control$directional_absolute_tolerance,
      directional_relative_tolerance = control$directional_relative_tolerance,
      directional_roundoff_multiplier = control$directional_roundoff_multiplier,
      directional_roundoff_allowance = directional_roundoff_allowance,
      directional_threshold = directional_threshold,
      directional_scaled_threshold = directional_scaled_threshold,
      directional_products_finite = directional_products_finite,
      directional_values_finite = directional_values_finite,
      directional_pass = directional_pass,
      pass = direction_geometry_finite &&
        probe_points_finite &&
        hessian_shape_ok &&
        hessian_finite &&
        gradient_shape_ok &&
        gradient_finite &&
        symmetry_values_finite &&
        symmetry_pass &&
        directional_values_finite &&
        directional_pass,
      callback_error = paste(callback_errors, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  result$probe_pass <- all(result$pass)
  result$classification <- if (all(result$pass)) {
    "passes_at_probe_point"
  } else {
    "callback_or_adapter_evidence"
  }
  result
}

hessian_probe_spd_case <- function(n) {
  diagonal <- seq(0.5, 2, length.out = n)
  hessian <- diag(diagonal, nrow = n)
  linear <- seq(-1, 1, length.out = n)
  point <- seq(-0.75, 0.75, length.out = n)
  minimizer <- linear / diagonal

  list(
    name = paste0("spd-quadratic-n", n),
    source = "built-in",
    par = point,
    fg = list(
      fn = function(x) {
        drop(0.5 * t(x) %*% hessian %*% x - t(linear) %*% x)
      },
      gr = function(x) {
        as.vector(hessian %*% x - linear)
      },
      hs = function(x) {
        hessian
      }
    ),
    start = list(
      resolution = "fixed",
      requested_dimension = n,
      requested_dimension_accepted = TRUE,
      actual_dimension = n
    ),
    reference = list(
      xmin = minimizer,
      xmin_applicable = TRUE,
      xmin_basis = "exact_spd_minimizer"
    )
  )
}

hessian_probe_point_values <- function(point) {
  if (is.null(point)) {
    return("")
  }

  paste(
    format(point, digits = 17, scientific = TRUE, trim = TRUE),
    collapse = ";"
  )
}

hessian_probe_start_candidate <- function(case) {
  list(
    point_id = "resolved_start",
    point_label = "resolved_start",
    point_provenance = "adapter_resolved_start",
    point = case$par,
    candidate_eligible = TRUE,
    candidate_basis = "resolved_start",
    optimizer = NA_character_,
    method = NA_character_,
    line_search = NA_character_,
    step0 = NA_character_,
    callback_mode = NA_character_,
    requested_max_iter = NA_integer_,
    actual_iter = NA_integer_,
    termination_what = NA_character_,
    termination_value = NA_character_,
    reference_xmin_applicable = NA,
    reference_xmin_basis = NA_character_,
    reference_dimension_match = NA,
    error = ""
  )
}

hessian_probe_newton_candidates <- function(case, optimizer) {
  lapply(c(1L, 2L), function(max_iter) {
    optimizer_fg <- case$fg
    optimizer_fg$fg <- NULL
    result <- tryCatch(
      optimizer(
        par = case$par,
        fg = optimizer_fg,
        method = "NEWTON",
        line_search = "More-Thuente",
        step0 = NULL,
        max_iter = max_iter
      ),
      error = identity
    )
    failed <- inherits(result, "error")
    point <- if (failed) NULL else result$par
    point_valid <- !failed &&
      is.numeric(point) &&
      is.null(dim(point)) &&
      length(point) == length(case$par) &&
      all(is.finite(point))

    termination_what <- if (
      failed ||
        !is.list(result$terminate) ||
        is.null(result$terminate$what)
    ) {
      NA_character_
    } else {
      as.character(result$terminate$what)
    }
    termination_value <- if (
      failed ||
        !is.list(result$terminate) ||
        is.null(result$terminate$val)
    ) {
      NA_character_
    } else {
      paste(as.character(result$terminate$val), collapse = ",")
    }
    actual_iter <- if (
      failed ||
        length(result$iter) != 1L ||
        !is.numeric(result$iter) ||
        !is.finite(result$iter)
    ) {
      NA_integer_
    } else {
      as.integer(result$iter)
    }

    list(
      point_id = paste0("newton_max_iter_", max_iter),
      point_label = paste0("newton_max_iter_", max_iter),
      point_provenance = "mize_returned_parameters",
      point = if (point_valid) point else NULL,
      candidate_eligible = point_valid,
      candidate_basis = if (failed) {
        "optimizer_error"
      } else if (!point_valid) {
        "invalid_returned_parameters"
      } else {
        "mize_returned_parameters"
      },
      optimizer = "mize",
      method = "NEWTON",
      line_search = "More-Thuente",
      step0 = "default",
      callback_mode = "separate",
      requested_max_iter = max_iter,
      actual_iter = actual_iter,
      termination_what = termination_what,
      termination_value = termination_value,
      reference_xmin_applicable = NA,
      reference_xmin_basis = NA_character_,
      reference_dimension_match = NA,
      error = if (failed) conditionMessage(result) else ""
    )
  })
}

hessian_probe_reference_candidate <- function(case) {
  reference <- case$reference
  if (is.null(reference)) {
    reference <- list(
      xmin = NULL,
      xmin_applicable = FALSE,
      xmin_basis = "unavailable"
    )
  }

  xmin <- reference$xmin
  applicable <- reference$xmin_applicable
  basis <- reference$xmin_basis
  if (is.null(basis) || length(basis) != 1L || is.na(basis)) {
    basis <- "not_encoded"
  }
  dimension_match <- is.numeric(xmin) &&
    is.null(dim(xmin)) &&
    length(xmin) == length(case$par)
  point_valid <- dimension_match && all(is.finite(xmin))
  candidate_eligible <- identical(applicable, TRUE) && point_valid
  candidate_basis <- if (!identical(applicable, TRUE)) {
    basis
  } else if (!dimension_match) {
    "dimension_mismatch"
  } else if (!point_valid) {
    "invalid_reference_point"
  } else {
    "xmin_applicable_true_and_dimension_match"
  }

  list(
    point_id = "reference_xmin",
    point_label = "reference_xmin",
    point_provenance = "documented_xmin_reference",
    point = if (is.numeric(xmin) && is.null(dim(xmin))) xmin else NULL,
    candidate_eligible = candidate_eligible,
    candidate_basis = candidate_basis,
    optimizer = NA_character_,
    method = NA_character_,
    line_search = NA_character_,
    step0 = NA_character_,
    callback_mode = NA_character_,
    requested_max_iter = NA_integer_,
    actual_iter = NA_integer_,
    termination_what = NA_character_,
    termination_value = NA_character_,
    reference_xmin_applicable = applicable,
    reference_xmin_basis = basis,
    reference_dimension_match = dimension_match,
    error = ""
  )
}

hessian_probe_select_points <- function(case, optimizer) {
  candidates <- c(
    list(hessian_probe_start_candidate(case)),
    hessian_probe_newton_candidates(case, optimizer),
    list(hessian_probe_reference_candidate(case))
  )
  selected <- list()

  for (i in seq_along(candidates)) {
    candidate <- candidates[[i]]
    candidate$selection_status <- "excluded"
    candidate$selection_basis <- candidate$candidate_basis
    candidate$deduplicated_to <- NA_character_

    if (candidate$candidate_eligible) {
      duplicate_index <- which(vapply(
        selected,
        function(existing) {
          identical(as.numeric(existing$point), as.numeric(candidate$point))
        },
        logical(1)
      ))
      if (length(duplicate_index) > 0L) {
        candidate$selection_status <- "duplicate"
        candidate$selection_basis <- "identical_to_selected_point"
        candidate$deduplicated_to <- selected[[duplicate_index[[1L]]]]$point_id
      } else {
        candidate$selection_status <- "selected"
        candidate$selection_basis <- "unique_eligible_point"
        selected[[length(selected) + 1L]] <- candidate
      }
    }
    candidates[[i]] <- candidate
  }

  list(candidates = candidates, selected = selected)
}

hessian_probe_candidate_row <- function(case, candidate) {
  start <- case$start
  if (is.null(start)) {
    start <- list(
      resolution = "fixed",
      requested_dimension = length(case$par),
      requested_dimension_accepted = NA,
      actual_dimension = length(case$par)
    )
  }

  data.frame(
    case = case$name,
    source = case$source,
    start_resolution = start$resolution,
    requested_dimension = start$requested_dimension,
    requested_dimension_accepted = start$requested_dimension_accepted,
    actual_dimension = start$actual_dimension,
    point_id = candidate$point_id,
    point_provenance = candidate$point_provenance,
    point_values = hessian_probe_point_values(candidate$point),
    candidate_dimension = if (is.null(candidate$point)) {
      NA_integer_
    } else {
      length(candidate$point)
    },
    candidate_eligible = candidate$candidate_eligible,
    candidate_basis = candidate$candidate_basis,
    selection_status = candidate$selection_status,
    selection_basis = candidate$selection_basis,
    deduplicated_to = candidate$deduplicated_to,
    optimizer = candidate$optimizer,
    method = candidate$method,
    line_search = candidate$line_search,
    step0 = candidate$step0,
    callback_mode = candidate$callback_mode,
    requested_max_iter = candidate$requested_max_iter,
    actual_iter = candidate$actual_iter,
    termination_what = candidate$termination_what,
    termination_value = candidate$termination_value,
    reference_xmin_applicable = candidate$reference_xmin_applicable,
    reference_xmin_basis = candidate$reference_xmin_basis,
    reference_dimension_match = candidate$reference_dimension_match,
    error = candidate$error,
    stringsAsFactors = FALSE
  )
}

hessian_probe_selected_point <- function(case, candidate) {
  result <- probe_hessian_integrity(
    fg = case$fg,
    point = candidate$point,
    point_label = candidate$point_label
  )
  metadata <- hessian_probe_candidate_row(case, candidate)
  metadata <- metadata[rep(1L, nrow(result)), , drop = FALSE]
  rownames(metadata) <- NULL

  data.frame(
    metadata,
    result,
    stringsAsFactors = FALSE
  )
}

hessian_probe_spectrum_point <- function(case, candidate, integrity_result) {
  integrity_pass <- nrow(integrity_result) > 0L &&
    all(integrity_result$probe_pass)
  if (!integrity_pass) {
    stop(
      "Hessian spectra require a passing directional integrity result",
      call. = FALSE
    )
  }

  metadata <- hessian_probe_candidate_row(case, candidate)
  spectrum <- probe_hessian_spectrum(
    fg = case$fg,
    point = candidate$point
  )
  data.frame(
    metadata,
    integrity_gate = "probe_hessian_integrity",
    integrity_direction_count = nrow(integrity_result),
    integrity_probe_pass = TRUE,
    spectrum,
    stringsAsFactors = FALSE
  )
}

hessian_probe_raw_newton_point <- function(
  case,
  candidate,
  integrity_result,
  spectrum_result,
  linear_solver = base::solve,
  linear_solver_name = "base_solve"
) {
  integrity_pass <- nrow(integrity_result) > 0L &&
    all(integrity_result$probe_pass)
  if (!integrity_pass) {
    stop(
      "raw Newton directions require a passing directional integrity result",
      call. = FALSE
    )
  }
  if (!is.data.frame(spectrum_result) || nrow(spectrum_result) != 1L) {
    stop(
      "raw Newton directions require one Hessian spectrum row",
      call. = FALSE
    )
  }

  metadata <- hessian_probe_candidate_row(case, candidate)
  spectrum_evidence <- spectrum_result[,
    c(
      "spectral_classification",
      "singularity",
      "spectral_condition_estimate",
      "eigen_min",
      "eigen_max",
      "eigen_abs_min",
      "eigen_abs_max",
      "spectral_scale",
      "spectral_absolute_tolerance",
      "spectral_relative_tolerance",
      "spectral_absolute_tolerance_term",
      "spectral_relative_tolerance_term",
      "spectral_sign_tolerance",
      "inertia_positive",
      "inertia_zero",
      "inertia_negative",
      "inertia_unresolved"
    ),
    drop = FALSE
  ]
  direction <- probe_raw_newton_direction(
    fg = case$fg,
    point = candidate$point,
    spectrum = spectrum_result,
    linear_solver = linear_solver,
    linear_solver_name = linear_solver_name
  )
  data.frame(
    metadata,
    integrity_gate = "probe_hessian_integrity",
    integrity_direction_count = nrow(integrity_result),
    integrity_probe_pass = TRUE,
    spectrum_gate = "probe_hessian_spectrum",
    spectrum_evidence,
    direction,
    stringsAsFactors = FALSE
  )
}

hessian_probe_public_newton_point <- function(
  case,
  candidate,
  integrity_result,
  spectrum_result,
  raw_newton_result,
  direction_constructor
) {
  integrity_pass <- nrow(integrity_result) > 0L &&
    all(integrity_result$probe_pass)
  if (!integrity_pass) {
    stop(
      "public Newton directions require a passing directional integrity result",
      call. = FALSE
    )
  }
  if (!is.data.frame(spectrum_result) || nrow(spectrum_result) != 1L) {
    stop(
      "public Newton directions require one Hessian spectrum row",
      call. = FALSE
    )
  }
  if (!is.data.frame(raw_newton_result) || nrow(raw_newton_result) != 1L) {
    stop(
      "public Newton directions require one raw exact-Newton row",
      call. = FALSE
    )
  }

  metadata <- hessian_probe_candidate_row(case, candidate)
  spectrum_evidence <- spectrum_result[,
    c(
      "spectral_classification",
      "singularity",
      "spectral_condition_estimate",
      "eigen_min",
      "eigen_max",
      "eigen_abs_min",
      "eigen_abs_max",
      "spectral_scale",
      "spectral_absolute_tolerance",
      "spectral_relative_tolerance",
      "spectral_absolute_tolerance_term",
      "spectral_relative_tolerance_term",
      "spectral_sign_tolerance",
      "inertia_positive",
      "inertia_zero",
      "inertia_negative",
      "inertia_unresolved"
    ),
    drop = FALSE
  ]
  direction <- probe_public_newton_direction(
    fg = case$fg,
    point = candidate$point,
    raw_newton = raw_newton_result,
    direction_constructor = direction_constructor
  )
  data.frame(
    metadata,
    integrity_gate = "probe_hessian_integrity",
    integrity_direction_count = nrow(integrity_result),
    integrity_probe_pass = TRUE,
    spectrum_gate = "probe_hessian_spectrum",
    spectrum_evidence,
    raw_direction_gate = "probe_raw_newton_direction",
    direction,
    stringsAsFactors = FALSE
  )
}

hessian_probe_safeguarded_newton_point <- function(
  case,
  candidate,
  integrity_result,
  spectrum_result,
  raw_newton_result,
  current_public_newton_result,
  direction_constructor
) {
  integrity_pass <- nrow(integrity_result) > 0L &&
    all(integrity_result$probe_pass)
  if (!integrity_pass) {
    stop(
      "safeguarded Newton directions require passing integrity evidence",
      call. = FALSE
    )
  }
  if (!is.data.frame(spectrum_result) || nrow(spectrum_result) != 1L) {
    stop(
      "safeguarded Newton directions require one Hessian spectrum row",
      call. = FALSE
    )
  }
  if (!is.data.frame(raw_newton_result) || nrow(raw_newton_result) != 1L) {
    stop(
      "safeguarded Newton directions require one raw exact-Newton row",
      call. = FALSE
    )
  }
  if (
    !is.data.frame(current_public_newton_result) ||
      nrow(current_public_newton_result) != 1L
  ) {
    stop(
      "safeguarded Newton directions require one current-public row",
      call. = FALSE
    )
  }

  metadata <- hessian_probe_candidate_row(case, candidate)
  spectrum_evidence <- spectrum_result[,
    c(
      "spectral_classification",
      "singularity",
      "spectral_condition_estimate",
      "eigen_min",
      "eigen_max",
      "eigen_abs_min",
      "eigen_abs_max",
      "spectral_scale",
      "spectral_absolute_tolerance",
      "spectral_relative_tolerance",
      "spectral_absolute_tolerance_term",
      "spectral_relative_tolerance_term",
      "spectral_sign_tolerance",
      "inertia_positive",
      "inertia_zero",
      "inertia_negative",
      "inertia_unresolved"
    ),
    drop = FALSE
  ]
  direction <- probe_safeguarded_newton_direction(
    fg = case$fg,
    point = candidate$point,
    raw_newton = raw_newton_result,
    current_public_newton = current_public_newton_result,
    direction_constructor = direction_constructor
  )
  data.frame(
    metadata,
    integrity_gate = "probe_hessian_integrity",
    integrity_direction_count = nrow(integrity_result),
    integrity_probe_pass = TRUE,
    spectrum_gate = "probe_hessian_spectrum",
    spectrum_evidence,
    raw_direction_gate = "probe_raw_newton_direction",
    current_public_direction_gate = "probe_public_newton_direction",
    direction,
    stringsAsFactors = FALSE
  )
}

hessian_probe_tn_point <- function(
  case,
  candidate,
  integrity_result,
  spectrum_result,
  raw_newton_result,
  current_public_newton_result,
  safeguarded_newton_result,
  direction_constructor
) {
  integrity_pass <- nrow(integrity_result) > 0L &&
    all(integrity_result$probe_pass)
  if (!integrity_pass) {
    stop(
      "TN directions require a passing directional integrity result",
      call. = FALSE
    )
  }
  inputs <- list(
    spectrum = spectrum_result,
    raw = raw_newton_result,
    current_public = current_public_newton_result,
    safeguarded = safeguarded_newton_result
  )
  for (input_name in names(inputs)) {
    if (
      !is.data.frame(inputs[[input_name]]) || nrow(inputs[[input_name]]) != 1L
    ) {
      stop(
        "TN directions require one ",
        input_name,
        " evidence row",
        call. = FALSE
      )
    }
  }

  metadata <- hessian_probe_candidate_row(case, candidate)
  spectrum_evidence <- spectrum_result[,
    c(
      "spectral_classification",
      "singularity",
      "spectral_condition_estimate",
      "eigen_min",
      "eigen_max",
      "eigen_abs_min",
      "eigen_abs_max",
      "spectral_scale",
      "spectral_absolute_tolerance",
      "spectral_relative_tolerance",
      "spectral_absolute_tolerance_term",
      "spectral_relative_tolerance_term",
      "spectral_sign_tolerance",
      "inertia_positive",
      "inertia_zero",
      "inertia_negative",
      "inertia_unresolved"
    ),
    drop = FALSE
  ]
  direction <- probe_tn_direction(
    fg = case$fg,
    point = candidate$point,
    raw_newton = raw_newton_result,
    current_public_newton = current_public_newton_result,
    safeguarded_newton = safeguarded_newton_result,
    direction_constructor = direction_constructor
  )
  data.frame(
    metadata,
    integrity_gate = "probe_hessian_integrity",
    integrity_direction_count = nrow(integrity_result),
    integrity_probe_pass = TRUE,
    spectrum_gate = "probe_hessian_spectrum",
    spectrum_evidence,
    raw_direction_gate = "probe_raw_newton_direction",
    current_public_direction_gate = "probe_public_newton_direction",
    safeguarded_direction_gate = "probe_safeguarded_newton_direction",
    direction,
    stringsAsFactors = FALSE
  )
}

hessian_probe_extended_cases <- function(
  cases,
  optimizer,
  include_spectrum = FALSE,
  include_raw_newton = FALSE,
  include_public_newton = FALSE,
  include_safeguarded_newton = FALSE,
  include_tn = FALSE,
  raw_newton_solver = base::solve,
  raw_newton_solver_name = "base_solve",
  public_newton_direction = NULL,
  safeguarded_newton_direction = NULL,
  tn_direction = NULL
) {
  if (include_raw_newton && !include_spectrum) {
    stop(
      "raw Newton directions require Hessian spectrum evidence",
      call. = FALSE
    )
  }
  if (include_public_newton && !include_raw_newton) {
    stop(
      "public Newton directions require raw exact-Newton evidence",
      call. = FALSE
    )
  }
  if (include_public_newton && !is.function(public_newton_direction)) {
    stop(
      "public_newton_direction must be a function when public directions are selected",
      call. = FALSE
    )
  }
  if (include_safeguarded_newton && !include_public_newton) {
    stop(
      "safeguarded Newton directions require current-public evidence",
      call. = FALSE
    )
  }
  if (
    include_safeguarded_newton &&
      !is.function(safeguarded_newton_direction)
  ) {
    stop(
      paste(
        "safeguarded_newton_direction must be a function when safeguarded",
        "directions are selected"
      ),
      call. = FALSE
    )
  }
  if (include_tn && !include_safeguarded_newton) {
    stop(
      "TN directions require safeguarded exact-Newton evidence",
      call. = FALSE
    )
  }
  if (include_tn && !is.function(tn_direction)) {
    stop(
      "tn_direction must be a function when TN directions are selected",
      call. = FALSE
    )
  }
  selections <- lapply(
    cases,
    hessian_probe_select_points,
    optimizer = optimizer
  )
  results <- Map(
    function(case, selection) {
      do.call(
        rbind,
        lapply(selection$selected, hessian_probe_selected_point, case = case)
      )
    },
    cases,
    selections
  )
  spectra_by_case <- if (include_spectrum) {
    Map(
      function(case, selection, case_result) {
        point_spectra <- lapply(selection$selected, function(candidate) {
          integrity_result <- case_result[
            case_result$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          if (
            nrow(integrity_result) == 0L ||
              !all(integrity_result$probe_pass)
          ) {
            return(NULL)
          }
          hessian_probe_spectrum_point(
            case = case,
            candidate = candidate,
            integrity_result = integrity_result
          )
        })
        point_spectra <- Filter(Negate(is.null), point_spectra)
        if (length(point_spectra) == 0L) {
          return(NULL)
        }
        do.call(rbind, point_spectra)
      },
      cases,
      selections,
      results
    )
  } else {
    NULL
  }
  raw_newton_by_case <- if (include_raw_newton) {
    Map(
      function(
        case,
        selection,
        case_result,
        case_spectrum
      ) {
        point_directions <- lapply(selection$selected, function(candidate) {
          integrity_result <- case_result[
            case_result$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          if (
            nrow(integrity_result) == 0L ||
              !all(integrity_result$probe_pass)
          ) {
            return(NULL)
          }
          spectrum_result <- case_spectrum[
            case_spectrum$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          hessian_probe_raw_newton_point(
            case = case,
            candidate = candidate,
            integrity_result = integrity_result,
            spectrum_result = spectrum_result,
            linear_solver = raw_newton_solver,
            linear_solver_name = raw_newton_solver_name
          )
        })
        point_directions <- Filter(Negate(is.null), point_directions)
        if (length(point_directions) == 0L) {
          return(NULL)
        }
        do.call(rbind, point_directions)
      },
      cases,
      selections,
      results,
      spectra_by_case
    )
  } else {
    NULL
  }
  public_newton_by_case <- if (include_public_newton) {
    Map(
      function(
        case,
        selection,
        case_result,
        case_spectrum,
        case_raw_newton
      ) {
        point_directions <- lapply(selection$selected, function(candidate) {
          integrity_result <- case_result[
            case_result$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          if (
            nrow(integrity_result) == 0L ||
              !all(integrity_result$probe_pass)
          ) {
            return(NULL)
          }
          spectrum_result <- case_spectrum[
            case_spectrum$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          raw_newton_result <- case_raw_newton[
            case_raw_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          hessian_probe_public_newton_point(
            case = case,
            candidate = candidate,
            integrity_result = integrity_result,
            spectrum_result = spectrum_result,
            raw_newton_result = raw_newton_result,
            direction_constructor = public_newton_direction
          )
        })
        point_directions <- Filter(Negate(is.null), point_directions)
        if (length(point_directions) == 0L) {
          return(NULL)
        }
        do.call(rbind, point_directions)
      },
      cases,
      selections,
      results,
      spectra_by_case,
      raw_newton_by_case
    )
  } else {
    NULL
  }
  safeguarded_newton_by_case <- if (include_safeguarded_newton) {
    Map(
      function(
        case,
        selection,
        case_result,
        case_spectrum,
        case_raw_newton,
        case_public_newton
      ) {
        point_directions <- lapply(selection$selected, function(candidate) {
          integrity_result <- case_result[
            case_result$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          if (
            nrow(integrity_result) == 0L ||
              !all(integrity_result$probe_pass)
          ) {
            return(NULL)
          }
          spectrum_result <- case_spectrum[
            case_spectrum$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          raw_newton_result <- case_raw_newton[
            case_raw_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          current_public_newton_result <- case_public_newton[
            case_public_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          hessian_probe_safeguarded_newton_point(
            case = case,
            candidate = candidate,
            integrity_result = integrity_result,
            spectrum_result = spectrum_result,
            raw_newton_result = raw_newton_result,
            current_public_newton_result = current_public_newton_result,
            direction_constructor = safeguarded_newton_direction
          )
        })
        point_directions <- Filter(Negate(is.null), point_directions)
        if (length(point_directions) == 0L) {
          return(NULL)
        }
        do.call(rbind, point_directions)
      },
      cases,
      selections,
      results,
      spectra_by_case,
      raw_newton_by_case,
      public_newton_by_case
    )
  } else {
    NULL
  }
  tn_by_case <- if (include_tn) {
    Map(
      function(
        case,
        selection,
        case_result,
        case_spectrum,
        case_raw_newton,
        case_public_newton,
        case_safeguarded_newton
      ) {
        point_directions <- lapply(selection$selected, function(candidate) {
          integrity_result <- case_result[
            case_result$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          if (
            nrow(integrity_result) == 0L ||
              !all(integrity_result$probe_pass)
          ) {
            return(NULL)
          }
          spectrum_result <- case_spectrum[
            case_spectrum$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          raw_newton_result <- case_raw_newton[
            case_raw_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          current_public_newton_result <- case_public_newton[
            case_public_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          safeguarded_newton_result <- case_safeguarded_newton[
            case_safeguarded_newton$point_id == candidate$point_id,
            ,
            drop = FALSE
          ]
          hessian_probe_tn_point(
            case = case,
            candidate = candidate,
            integrity_result = integrity_result,
            spectrum_result = spectrum_result,
            raw_newton_result = raw_newton_result,
            current_public_newton_result = current_public_newton_result,
            safeguarded_newton_result = safeguarded_newton_result,
            direction_constructor = tn_direction
          )
        })
        point_directions <- Filter(Negate(is.null), point_directions)
        if (length(point_directions) == 0L) {
          return(NULL)
        }
        do.call(rbind, point_directions)
      },
      cases,
      selections,
      results,
      spectra_by_case,
      raw_newton_by_case,
      public_newton_by_case,
      safeguarded_newton_by_case
    )
  } else {
    NULL
  }
  manifests <- Map(
    function(case, selection) {
      do.call(
        rbind,
        lapply(selection$candidates, hessian_probe_candidate_row, case = case)
      )
    },
    cases,
    selections
  )

  coverage <- list(
    results = do.call(rbind, results),
    selection = do.call(rbind, manifests)
  )
  if (include_spectrum) {
    spectra <- Filter(Negate(is.null), spectra_by_case)
    coverage$spectrum <- if (length(spectra) == 0L) {
      data.frame()
    } else {
      do.call(rbind, spectra)
    }
  }
  if (include_raw_newton) {
    raw_newton <- Filter(Negate(is.null), raw_newton_by_case)
    coverage$raw_newton <- if (length(raw_newton) == 0L) {
      data.frame()
    } else {
      do.call(rbind, raw_newton)
    }
  }
  if (include_public_newton) {
    public_newton <- Filter(Negate(is.null), public_newton_by_case)
    coverage$public_newton <- if (length(public_newton) == 0L) {
      data.frame()
    } else {
      do.call(rbind, public_newton)
    }
  }
  if (include_safeguarded_newton) {
    safeguarded_newton <- Filter(
      Negate(is.null),
      safeguarded_newton_by_case
    )
    coverage$safeguarded_newton <- if (length(safeguarded_newton) == 0L) {
      data.frame()
    } else {
      do.call(rbind, safeguarded_newton)
    }
  }
  if (include_tn) {
    tn <- Filter(Negate(is.null), tn_by_case)
    coverage$tn <- if (length(tn) == 0L) {
      data.frame()
    } else {
      do.call(rbind, tn)
    }
  }
  coverage
}

hessian_probe_case <- function(case) {
  result <- probe_hessian_integrity(
    fg = case$fg,
    point = case$par,
    point_label = "resolved_start"
  )
  start <- case$start
  if (is.null(start)) {
    start <- list(
      resolution = "fixed",
      requested_dimension = length(case$par),
      requested_dimension_accepted = NA,
      actual_dimension = length(case$par)
    )
  }

  data.frame(
    case = case$name,
    source = case$source,
    start_resolution = start$resolution,
    requested_dimension = start$requested_dimension,
    requested_dimension_accepted = start$requested_dimension_accepted,
    actual_dimension = start$actual_dimension,
    result,
    stringsAsFactors = FALSE
  )
}

hessian_probe_split_arg <- function(value) {
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

hessian_probe_parse_int <- function(value, name) {
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  parsed
}

hessian_probe_resolved_output_path <- function(path, option) {
  if (
    !is.character(path) ||
      length(path) != 1L ||
      is.na(path) ||
      !nzchar(path)
  ) {
    stop(option, " must be a non-empty path", call. = FALSE)
  }

  expanded <- path.expand(path)
  link_target <- Sys.readlink(expanded)
  if (!is.na(link_target) && nzchar(link_target)) {
    stop(option, " must not be a symbolic link", call. = FALSE)
  }
  if (file.exists(expanded)) {
    return(normalizePath(expanded, winslash = "/", mustWork = TRUE))
  }

  resolved_parent <- normalizePath(
    dirname(expanded),
    winslash = "/",
    mustWork = FALSE
  )
  file.path(resolved_parent, basename(expanded))
}

hessian_probe_existing_file_identity <- function(path, option) {
  if (!file.exists(path)) {
    stop("existing-file identity requires an existing path", call. = FALSE)
  }

  stat_command <- Sys.which("stat")
  system_name <- Sys.info()[["sysname"]]
  stat_args <- if (
    nzchar(stat_command) &&
      system_name %in% c("Darwin", "FreeBSD", "OpenBSD", "NetBSD")
  ) {
    c("-f", "%d:%i", shQuote(path))
  } else if (nzchar(stat_command) && identical(system_name, "Linux")) {
    c("-Lc", "%d:%i", shQuote(path))
  } else {
    NULL
  }
  if (!is.null(stat_args)) {
    stat_output <- suppressWarnings(system2(
      stat_command,
      args = stat_args,
      stdout = TRUE,
      stderr = TRUE
    ))
    stat_status <- attr(stat_output, "status")
    if (is.null(stat_status)) {
      stat_status <- 0L
    }
    if (
      identical(stat_status, 0L) &&
        length(stat_output) == 1L &&
        grepl("^[0-9]+:[0-9]+$", stat_output)
    ) {
      return(stat_output)
    }
  }

  if (requireNamespace("fs", quietly = TRUE)) {
    info <- fs::file_info(path)
    device <- as.character(info$device_id[[1L]])
    inode <- as.character(info$inode[[1L]])
    if (!is.na(device) && !is.na(inode)) {
      return(paste(device, inode, sep = ":"))
    }
  }

  stop(
    "could not establish existing-file identity for ",
    option,
    "; install the fs package or provide unused output paths",
    call. = FALSE
  )
}

hessian_probe_directory_case_sensitive <- function(directory) {
  if (!dir.exists(directory)) {
    stop(
      "output parent directory does not exist: ",
      directory,
      call. = FALSE
    )
  }

  probe <- tempfile(
    pattern = ".mize-case-sensitivity-",
    tmpdir = directory
  )
  case_variant <- file.path(
    directory,
    toupper(basename(probe))
  )
  if (file.exists(case_variant)) {
    stop(
      "could not choose an unused output case-sensitivity probe",
      call. = FALSE
    )
  }
  if (!file.create(probe, showWarnings = FALSE)) {
    stop(
      "could not probe output-directory case sensitivity: ",
      directory,
      call. = FALSE
    )
  }
  on.exit(unlink(probe, force = TRUE), add = TRUE)

  if (!file.exists(case_variant)) {
    return(TRUE)
  }
  probe_identity <- hessian_probe_existing_file_identity(
    probe,
    "case-sensitivity probe"
  )
  variant_identity <- hessian_probe_existing_file_identity(
    case_variant,
    "case-sensitivity probe"
  )
  !identical(probe_identity, variant_identity)
}

hessian_probe_case_aliased_paths <- function(resolved) {
  case_aliased <- rep(FALSE, length(resolved))
  parents <- dirname(resolved)
  for (parent in unique(parents)) {
    indices <- which(parents == parent)
    folded <- tolower(basename(resolved[indices]))
    candidates <- duplicated(folded) |
      duplicated(folded, fromLast = TRUE)
    if (
      any(candidates) &&
        !hessian_probe_directory_case_sensitive(parent)
    ) {
      case_aliased[indices[candidates]] <- TRUE
    }
  }
  case_aliased
}

hessian_probe_preflight_output_paths <- function(paths) {
  active <- !vapply(paths, is.null, logical(1))
  paths <- paths[active]
  if (length(paths) == 0L) {
    return(invisible(character()))
  }
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("output paths must be named by their command-line options")
  }

  resolved <- Map(
    hessian_probe_resolved_output_path,
    path = paths,
    option = names(paths)
  )
  resolved <- vapply(resolved, identity, character(1))
  path_aliased <- duplicated(resolved) |
    duplicated(resolved, fromLast = TRUE)
  if (any(path_aliased)) {
    stop(
      "active output paths must resolve to different destinations: ",
      paste(names(resolved)[path_aliased], collapse = ", "),
      call. = FALSE
    )
  }

  case_aliased <- hessian_probe_case_aliased_paths(resolved)
  if (any(case_aliased)) {
    stop(
      "active output paths must not be case-equivalent aliases: ",
      paste(names(resolved)[case_aliased], collapse = ", "),
      call. = FALSE
    )
  }

  existing <- which(file.exists(resolved))
  if (length(existing) > 1L) {
    identities <- vapply(
      existing,
      function(index) {
        hessian_probe_existing_file_identity(
          resolved[[index]],
          names(resolved)[[index]]
        )
      },
      character(1)
    )
    identity_aliased <- duplicated(identities) |
      duplicated(identities, fromLast = TRUE)
    if (any(identity_aliased)) {
      stop(
        "active output paths must identify different existing files: ",
        paste(
          names(resolved)[existing[identity_aliased]],
          collapse = ", "
        ),
        call. = FALSE
      )
    }
  }

  invisible(resolved)
}

hessian_probe_usage <- function(status = 0L) {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/hessian-integrity-probe.R [options]",
      "",
      "Options:",
      "  --dimension N             SPD and requested variable-case dimension.",
      "  --funconstrain-cases NAMES Comma-separated funconstrain factories.",
      "  --no-funconstrain         Run only the exact SPD control.",
      "  --point-set MODE          resolved-start (default) or extended.",
      "  --selection-out PATH      Required manifest path for extended points.",
      "  --spectrum-out PATH       Optional spectrum CSV for extended points.",
      "  --direction-out PATH      Optional raw Newton CSV; requires spectrum.",
      paste(
        "  --public-direction-out PATH",
        "Optional current-public Newton CSV; requires raw."
      ),
      paste(
        "  --safeguarded-direction-out PATH",
        "Optional safeguarded Newton CSV; requires current-public."
      ),
      paste(
        "  --tn-direction-out PATH",
        "Optional TN direction CSV; requires safeguarded."
      ),
      "  --out PATH                Write CSV results to PATH instead of stdout.",
      "  --help                    Show this help.",
      sep = "\n"
    )
  )
  quit(save = "no", status = status)
}

hessian_probe_parse_args <- function(args) {
  config <- list(
    dimension = 5L,
    funconstrain_cases = c("rosen", "brown_bs", "var_dim", "chebyquad"),
    include_funconstrain = TRUE,
    point_set = "resolved-start",
    selection_out = NULL,
    spectrum_out = NULL,
    direction_out = NULL,
    public_direction_out = NULL,
    safeguarded_direction_out = NULL,
    tn_direction_out = NULL,
    out = NULL
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    read_value <- function() {
      if (i == length(args)) {
        stop(arg, " requires a value", call. = FALSE)
      }
      args[[i + 1L]]
    }

    if (arg %in% c("--help", "-h")) {
      hessian_probe_usage()
    } else if (arg == "--dimension") {
      config$dimension <- hessian_probe_parse_int(
        read_value(),
        "--dimension"
      )
      i <- i + 1L
    } else if (arg == "--funconstrain-cases") {
      config$funconstrain_cases <- hessian_probe_split_arg(read_value())
      i <- i + 1L
    } else if (arg == "--no-funconstrain") {
      config$include_funconstrain <- FALSE
    } else if (arg == "--point-set") {
      config$point_set <- tolower(read_value())
      i <- i + 1L
    } else if (arg == "--selection-out") {
      config$selection_out <- read_value()
      i <- i + 1L
    } else if (arg == "--spectrum-out") {
      config$spectrum_out <- read_value()
      i <- i + 1L
    } else if (arg == "--direction-out") {
      config$direction_out <- read_value()
      i <- i + 1L
    } else if (arg == "--public-direction-out") {
      config$public_direction_out <- read_value()
      i <- i + 1L
    } else if (arg == "--safeguarded-direction-out") {
      config$safeguarded_direction_out <- read_value()
      i <- i + 1L
    } else if (arg == "--tn-direction-out") {
      config$tn_direction_out <- read_value()
      i <- i + 1L
    } else if (arg == "--out") {
      config$out <- read_value()
      i <- i + 1L
    } else {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  if (!config$point_set %in% c("resolved-start", "extended")) {
    stop(
      "--point-set must be resolved-start or extended",
      call. = FALSE
    )
  }
  if (config$point_set == "extended" && is.null(config$selection_out)) {
    stop(
      "--selection-out is required when --point-set extended is selected",
      call. = FALSE
    )
  }
  if (
    config$point_set != "extended" &&
      (!is.null(config$spectrum_out) ||
        !is.null(config$direction_out) ||
        !is.null(config$public_direction_out) ||
        !is.null(config$safeguarded_direction_out) ||
        !is.null(config$tn_direction_out))
  ) {
    stop(
      paste(
        paste(
          "--spectrum-out, --direction-out, --public-direction-out, and",
          "--safeguarded-direction-out and --tn-direction-out"
        ),
        "require --point-set extended"
      ),
      call. = FALSE
    )
  }
  if (
    !is.null(config$direction_out) &&
      is.null(config$spectrum_out)
  ) {
    stop(
      "--direction-out requires --spectrum-out",
      call. = FALSE
    )
  }
  if (
    !is.null(config$public_direction_out) &&
      is.null(config$direction_out)
  ) {
    stop(
      "--public-direction-out requires --direction-out",
      call. = FALSE
    )
  }
  if (
    !is.null(config$safeguarded_direction_out) &&
      is.null(config$public_direction_out)
  ) {
    stop(
      "--safeguarded-direction-out requires --public-direction-out",
      call. = FALSE
    )
  }
  if (
    !is.null(config$tn_direction_out) &&
      is.null(config$safeguarded_direction_out)
  ) {
    stop(
      "--tn-direction-out requires --safeguarded-direction-out",
      call. = FALSE
    )
  }
  output_paths <- list("--out" = config$out)
  if (config$point_set == "extended") {
    output_paths[["--selection-out"]] <- config$selection_out
  }
  if (!is.null(config$spectrum_out)) {
    output_paths[["--spectrum-out"]] <- config$spectrum_out
  }
  if (!is.null(config$direction_out)) {
    output_paths[["--direction-out"]] <- config$direction_out
  }
  if (!is.null(config$public_direction_out)) {
    output_paths[["--public-direction-out"]] <- config$public_direction_out
  }
  if (!is.null(config$safeguarded_direction_out)) {
    output_paths[["--safeguarded-direction-out"]] <-
      config$safeguarded_direction_out
  }
  if (!is.null(config$tn_direction_out)) {
    output_paths[["--tn-direction-out"]] <- config$tn_direction_out
  }
  hessian_probe_preflight_output_paths(output_paths)

  config
}

hessian_probe_cases <- function(config) {
  spd <- hessian_probe_spd_case(config$dimension)
  cases <- setNames(list(spd), spd$name)
  if (config$include_funconstrain) {
    cases <- c(
      cases,
      optional_funconstrain_cases(
        problem_names = config$funconstrain_cases,
        n = config$dimension
      )
    )
  }
  cases
}

write_hessian_probe_results <- function(results, out) {
  if (is.null(out)) {
    utils::write.csv(results, file = stdout(), row.names = FALSE, na = "")
    return(invisible(NULL))
  }

  utils::write.csv(results, file = out, row.names = FALSE, na = "")
  message("Wrote Hessian-integrity results to ", out)
}

hessian_probe_load_mize <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
    return(getExportedValue("mize", "mize"))
  }
  if (requireNamespace("mize", quietly = TRUE)) {
    return(getExportedValue("mize", "mize"))
  }
  stop(
    "Install pkgload or install mize before selecting extended points.",
    call. = FALSE
  )
}

hessian_probe_load_public_newton_direction <- function() {
  namespace <- asNamespace("mize")
  if (!exists("newton_direction", envir = namespace, inherits = FALSE)) {
    stop(
      "The loaded mize package does not provide internal newton_direction().",
      call. = FALSE
    )
  }
  constructor <- get(
    "newton_direction",
    envir = namespace,
    inherits = FALSE
  )
  if (!is.function(constructor)) {
    stop(
      "The loaded mize newton_direction object is not a function.",
      call. = FALSE
    )
  }
  constructor
}

hessian_probe_load_tn_direction <- function() {
  namespace <- asNamespace("mize")
  if (!exists("tn_direction", envir = namespace, inherits = FALSE)) {
    stop(
      "The loaded mize package does not provide internal tn_direction().",
      call. = FALSE
    )
  }
  constructor <- get("tn_direction", envir = namespace, inherits = FALSE)
  if (!is.function(constructor)) {
    stop(
      "The loaded mize tn_direction object is not a function.",
      call. = FALSE
    )
  }
  constructor
}

write_hessian_probe_selection <- function(selection, out) {
  utils::write.csv(selection, file = out, row.names = FALSE, na = "")
  message("Wrote Hessian point-selection manifest to ", out)
}

write_hessian_probe_spectrum <- function(spectrum, out) {
  utils::write.csv(spectrum, file = out, row.names = FALSE, na = "")
  message("Wrote Hessian eigenspectrum results to ", out)
}

write_hessian_probe_raw_newton <- function(raw_newton, out) {
  utils::write.csv(raw_newton, file = out, row.names = FALSE, na = "")
  message("Wrote raw exact-Newton direction results to ", out)
}

write_hessian_probe_public_newton <- function(public_newton, out) {
  utils::write.csv(public_newton, file = out, row.names = FALSE, na = "")
  message("Wrote current-public exact-Newton direction results to ", out)
}

write_hessian_probe_safeguarded_newton <- function(
  safeguarded_newton,
  out
) {
  utils::write.csv(safeguarded_newton, file = out, row.names = FALSE, na = "")
  message("Wrote safeguarded exact-Newton direction results to ", out)
}

write_hessian_probe_tn <- function(tn, out) {
  utils::write.csv(tn, file = out, row.names = FALSE, na = "")
  message("Wrote TN direction results to ", out)
}

hessian_probe_main <- function() {
  config <- hessian_probe_parse_args(commandArgs(trailingOnly = TRUE))
  cases <- hessian_probe_cases(config)
  if (config$point_set == "resolved-start") {
    results <- do.call(rbind, lapply(cases, hessian_probe_case))
    write_hessian_probe_results(results, config$out)
    return(invisible(NULL))
  }

  optimizer <- hessian_probe_load_mize()
  public_newton_direction <- if (!is.null(config$public_direction_out)) {
    hessian_probe_load_public_newton_direction()
  } else {
    NULL
  }
  tn_direction <- if (!is.null(config$tn_direction_out)) {
    hessian_probe_load_tn_direction()
  } else {
    NULL
  }
  coverage <- hessian_probe_extended_cases(
    cases,
    optimizer = optimizer,
    include_spectrum = !is.null(config$spectrum_out),
    include_raw_newton = !is.null(config$direction_out),
    include_public_newton = !is.null(config$public_direction_out),
    include_safeguarded_newton = !is.null(config$safeguarded_direction_out),
    include_tn = !is.null(config$tn_direction_out),
    public_newton_direction = public_newton_direction,
    safeguarded_newton_direction = public_newton_direction,
    tn_direction = tn_direction
  )
  write_hessian_probe_results(coverage$results, config$out)
  write_hessian_probe_selection(coverage$selection, config$selection_out)
  if (!is.null(config$spectrum_out)) {
    write_hessian_probe_spectrum(coverage$spectrum, config$spectrum_out)
  }
  if (!is.null(config$direction_out)) {
    write_hessian_probe_raw_newton(
      coverage$raw_newton,
      config$direction_out
    )
  }
  if (!is.null(config$public_direction_out)) {
    write_hessian_probe_public_newton(
      coverage$public_newton,
      config$public_direction_out
    )
  }
  if (!is.null(config$safeguarded_direction_out)) {
    write_hessian_probe_safeguarded_newton(
      coverage$safeguarded_newton,
      config$safeguarded_direction_out
    )
  }
  if (!is.null(config$tn_direction_out)) {
    write_hessian_probe_tn(coverage$tn, config$tn_direction_out)
  }
}

if (sys.nframe() == 0L) {
  hessian_probe_main()
}
