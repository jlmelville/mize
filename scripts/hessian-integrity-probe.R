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
    )
  )
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
    } else if (arg == "--out") {
      config$out <- read_value()
      i <- i + 1L
    } else {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

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

hessian_probe_main <- function() {
  config <- hessian_probe_parse_args(commandArgs(trailingOnly = TRUE))
  cases <- hessian_probe_cases(config)
  results <- do.call(rbind, lapply(cases, hessian_probe_case))
  write_hessian_probe_results(results, config$out)
}

if (sys.nframe() == 0L) {
  hessian_probe_main()
}
