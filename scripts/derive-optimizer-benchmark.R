#!/usr/bin/env Rscript

source(file.path("scripts", "benchmark-optimizers.R"))

derivation_usage <- function(status = 0L) {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/derive-optimizer-benchmark.R [options]",
      "",
      "Options:",
      "  --summary PATH            Persisted resource-summary CSV input.",
      "  --cell-out PATH           Write deterministic cell medians.",
      "  --gate-out PATH           Write frozen material-cost gate rows.",
      "  --profile-case-out PATH   Write combined-mode profile/case summary.",
      "  --help                    Show this help.",
      sep = "\n"
    )
  )
  quit(save = "no", status = status)
}

parse_derivation_args <- function(args) {
  config <- list(
    summary = NULL,
    cell_out = NULL,
    gate_out = NULL,
    profile_case_out = NULL
  )
  options <- c(
    `--summary` = "summary",
    `--cell-out` = "cell_out",
    `--gate-out` = "gate_out",
    `--profile-case-out` = "profile_case_out"
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      derivation_usage()
    }
    if (!arg %in% names(options)) {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop(arg, " requires a value", call. = FALSE)
    }
    config[[options[[arg]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  missing <- names(config)[vapply(config, is.null, logical(1))]
  if (length(missing) > 0L) {
    stop(
      "Missing required option(s): ",
      paste(gsub("_", "-", paste0("--", missing)), collapse = ", "),
      call. = FALSE
    )
  }
  config
}

benchmark_read_resource_summary <- function(path) {
  summary <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c(
    "case",
    "profile",
    "line_search",
    "step0",
    "callback_mode",
    "rep",
    "final_f",
    "grad_norm",
    "iter",
    "n_fn_call",
    "n_gr_call",
    "n_fg_call",
    "n_hs_call",
    "elapsed_sec",
    "objective_reduction_hit",
    "gradient_reduction_hit",
    "reference_gap_target_applicable",
    "reference_gap_hit"
  )
  missing <- setdiff(required, names(summary))
  if (length(missing) > 0L) {
    stop(
      "resource summary lacks required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(summary) == 0L) {
    stop("resource summary must contain at least one row", call. = FALSE)
  }

  grouping <- c("case", "profile", "line_search", "step0", "callback_mode")
  missing_group <- vapply(
    summary[grouping],
    function(value) any(is.na(value) | !nzchar(as.character(value))),
    logical(1)
  )
  if (any(missing_group)) {
    stop(
      "resource summary has unavailable grouping values in: ",
      paste(grouping[missing_group], collapse = ", "),
      call. = FALSE
    )
  }
  if (any(is.na(summary$rep))) {
    stop("resource summary has unavailable repetitions", call. = FALSE)
  }

  expected_profiles <- c(
    "mize-newton",
    "mize-newton-safeguarded",
    "mize-tn",
    "mize-tn-l-bfgs",
    "mize-bfgs",
    "mize-l-bfgs",
    "mize-sr1",
    "mize-cg-pr+",
    "mize-cg-hz+",
    "mize-cg-pr+-l-bfgs"
  )
  expected_searches <- c(
    "More-Thuente",
    "Rasmussen",
    "Schmidt",
    "Hager-Zhang"
  )
  expected_callbacks <- c("combined", "separate")
  if (!setequal(summary$profile, expected_profiles)) {
    stop("resource summary does not contain the frozen Packet 4 profiles")
  }
  if (!setequal(summary$line_search, expected_searches)) {
    stop("resource summary does not contain the frozen Packet 4 searches")
  }
  if (!setequal(summary$callback_mode, expected_callbacks)) {
    stop("resource summary does not contain both Packet 4 callback modes")
  }
  if (!identical(unique(summary$step0), "default")) {
    stop("resource summary does not use only Packet 4 step0 = default")
  }

  run_key <- do.call(
    paste,
    c(summary[c(grouping, "rep")], sep = "\r")
  )
  if (anyDuplicated(run_key)) {
    stop("resource summary has duplicate run keys", call. = FALSE)
  }
  cell_key <- do.call(paste, c(summary[grouping], sep = "\r"))
  repetitions <- split(summary$rep, cell_key)
  complete_repetitions <- vapply(
    repetitions,
    function(value) identical(sort(as.integer(value)), 1:3),
    logical(1)
  )
  if (!all(complete_repetitions)) {
    stop("resource summary does not have repetitions 1, 2, and 3 per cell")
  }
  summary
}

benchmark_logical_code <- function(value) {
  if (length(value) != 1L) {
    stop("target state must be scalar", call. = FALSE)
  }
  if (is.na(value)) {
    return("u")
  }
  if (isTRUE(value)) "1" else "0"
}

benchmark_constant_logical <- function(values, label) {
  codes <- vapply(values, benchmark_logical_code, character(1))
  codes <- unique(codes)
  if (length(codes) != 1L) {
    stop(label, " differs across repetitions", call. = FALSE)
  }
  switch(codes, `1` = TRUE, `0` = FALSE, u = NA)
}

benchmark_cell_median_row <- function(rows) {
  objective_hit <- benchmark_constant_logical(
    rows$objective_reduction_hit,
    "objective target state"
  )
  gradient_hit <- benchmark_constant_logical(
    rows$gradient_reduction_hit,
    "gradient target state"
  )
  reference_applicable <- benchmark_constant_logical(
    rows$reference_gap_target_applicable,
    "reference-target applicability"
  )
  reference_hit <- benchmark_constant_logical(
    rows$reference_gap_hit,
    "reference target state"
  )
  reference_code <- if (isFALSE(reference_applicable)) {
    "x"
  } else {
    benchmark_logical_code(reference_hit)
  }

  data.frame(
    case = rows$case[[1L]],
    profile = rows$profile[[1L]],
    line_search = rows$line_search[[1L]],
    callback_mode = rows$callback_mode[[1L]],
    objective_hit = objective_hit,
    gradient_hit = gradient_hit,
    reference_applicable = reference_applicable,
    reference_hit = reference_hit,
    n_fn_call = stats::median(rows$n_fn_call),
    n_gr_call = stats::median(rows$n_gr_call),
    n_fg_call = stats::median(rows$n_fg_call),
    n_hs_call = stats::median(rows$n_hs_call),
    elapsed_sec = stats::median(rows$elapsed_sec),
    final_f = stats::median(rows$final_f),
    grad_norm = stats::median(rows$grad_norm),
    iter = stats::median(rows$iter),
    target_signature = paste(
      benchmark_logical_code(objective_hit),
      benchmark_logical_code(gradient_hit),
      reference_code,
      sep = "/"
    ),
    any_target = isTRUE(objective_hit) ||
      isTRUE(gradient_hit) ||
      (isTRUE(reference_applicable) && isTRUE(reference_hit)),
    stringsAsFactors = FALSE
  )
}

benchmark_cell_medians <- function(summary) {
  grouping <- c("case", "profile", "line_search", "callback_mode", "step0")
  ordering <- do.call(order, c(summary[grouping], list(summary$rep)))
  summary <- summary[ordering, , drop = FALSE]
  keys <- do.call(paste, c(summary[grouping], sep = "\r"))
  groups <- split(seq_len(nrow(summary)), keys)
  cells <- lapply(
    groups,
    function(indices) {
      benchmark_cell_median_row(summary[indices, , drop = FALSE])
    }
  )
  cells <- do.call(rbind, cells)
  rownames(cells) <- NULL
  ordering <- order(
    cells$case,
    cells$profile,
    cells$line_search,
    cells$callback_mode
  )
  cells[ordering, , drop = FALSE]
}

benchmark_material_cost_gates <- function(cells) {
  group_names <- c("case", "line_search", "callback_mode")
  keys <- do.call(paste, c(cells[group_names], sep = "\r"))
  groups <- split(seq_len(nrow(cells)), keys)
  resources <- c(
    "n_fn_call",
    "n_gr_call",
    "n_fg_call",
    "n_hs_call",
    "elapsed_sec"
  )
  gates <- list()

  for (indices in groups) {
    group <- cells[indices, , drop = FALSE]
    for (candidate in seq_len(nrow(group))) {
      if (!isTRUE(group$any_target[[candidate]])) {
        next
      }
      for (comparator in seq_len(nrow(group))) {
        if (
          candidate == comparator ||
            group$target_signature[[candidate]] !=
              group$target_signature[[comparator]]
        ) {
          next
        }
        for (resource in resources) {
          candidate_value <- group[[resource]][[candidate]]
          comparator_value <- group[[resource]][[comparator]]
          minimum_difference <- if (resource == "elapsed_sec") 0.05 else 50
          if (
            is.finite(candidate_value) &&
              is.finite(comparator_value) &&
              comparator_value > 0 &&
              candidate_value > 10 * comparator_value &&
              candidate_value - comparator_value >= minimum_difference
          ) {
            gates[[length(gates) + 1L]] <- data.frame(
              case = group$case[[candidate]],
              line_search = group$line_search[[candidate]],
              callback_mode = group$callback_mode[[candidate]],
              target_signature = group$target_signature[[candidate]],
              resource = resource,
              candidate_profile = group$profile[[candidate]],
              candidate_median = candidate_value,
              comparator_profile = group$profile[[comparator]],
              comparator_median = comparator_value,
              ratio = candidate_value / comparator_value,
              difference = candidate_value - comparator_value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(gates) == 0L) {
    return(data.frame(
      case = character(),
      line_search = character(),
      callback_mode = character(),
      target_signature = character(),
      resource = character(),
      candidate_profile = character(),
      candidate_median = numeric(),
      comparator_profile = character(),
      comparator_median = numeric(),
      ratio = numeric(),
      difference = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  gates <- do.call(rbind, gates)
  rownames(gates) <- NULL
  ordering <- order(
    gates$case,
    gates$line_search,
    gates$callback_mode,
    gates$target_signature,
    gates$resource,
    gates$candidate_profile,
    gates$comparator_profile
  )
  gates[ordering, , drop = FALSE]
}

benchmark_profile_case_summary <- function(cells) {
  cells <- cells[cells$callback_mode == "combined", , drop = FALSE]
  keys <- paste(cells$case, cells$profile, sep = "\r")
  groups <- split(seq_len(nrow(cells)), keys)
  rows <- lapply(groups, function(indices) {
    group <- cells[indices, , drop = FALSE]
    data.frame(
      case = group$case[[1L]],
      profile = group$profile[[1L]],
      objective_hits = sum(group$objective_hit %in% TRUE),
      gradient_hits = sum(group$gradient_hit %in% TRUE),
      reference_hits = sum(
        group$reference_applicable %in% TRUE & group$reference_hit %in% TRUE
      ),
      best_f = min(group$final_f),
      worst_f = max(group$final_f),
      best_grad_norm = min(group$grad_norm),
      worst_grad_norm = max(group$grad_norm),
      median_fn = stats::median(group$n_fn_call),
      median_gr = stats::median(group$n_gr_call),
      median_fg = stats::median(group$n_fg_call),
      median_hs = stats::median(group$n_hs_call),
      median_elapsed_sec = stats::median(group$elapsed_sec),
      median_iter = stats::median(group$iter),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  ordering <- order(result$case, result$profile)
  result[ordering, , drop = FALSE]
}

write_decision_artifact <- function(value, path, label) {
  utils::write.csv(
    value,
    file = path,
    row.names = FALSE,
    na = "",
    eol = "\n",
    fileEncoding = "UTF-8"
  )
  message("Wrote ", label, " to ", path)
}

derive_optimizer_benchmark <- function(config) {
  benchmark_preflight_output_paths(c(
    `--summary` = config$summary,
    `--cell-out` = config$cell_out,
    `--gate-out` = config$gate_out,
    `--profile-case-out` = config$profile_case_out
  ))
  summary <- benchmark_read_resource_summary(config$summary)
  cells <- benchmark_cell_medians(summary)
  gates <- benchmark_material_cost_gates(cells)
  profile_cases <- benchmark_profile_case_summary(cells)
  write_decision_artifact(cells, config$cell_out, "cell medians")
  write_decision_artifact(gates, config$gate_out, "material-cost gates")
  write_decision_artifact(
    profile_cases,
    config$profile_case_out,
    "profile/case summary"
  )
  invisible(list(cells = cells, gates = gates, profile_cases = profile_cases))
}

derivation_main <- function() {
  config <- parse_derivation_args(commandArgs(trailingOnly = TRUE))
  derive_optimizer_benchmark(config)
}

if (sys.nframe() == 0L) {
  derivation_main()
}
