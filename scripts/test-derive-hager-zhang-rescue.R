#!/usr/bin/env Rscript

source(file.path("scripts", "derive-hager-zhang-rescue.R"))

fixture_dir <- tempfile("mize-packet5a-derive-")
dir.create(fixture_dir)
on.exit(unlink(fixture_dir, recursive = TRUE, force = TRUE), add = TRUE)

cases <- c(
  "spd-quadratic-n20",
  "funconstrain-brown_bs",
  "funconstrain-gulf",
  "funconstrain-meyer"
)
summary <- expand.grid(
  case = cases,
  profile = packet5a_profiles(),
  policy = names(packet5a_policies()),
  callback_mode = c("combined", "separate"),
  rep = 1:3,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
summary$dimension <- 2L
summary$initial_f <- 1
summary$final_f <- 0
summary$grad_norm <- 0
summary$nf <- 2L
summary$ng <- 2L
summary$n_fn_call <- 0L
summary$n_gr_call <- 0L
summary$n_fg_call <- 2L
summary$n_hs_call <- 0L
summary$elapsed_sec <- 0.01
summary$iter <- 1L
summary$termination <- "abs_tol=0"
summary$failure_mode <- "ok"
summary$warnings <- ""
summary$error <- ""
summary$objective_reduction_hit <- TRUE
summary$gradient_reduction_hit <- TRUE
summary$reference_gap_target_applicable <- FALSE
summary$reference_gap_hit <- NA
summary$nonstationary_no_step_count <- 0L
summary$false_convergence <- FALSE
summary$objective_worse_than_start <- FALSE
summary$non_descent_escape_count <- 0L
summary$final_metric_nonfinite <- FALSE
summary$callback_or_optimizer_error <- FALSE

progress <- summary[c("case", "profile", "policy", "callback_mode", "rep")]
progress$progress_iter <- 0L
progress$alpha_init <- NA_real_
progress$slope_init <- NA_real_
progress$ls_nf <- NA_integer_
progress$ls_ng <- NA_integer_
progress$ls_reason <- ""
progress$ls_outcome <- ""
progress$g2n <- 0

initializers <- data.frame(
  case = character(),
  profile = character(),
  policy = character(),
  callback_mode = character(),
  rep = integer(),
  iter = integer(),
  unrescued_alpha = numeric(),
  scale_estimate = numeric(),
  trial_evaluation_budget = numeric(),
  required_contractions = numeric(),
  available_contractions = numeric(),
  selected_initial_alpha = numeric(),
  rescue_applied = logical(),
  observer_provider_callbacks = integer(),
  probe_function_evaluations = integer(),
  stringsAsFactors = FALSE
)
packet4 <- data.frame(dummy = TRUE)

write_fixture <- function(name, value) {
  path <- file.path(fixture_dir, paste0(name, ".csv"))
  utils::write.csv(value, path, row.names = FALSE, na = "")
  path
}

config <- list(
  summary = write_fixture("summary", summary),
  progress = write_fixture("progress", progress),
  initializers = write_fixture("initializers", initializers),
  packet4_summary = write_fixture("packet4", packet4)
)

complete <- packet5a_read_inputs(config)
stopifnot(nrow(complete$summary) == 1200L)
stopifnot(identical(complete$scope, "tranche"))
stopifnot(nrow(packet5a_correctness_gates(complete$summary)) == 0L)

packet4_parity <- complete$summary[
  complete$summary$policy == "current-hz-20",
  ,
  drop = FALSE
]
packet4_parity$line_search <- "Hager-Zhang"
packet4_parity$step0 <- "default"
packet4_parity <- packet4_parity[
  rev(seq_len(nrow(packet4_parity))),
  ,
  drop = FALSE
]
rownames(packet4_parity) <- seq.int(1001L, length.out = nrow(packet4_parity))
stopifnot(
  nrow(
    packet5a_packet4_parity_gates(complete$summary, packet4_parity)
  ) ==
    0L
)

selected_progress <- packet5a_select_run_progress(
  progress,
  summary[137L, , drop = FALSE]
)
stopifnot(nrow(selected_progress) == 1L)
stopifnot(selected_progress$case == summary$case[[137L]])
stopifnot(selected_progress$profile == summary$profile[[137L]])
stopifnot(selected_progress$policy == summary$policy[[137L]])
stopifnot(selected_progress$callback_mode == summary$callback_mode[[137L]])
stopifnot(selected_progress$rep == summary$rep[[137L]])

missing_cell <- summary[
  !(summary$case == "funconstrain-brown_bs" &
    summary$profile == "mize-cg-pr+" &
    summary$policy == "rescue-hz-20" &
    summary$callback_mode == "combined"),
  ,
  drop = FALSE
]
config$summary <- write_fixture("summary-missing-cell", missing_cell)
missing_error <- tryCatch(packet5a_read_inputs(config), error = identity)
stopifnot(inherits(missing_error, "error"))
stopifnot(grepl("complete Packet 5A run grid", conditionMessage(missing_error)))

primary_paths <- c(
  runs = write_fixture("runs", summary[1L, , drop = FALSE]),
  progress = write_fixture("progress-primary", progress),
  cases = write_fixture("cases", data.frame(case = cases)),
  summary = write_fixture("summary-primary", summary),
  initializers = write_fixture("initializers-primary", initializers)
)
primary_manifest <- data.frame(
  scope = "tranche",
  artifact = names(primary_paths),
  path = unname(primary_paths),
  data_rows = vapply(
    primary_paths,
    function(path) max(0L, length(readLines(path)) - 1L),
    integer(1)
  ),
  sha256 = vapply(primary_paths, packet5a_sha256, character(1)),
  mize_head = "benchmark-head",
  stringsAsFactors = FALSE
)
config$summary <- primary_paths[["summary"]]
config$progress <- primary_paths[["progress"]]
config$initializers <- primary_paths[["initializers"]]
config$cell_out <- write_fixture("cell-medians", data.frame(value = 1))
config$comparison_out <- write_fixture(
  "policy-comparisons",
  data.frame(value = 2)
)
config$gate_out <- write_fixture("gates", data.frame(value = 3))
config$manifest_out <- write_fixture("manifest", primary_manifest)
config$benchmark_command <- "R_LIBS=/tmp/lib Rscript benchmark.R --scope tranche"
config$derivation_command <- "Rscript derive.R --scope tranche"
final_manifest <- packet5a_finalize_manifest(config)
stopifnot(nrow(final_manifest) == 9L)
stopifnot(
  identical(
    final_manifest$artifact,
    c(
      "runs",
      "progress",
      "cases",
      "summary",
      "initializers",
      "packet4-summary-input",
      "cell-medians",
      "policy-comparisons",
      "gates"
    )
  )
)
stopifnot(all(final_manifest$benchmark_command == config$benchmark_command))
stopifnot(all(final_manifest$derivation_command == config$derivation_command))
stopifnot(all(nzchar(final_manifest$derivation_head)))

message("derive-hager-zhang-rescue completeness checks passed")
