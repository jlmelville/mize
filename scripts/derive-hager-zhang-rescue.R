#!/usr/bin/env Rscript

source(file.path("scripts", "derive-optimizer-benchmark.R"))
source(file.path("scripts", "benchmark-hager-zhang-rescue.R"))

packet5a_derivation_usage <- function(status = 0L) {
  cat(
    paste(
      "Usage:",
      "  Rscript scripts/derive-hager-zhang-rescue.R [options]",
      "",
      "Options:",
      "  --summary PATH            Persisted Packet 5A summary CSV.",
      "  --progress PATH           Persisted Packet 5A progress CSV.",
      "  --initializers PATH       Persisted initializer evidence CSV.",
      "  --packet4-summary PATH    Accepted Packet 4 summary CSV.",
      "  --cell-out PATH           Deterministic cell-median CSV.",
      "  --comparison-out PATH     Candidate/primary comparison CSV.",
      "  --gate-out PATH           Frozen-gate CSV.",
      "  --help                    Show this help.",
      sep = "\n"
    )
  )
  quit(save = "no", status = status)
}

packet5a_parse_derivation_args <- function(args) {
  option_names <- c(
    `--summary` = "summary",
    `--progress` = "progress",
    `--initializers` = "initializers",
    `--packet4-summary` = "packet4_summary",
    `--cell-out` = "cell_out",
    `--comparison-out` = "comparison_out",
    `--gate-out` = "gate_out"
  )
  config <- as.list(rep(NA_character_, length(option_names)))
  names(config) <- unname(option_names)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      packet5a_derivation_usage()
    }
    if (!arg %in% names(option_names)) {
      stop("Unknown option: ", arg, call. = FALSE)
    }
    if (i == length(args)) {
      stop(arg, " requires a value", call. = FALSE)
    }
    config[[option_names[[arg]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  missing <- names(config)[is.na(unlist(config)) | !nzchar(unlist(config))]
  if (length(missing) > 0L) {
    stop(
      "Missing required option(s): ",
      paste(gsub("_", "-", paste0("--", missing)), collapse = ", "),
      call. = FALSE
    )
  }
  config
}

packet5a_required_columns <- function(value, required, label) {
  missing <- setdiff(required, names(value))
  if (length(missing) > 0L) {
    stop(
      label,
      " lacks required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

packet5a_run_key <- function(value, include_rep = TRUE) {
  columns <- c("case", "profile", "policy", "callback_mode")
  if (include_rep) {
    columns <- c(columns, "rep")
  }
  do.call(paste, c(value[columns], sep = "\r"))
}

packet5a_read_inputs <- function(config) {
  summary <- utils::read.csv(
    config$summary,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  progress <- utils::read.csv(
    config$progress,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  initializers <- utils::read.csv(
    config$initializers,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  packet4 <- utils::read.csv(
    config$packet4_summary,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  packet5a_required_columns(
    summary,
    c(
      "case",
      "profile",
      "policy",
      "callback_mode",
      "rep",
      "dimension",
      "initial_f",
      "final_f",
      "grad_norm",
      "nf",
      "ng",
      "n_fn_call",
      "n_gr_call",
      "n_fg_call",
      "n_hs_call",
      "elapsed_sec",
      "iter",
      "termination",
      "failure_mode",
      "warnings",
      "error",
      "objective_reduction_hit",
      "gradient_reduction_hit",
      "reference_gap_target_applicable",
      "reference_gap_hit",
      "nonstationary_no_step_count",
      "false_convergence",
      "objective_worse_than_start",
      "non_descent_escape_count",
      "final_metric_nonfinite",
      "callback_or_optimizer_error"
    ),
    "Packet 5A summary"
  )
  packet5a_required_columns(
    progress,
    c(
      "case",
      "profile",
      "policy",
      "callback_mode",
      "rep",
      "progress_iter",
      "alpha_init",
      "slope_init",
      "ls_nf",
      "ls_ng",
      "ls_reason",
      "ls_outcome",
      "g2n"
    ),
    "Packet 5A progress"
  )
  packet5a_required_columns(
    initializers,
    c(
      "case",
      "profile",
      "policy",
      "callback_mode",
      "rep",
      "iter",
      "unrescued_alpha",
      "scale_estimate",
      "trial_evaluation_budget",
      "required_contractions",
      "available_contractions",
      "selected_initial_alpha",
      "rescue_applied",
      "observer_provider_callbacks",
      "probe_function_evaluations"
    ),
    "Packet 5A initializer evidence"
  )

  expected_profiles <- packet5a_profiles()
  expected_policies <- names(packet5a_policies())
  expected_callbacks <- c("combined", "separate")
  if (!setequal(summary$profile, expected_profiles)) {
    stop("summary does not contain the frozen Packet 5A profiles")
  }
  if (!setequal(summary$policy, expected_policies)) {
    stop("summary does not contain the frozen Packet 5A policies")
  }
  if (!setequal(summary$callback_mode, expected_callbacks)) {
    stop("summary does not contain both Packet 5A callback modes")
  }
  expected_case_count <- if (
    setequal(
      summary$case,
      c(
        "spd-quadratic-n20",
        "funconstrain-brown_bs",
        "funconstrain-gulf",
        "funconstrain-meyer"
      )
    )
  ) {
    4L
  } else {
    36L
  }
  if (length(unique(summary$case)) != expected_case_count) {
    stop("summary case set is neither the tranche nor all-35 scope")
  }
  expected <- expand.grid(
    case = sort(unique(summary$case)),
    profile = expected_profiles,
    policy = expected_policies,
    callback_mode = expected_callbacks,
    rep = 1:3,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  summary_key <- packet5a_run_key(summary)
  expected_key <- packet5a_run_key(expected)
  if (anyDuplicated(summary_key) || !setequal(summary_key, expected_key)) {
    stop("summary does not contain the complete Packet 5A run grid")
  }
  progress_key <- unique(packet5a_run_key(progress))
  if (!setequal(progress_key, expected_key)) {
    stop("progress does not cover the complete Packet 5A run grid")
  }
  initializer_key <- packet5a_run_key(initializers)
  summary_initializer_keys <- summary_key[
    summary$policy %in% c("current-hz-20", "rescue-hz-20", "current-hz-100")
  ]
  if (
    any(!initializer_key %in% summary_initializer_keys) ||
      anyDuplicated(paste(initializer_key, initializers$iter, sep = "\r"))
  ) {
    stop("initializer evidence has invalid or duplicate run/iteration keys")
  }
  if (any(initializers$observer_provider_callbacks != 0L)) {
    stop("initializer observer records provider callbacks")
  }
  list(
    summary = summary,
    progress = progress,
    initializers = initializers,
    packet4 = packet4,
    scope = if (expected_case_count == 4L) "tranche" else "all35"
  )
}

packet5a_constant_character <- function(values, label) {
  values <- unique(as.character(values))
  if (length(values) != 1L) {
    stop(label, " differs across repetitions", call. = FALSE)
  }
  values[[1L]]
}

packet5a_target_signature <- function(rows) {
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
    "reference applicability"
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
  list(
    objective_hit = objective_hit,
    gradient_hit = gradient_hit,
    reference_applicable = reference_applicable,
    reference_hit = reference_hit,
    signature = paste(
      benchmark_logical_code(objective_hit),
      benchmark_logical_code(gradient_hit),
      reference_code,
      sep = "/"
    )
  )
}

packet5a_cell_medians <- function(inputs) {
  summary <- inputs$summary
  grouping <- c("case", "profile", "policy", "callback_mode")
  ordering <- do.call(order, c(summary[grouping], list(summary$rep)))
  summary <- summary[ordering, , drop = FALSE]
  keys <- do.call(paste, c(summary[grouping], sep = "\r"))
  groups <- split(seq_len(nrow(summary)), keys)
  rows <- lapply(groups, function(indices) {
    group <- summary[indices, , drop = FALSE]
    targets <- packet5a_target_signature(group)
    initializer <- inputs$initializers[
      inputs$initializers$case == group$case[[1L]] &
        inputs$initializers$profile == group$profile[[1L]] &
        inputs$initializers$policy == group$policy[[1L]] &
        inputs$initializers$callback_mode == group$callback_mode[[1L]],
      ,
      drop = FALSE
    ]
    trigger_counts <- if (nrow(initializer) == 0L) {
      rep(0L, 3L)
    } else {
      vapply(
        1:3,
        function(rep) {
          sum(initializer$rep == rep & initializer$rescue_applied %in% TRUE)
        },
        integer(1)
      )
    }
    if (length(unique(trigger_counts)) != 1L) {
      stop("rescue trigger count differs across repetitions")
    }
    data.frame(
      case = group$case[[1L]],
      profile = group$profile[[1L]],
      policy = group$policy[[1L]],
      callback_mode = group$callback_mode[[1L]],
      dimension = group$dimension[[1L]],
      initial_f = group$initial_f[[1L]],
      termination = packet5a_constant_character(
        group$termination,
        "termination"
      ),
      failure_mode = packet5a_constant_character(
        group$failure_mode,
        "failure mode"
      ),
      objective_hit = targets$objective_hit,
      gradient_hit = targets$gradient_hit,
      reference_applicable = targets$reference_applicable,
      reference_hit = targets$reference_hit,
      target_signature = targets$signature,
      final_f_median = stats::median(group$final_f),
      final_f_min = min(group$final_f),
      final_f_max = max(group$final_f),
      grad_norm_median = stats::median(group$grad_norm),
      grad_norm_min = min(group$grad_norm),
      grad_norm_max = max(group$grad_norm),
      iter_median = stats::median(group$iter),
      n_fn_call_median = stats::median(group$n_fn_call),
      n_gr_call_median = stats::median(group$n_gr_call),
      n_fg_call_median = stats::median(group$n_fg_call),
      n_hs_call_median = stats::median(group$n_hs_call),
      elapsed_sec_median = stats::median(group$elapsed_sec),
      elapsed_sec_min = min(group$elapsed_sec),
      elapsed_sec_max = max(group$elapsed_sec),
      nonstationary_no_step = any(group$nonstationary_no_step_count > 0L),
      nonstationary_no_step_count = group$nonstationary_no_step_count[[1L]],
      rescue_trigger_count = trigger_counts[[1L]],
      stringsAsFactors = FALSE
    )
  })
  cells <- do.call(rbind, rows)
  rownames(cells) <- NULL
  cells[do.call(order, cells[grouping]), , drop = FALSE]
}

packet5a_policy_comparisons <- function(cells) {
  primary <- cells[cells$policy == "current-hz-20", , drop = FALSE]
  candidate <- cells[cells$policy == "rescue-hz-20", , drop = FALSE]
  key_columns <- c("case", "profile", "callback_mode")
  joined <- merge(
    primary,
    candidate,
    by = key_columns,
    suffixes = c("_current", "_rescue"),
    sort = FALSE
  )
  target_loss <- function(current, rescue) isTRUE(current) && !isTRUE(rescue)
  joined$objective_target_lost <- mapply(
    target_loss,
    joined$objective_hit_current,
    joined$objective_hit_rescue
  )
  joined$gradient_target_lost <- mapply(
    target_loss,
    joined$gradient_hit_current,
    joined$gradient_hit_rescue
  )
  joined$reference_target_lost <- mapply(
    target_loss,
    joined$reference_hit_current,
    joined$reference_hit_rescue
  )
  joined$new_nonstationary_no_step <-
    !joined$nonstationary_no_step_current &
    joined$nonstationary_no_step_rescue
  joined$objective_allowance <- 1e-8 *
    pmax(
      1,
      abs(joined$initial_f_current),
      abs(joined$final_f_median_current)
    )
  joined$objective_materially_worse <-
    joined$final_f_median_rescue >
      joined$final_f_median_current + joined$objective_allowance
  joined$final_f_difference <-
    joined$final_f_median_rescue - joined$final_f_median_current
  joined$grad_norm_difference <-
    joined$grad_norm_median_rescue - joined$grad_norm_median_current
  for (resource in c(
    "n_fn_call_median",
    "n_gr_call_median",
    "n_fg_call_median",
    "n_hs_call_median",
    "elapsed_sec_median"
  )) {
    current <- joined[[paste0(resource, "_current")]]
    rescue <- joined[[paste0(resource, "_rescue")]]
    joined[[paste0(resource, "_difference")]] <- rescue - current
    joined[[paste0(resource, "_ratio")]] <- ifelse(
      current > 0,
      rescue / current,
      NA_real_
    )
  }
  joined <- joined[do.call(order, joined[key_columns]), , drop = FALSE]
  rownames(joined) <- NULL
  joined
}

packet5a_empty_gates <- function() {
  data.frame(
    gate = character(),
    case = character(),
    profile = character(),
    policy = character(),
    callback_mode = character(),
    rep = integer(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}

packet5a_gate_row <- function(
  gate,
  case = "",
  profile = "",
  policy = "",
  callback_mode = "",
  rep = NA_integer_,
  detail = ""
) {
  data.frame(
    gate = gate,
    case = case,
    profile = profile,
    policy = policy,
    callback_mode = callback_mode,
    rep = rep,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

packet5a_correctness_gates <- function(summary) {
  gates <- list()
  add <- function(indices, gate, detail) {
    for (index in indices) {
      row <- summary[index, , drop = FALSE]
      gates[[length(gates) + 1L]] <<- packet5a_gate_row(
        gate,
        row$case,
        row$profile,
        row$policy,
        row$callback_mode,
        row$rep,
        detail
      )
    }
  }
  add(
    which(
      summary$callback_or_optimizer_error |
        summary$final_metric_nonfinite |
        (!is.na(summary$error) & nzchar(summary$error)) |
        (!is.na(summary$warnings) & nzchar(summary$warnings))
    ),
    "callback_or_numerical_failure",
    "error, warning, or nonfinite result"
  )
  add(which(summary$false_convergence), "false_convergence", "")
  add(which(summary$objective_worse_than_start), "objective_worse", "")
  add(
    which(summary$nonstationary_no_step_count >= 2L),
    "repeated_nonstationary_no_step",
    ""
  )
  add(
    which(summary$non_descent_escape_count > 0L),
    "non_descent_escape",
    ""
  )
  if (length(gates) == 0L) packet5a_empty_gates() else do.call(rbind, gates)
}

packet5a_packet4_parity_gates <- function(summary, packet4) {
  current <- summary[summary$policy == "current-hz-20", , drop = FALSE]
  accepted <- packet4[
    packet4$line_search == "Hager-Zhang" &
      packet4$case %in% current$case,
    ,
    drop = FALSE
  ]
  key <- c("case", "profile", "callback_mode", "rep")
  current <- current[do.call(order, current[key]), , drop = FALSE]
  accepted <- accepted[do.call(order, accepted[key]), , drop = FALSE]
  current_key <- do.call(paste, c(current[key], sep = "\r"))
  accepted_key <- do.call(paste, c(accepted[key], sep = "\r"))
  if (!setequal(current_key, accepted_key)) {
    return(packet5a_gate_row(
      "packet4_baseline_parity",
      detail = "accepted Packet 4 keys do not match"
    ))
  }
  accepted <- accepted[match(current_key, accepted_key), , drop = FALSE]
  rownames(current) <- NULL
  rownames(accepted) <- NULL
  common <- setdiff(
    intersect(names(current), names(accepted)),
    c(
      "elapsed_sec",
      "step0",
      "line_search",
      "later_initializer",
      "ls_max_fn",
      "scale_rescue",
      "policy"
    )
  )
  mismatch <- vapply(
    seq_len(nrow(current)),
    function(index) {
      !isTRUE(all.equal(
        current[index, common, drop = FALSE],
        accepted[index, common, drop = FALSE],
        tolerance = 0,
        check.attributes = TRUE
      ))
    },
    logical(1)
  )
  if (!any(mismatch)) {
    return(packet5a_empty_gates())
  }
  rows <- current[mismatch, , drop = FALSE]
  do.call(
    rbind,
    lapply(seq_len(nrow(rows)), function(index) {
      packet5a_gate_row(
        "packet4_baseline_parity",
        rows$case[[index]],
        rows$profile[[index]],
        rows$policy[[index]],
        rows$callback_mode[[index]],
        rows$rep[[index]],
        "non-timing fields differ"
      )
    })
  )
}

packet5a_comparison_gates <- function(comparisons, cells, scope) {
  gates <- list()
  add_comparison <- function(rows, gate, detail) {
    for (index in which(rows)) {
      row <- comparisons[index, , drop = FALSE]
      gates[[length(gates) + 1L]] <<- packet5a_gate_row(
        gate,
        row$case,
        row$profile,
        "rescue-hz-20",
        row$callback_mode,
        detail = detail
      )
    }
  }
  target_loss <- comparisons$objective_target_lost |
    comparisons$gradient_target_lost |
    comparisons$reference_target_lost
  if (scope == "all35") {
    add_comparison(
      target_loss,
      "target_loss",
      "candidate loses a primary target"
    )
    add_comparison(
      comparisons$new_nonstationary_no_step,
      "new_nonstationary_no_step",
      "candidate creates a no-step configuration"
    )
    add_comparison(
      comparisons$objective_materially_worse,
      "objective_materially_worse",
      "candidate exceeds the frozen objective allowance"
    )
  }

  tn_controls <- comparisons$profile == "mize-tn" &
    comparisons$case %in% c("funconstrain-brown_bs", "funconstrain-gulf")
  add_comparison(
    tn_controls &
      (target_loss |
        comparisons$new_nonstationary_no_step |
        comparisons$objective_materially_worse),
    "tn_regression_control",
    "Brown/Gulf TN hard gate"
  )

  spd <- comparisons$case == "spd-quadratic-n20"
  spd_fields <- c(
    "termination",
    "failure_mode",
    "target_signature",
    "final_f_median",
    "grad_norm_median",
    "iter_median",
    "n_fn_call_median",
    "n_gr_call_median",
    "n_fg_call_median",
    "n_hs_call_median",
    "nonstationary_no_step_count"
  )
  spd_mismatch <- rep(FALSE, nrow(comparisons))
  for (field in spd_fields) {
    spd_mismatch <- spd_mismatch |
      comparisons[[paste0(field, "_current")]] !=
        comparisons[[paste0(field, "_rescue")]]
  }
  spd_mismatch[is.na(spd_mismatch)] <- TRUE
  spd_mismatch <- spd_mismatch |
    comparisons$rescue_trigger_count_rescue != 0L
  add_comparison(
    spd & spd_mismatch,
    "spd_candidate_parity",
    "candidate changes the SPD control"
  )

  witness <- comparisons$case %in%
    c("funconstrain-brown_bs", "funconstrain-meyer") &
    comparisons$profile %in% c("mize-cg-pr+", "mize-cg-hz+")
  witness_failure <- witness &
    (comparisons$rescue_trigger_count_rescue < 1L |
      comparisons$nonstationary_no_step_rescue |
      comparisons$target_signature_current !=
        comparisons$target_signature_rescue)
  add_comparison(
    witness_failure,
    "primary_witness_not_rescued",
    "rescue must fire, remove no-step, and preserve targets"
  )

  witness_keys <- comparisons[witness, c("case", "profile", "callback_mode")]
  for (index in seq_len(nrow(witness_keys))) {
    key <- witness_keys[index, , drop = FALSE]
    alternatives <- cells[
      cells$case == key$case &
        cells$profile == key$profile &
        cells$callback_mode == key$callback_mode &
        cells$policy %in% c("slope-ratio-20", "current-hz-100"),
      ,
      drop = FALSE
    ]
    if (nrow(alternatives) != 2L || all(alternatives$nonstationary_no_step)) {
      gates[[length(gates) + 1L]] <- packet5a_gate_row(
        "mechanism_comparator_failure",
        key$case,
        key$profile,
        callback_mode = key$callback_mode,
        detail = "neither scale nor larger-budget comparator removes no-step"
      )
    }
  }

  same_target <- comparisons$target_signature_current ==
    comparisons$target_signature_rescue
  for (resource in c(
    "n_fn_call_median",
    "n_gr_call_median",
    "n_fg_call_median",
    "n_hs_call_median",
    "elapsed_sec_median"
  )) {
    ratio <- comparisons[[paste0(resource, "_ratio")]]
    difference <- comparisons[[paste0(resource, "_difference")]]
    minimum <- if (resource == "elapsed_sec_median") 0.05 else 50
    crossed <- same_target &
      is.finite(ratio) &
      ratio > 10 &
      difference >= minimum
    add_comparison(
      crossed,
      "comparable_quality_resource_cost",
      paste0(resource, " crosses the frozen resource gate")
    )
  }
  if (length(gates) == 0L) packet5a_empty_gates() else do.call(rbind, gates)
}

packet5a_derive <- function(config) {
  benchmark_preflight_output_paths(c(
    `--summary` = config$summary,
    `--progress` = config$progress,
    `--initializers` = config$initializers,
    `--packet4-summary` = config$packet4_summary,
    `--cell-out` = config$cell_out,
    `--comparison-out` = config$comparison_out,
    `--gate-out` = config$gate_out
  ))
  inputs <- packet5a_read_inputs(config)
  cells <- packet5a_cell_medians(inputs)
  comparisons <- packet5a_policy_comparisons(cells)
  gates <- packet5a_rbind(list(
    packet5a_correctness_gates(inputs$summary),
    packet5a_packet4_parity_gates(inputs$summary, inputs$packet4),
    packet5a_comparison_gates(comparisons, cells, inputs$scope)
  ))
  if (nrow(gates) == 0L) {
    gates <- packet5a_empty_gates()
  } else {
    gates <- unique(gates)
    gates <- gates[
      order(
        gates$gate,
        gates$case,
        gates$profile,
        gates$policy,
        gates$callback_mode,
        gates$rep,
        gates$detail,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    rownames(gates) <- NULL
  }
  write_decision_artifact(cells, config$cell_out, "Packet 5A cell medians")
  write_decision_artifact(
    comparisons,
    config$comparison_out,
    "Packet 5A policy comparisons"
  )
  write_decision_artifact(gates, config$gate_out, "Packet 5A gates")
  invisible(list(cells = cells, comparisons = comparisons, gates = gates))
}

packet5a_derivation_main <- function() {
  config <- packet5a_parse_derivation_args(commandArgs(trailingOnly = TRUE))
  packet5a_derive(config)
}

if (sys.nframe() == 0L) {
  packet5a_derivation_main()
}
