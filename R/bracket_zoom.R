# Shared bracket-and-zoom control for Rasmussen and Schmidt Wolfe searches.

run_bracket_zoom <- function(
  evaluator,
  initial_point,
  initial_alpha,
  condition_policy,
  search_direction,
  method_policy
) {
  if (!isTRUE(is.finite(initial_point$slope)) || initial_point$slope >= 0) {
    return(make_line_search_core_result("non_descent_direction"))
  }

  evaluator_state <- environment(evaluator)
  parameter_projection_witness <- if (
    isTRUE(evaluator_state$parameter_projection_witness_initialized)
  ) {
    evaluator_state$parameter_projection_witness
  } else {
    NULL
  }
  trial_alpha <- initial_alpha
  expansion_state <- list(
    initial_point = initial_point,
    previous_point = initial_point,
    iteration = 0L
  )
  bracket <- NULL

  while (evaluator_state$evaluation_count < evaluator_state$max_evaluations) {
    recovery <- recover_finite_line_point(
      evaluator,
      trial_alpha,
      min_alpha = method_policy$expansion_recovery_lower_bound(expansion_state),
      max_evaluations = Inf
    )
    trial_point <- recovery$line_point
    if (!recovery$succeeded) {
      return(make_line_search_core_result("nonfinite_recovery"))
    }

    expansion <- method_policy$classify_expansion(
      expansion_state,
      trial_point,
      condition_policy
    )
    if (expansion$accepted) {
      return(make_line_search_core_result("wolfe", trial_point))
    }
    if (evaluator_state$evaluation_count >= evaluator_state$max_evaluations) {
      return(make_line_search_core_result("budget_exhausted"))
    }
    if (!is.null(expansion$bracket)) {
      bracket <- expansion$bracket
      break
    }

    proposed_alpha <- method_policy$propose_expansion(
      expansion_state,
      trial_point
    )
    trial_alpha <- safeguard_expansion_alpha(
      proposed_alpha,
      trial_point$alpha
    )
    if (is.null(trial_alpha)) {
      return(make_line_search_core_result("progress_failure"))
    }
    expansion_state$previous_point <- trial_point
    expansion_state$iteration <- expansion_state$iteration + 1L
  }

  if (is.null(bracket)) {
    return(make_line_search_core_result("budget_exhausted"))
  }

  zoom_state <- method_policy$initialize_zoom(bracket)
  direction_scale <- if (is.null(search_direction)) {
    NULL
  } else {
    max(abs(search_direction))
  }

  while (evaluator_state$evaluation_count < evaluator_state$max_evaluations) {
    proposal <- method_policy$propose_zoom(zoom_state, initial_point)
    zoom_state <- proposal$state
    zoom_bounds <- method_policy$zoom_bounds(zoom_state)
    proposal$alpha <- safeguard_zoom_alpha(
      proposal$alpha,
      zoom_bounds[[1L]],
      zoom_bounds[[2L]]
    )
    if (is.null(proposal$alpha)) {
      return(make_line_search_core_result("progress_failure"))
    }
    zoom_endpoints <- method_policy$zoom_endpoints(zoom_state)
    recovery_lower_bound <- method_policy$zoom_recovery_lower_bound(zoom_state)
    proposal_is_definitely_novel <- FALSE
    if (
      !is.null(parameter_projection_witness) &&
        parameter_projection_witness$index <=
          length(zoom_endpoints[[1L]]$parameters) &&
        parameter_projection_witness$index <=
          length(zoom_endpoints[[2L]]$parameters)
    ) {
      witness_index <- parameter_projection_witness$index
      projected_witness <- parameter_projection_witness$initial_parameter +
        proposal$alpha * parameter_projection_witness$direction
      proposal_is_definitely_novel <- is.finite(projected_witness) &&
        projected_witness != zoom_endpoints[[1L]]$parameters[[witness_index]] &&
        projected_witness != zoom_endpoints[[2L]]$parameters[[witness_index]]
    }
    if (proposal_is_definitely_novel) {
      recovery <- recover_finite_line_point(
        evaluator,
        proposal$alpha,
        min_alpha = proposal$alpha,
        max_evaluations = 1
      )
      if (recovery$succeeded) {
        recovery$accepted <- FALSE
        recovery$termination_reason <- NULL
      } else {
        recovery <- continue_bracketed_line_recovery(
          evaluator = evaluator,
          recovery = recovery,
          failed_alpha = proposal$alpha,
          min_alpha = recovery_lower_bound,
          first_endpoint = zoom_endpoints[[1L]],
          second_endpoint = zoom_endpoints[[2L]],
          initial_point = initial_point,
          condition_policy = condition_policy,
          search_direction = search_direction,
          max_evaluations = Inf
        )
      }
    } else {
      recovery <- recover_bracketed_line_point(
        evaluator = evaluator,
        alpha = proposal$alpha,
        min_alpha = recovery_lower_bound,
        first_endpoint = zoom_endpoints[[1L]],
        second_endpoint = zoom_endpoints[[2L]],
        initial_point = initial_point,
        condition_policy = condition_policy,
        search_direction = search_direction,
        max_evaluations = Inf
      )
    }
    if (isTRUE(recovery$accepted)) {
      return(make_line_search_core_result("wolfe", recovery$line_point))
    }
    if (!is.null(recovery$termination_reason)) {
      return(make_line_search_core_result(recovery$termination_reason))
    }
    trial_point <- recovery$line_point
    if (!recovery$succeeded) {
      return(make_line_search_core_result("nonfinite_recovery"))
    }
    if (condition_policy$wolfe(initial_point, trial_point)) {
      return(make_line_search_core_result("wolfe", trial_point))
    }

    zoom_update <- method_policy$process_zoom_trial(
      zoom_state,
      trial_point,
      initial_point,
      condition_policy,
      direction_scale
    )
    zoom_state <- zoom_update$state
    if (!is.null(zoom_update$termination_reason)) {
      return(make_line_search_core_result(zoom_update$termination_reason))
    }
  }

  make_line_search_core_result("budget_exhausted")
}

safeguard_expansion_alpha <- function(proposed_alpha, current_alpha) {
  if (
    is.numeric(proposed_alpha) &&
      length(proposed_alpha) == 1L &&
      isTRUE(proposed_alpha == Inf)
  ) {
    proposed_alpha <- .Machine$double.xmax
  }
  if (
    !is.numeric(proposed_alpha) ||
      length(proposed_alpha) != 1L ||
      !isTRUE(is.finite(proposed_alpha)) ||
      !isTRUE(proposed_alpha > current_alpha)
  ) {
    return(NULL)
  }
  proposed_alpha
}

safeguard_zoom_alpha <- function(proposed_alpha, first_alpha, second_alpha) {
  lower_alpha <- min(first_alpha, second_alpha)
  upper_alpha <- max(first_alpha, second_alpha)
  is_interior <- is.numeric(proposed_alpha) &&
    length(proposed_alpha) == 1L &&
    isTRUE(is.finite(proposed_alpha)) &&
    isTRUE(proposed_alpha > lower_alpha) &&
    isTRUE(proposed_alpha < upper_alpha)
  if (!is_interior) {
    proposed_alpha <- lower_alpha + (upper_alpha - lower_alpha) / 2
  }
  if (
    !isTRUE(is.finite(proposed_alpha)) ||
      !isTRUE(proposed_alpha > lower_alpha) ||
      !isTRUE(proposed_alpha < upper_alpha)
  ) {
    return(NULL)
  }
  proposed_alpha
}
