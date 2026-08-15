# Repository review: `jlmelville/mize`

## Scope

I reviewed the repository snapshot at commit `d4a3901c152e457981a5e87bb73d1e9741f61176` (“more CG tests,” July 8, 2026), rather than a moving branch.

This was a source-level review of the package code, tests, documentation, and CI. R is unavailable in the execution environment, so I could not run `R CMD check`, the test suite, or numerical reproductions. The highest-priority findings below follow directly from control flow and indexing, but they should still receive focused regression tests before fixes are merged.

No files or code were changed.

# Overall assessment

This is **not a repository that needs a wholesale rewrite**. It already has several attributes of a professional R package:

* The exported API is deliberately small: seven functions.
* Runtime dependencies are minimal.
* CI covers release/devel/oldrel combinations across Linux, Windows, and macOS, with actions pinned to immutable SHAs.
* Recent tests include genuine numerical invariants: exact Newton behavior on SPD quadratics, all CG update variants under exact line search, and line-search condition verification.

The problems are concentrated in three places:

1. Rare but real failure paths in L-BFGS, SR1, TN, and line searches.
2. An underspecified and partly broken stateful-initialization contract.
3. An internal lifecycle system that is substantially more elaborate than the public API requires.

I would freeze algorithm additions temporarily, fix the correctness and accounting contracts, and only then simplify the lifecycle machinery.

# Fix-first findings

## 1. L-BFGS crashes when the first curvature pair is rejected

**Severity: high; likely release-blocking.**

`lbfgs_memory_update()` correctly rejects a pair that fails the curvature safeguard by returning the state unchanged. That state can still have zero stored pairs. The caller then unconditionally invokes `lbfgs_solve()`.

`lbfgs_solve()` assumes nonempty history:

* Its reverse loop is `length(rhos):1`.
* It indexes `rhos[[length(rhos)]]` and `yms[[length(yms)]]`.
* With zero history, the loop sequence is effectively `0:1`, and zero-element `[[` indexing is invalid.

A constant-gradient objective, or any first step for which (y^\top s) is zero, negative, or too small, should reach this path on the second iteration. This also affects the L-BFGS preconditioner used by CG: it updates the memory and immediately closes over `lbfgs_solve()` without checking whether a pair was accepted.  The TN preconditioner has the same assumption.

The solution should be conceptually simple:

* Define zero-history L-BFGS as the identity or user-provided `hi` initial approximation.
* Use `rev(seq_along(rhos))`, not `length(rhos):1`.
* Add separate regressions for ordinary L-BFGS, CG with L-BFGS preconditioning, and TN with L-BFGS preconditioning.

The existing invariant test checks directions **after accepted curvature**, so it does not cover this first-rejected-pair case.

## 2. `mize_init()` does not honor its advertised convergence overrides

**Severity: high for the stateful API.**

`make_mize()` always populates `opt$convergence`, including `max_iter`, `max_fn`, `max_gr`, and `max_fg`.

`mize_init()` accepts those same arguments but only assigns them when the corresponding existing field is `NULL`. In a normal optimizer produced by `make_mize()`, they are not `NULL`, so calls such as:

```r
opt <- make_mize()
opt <- mize_init(opt, par, fg, max_iter = 30)
```

normally ignore `max_iter = 30`, along with the evaluation budgets.

That contradicts the documented role of `mize_init()` and the stateful usage examples. This is not merely an awkward default: it means the apparent API does something different from what it says.

Arguments explicitly supplied to `mize_init()` should override the factory values. The implementation will probably need missing-argument detection rather than defaults such as `Inf`, so that it can distinguish “not supplied” from “explicitly set to `Inf`.”

## 3. Reinitializing a terminated optimizer leaves it terminated

**Severity: high for reusable optimizers.**

`mize_init()` resets the iteration and reruns initialization hooks, but it does not visibly clear:

* `terminate`
* `is_terminated`
* `converged`
* `status`
* `message`
* `ok`
* the previous function-convergence baseline
* the adaptive-restart timestamp

It only clears the ordinary value cache at the end.

On the next `mize_step()`, the old `terminate` value causes execution to break during the first stage.  In practical terms, a previously terminated optimizer cannot be cleanly reused even though the function is called `mize_init()`.

There should be one explicit state-reset contract:

* Clear transient and terminal state.
* Reinitialize all algorithm-specific state.
* Rebuild hooks and handlers from a clean registry.
* Decide separately whether evaluation counts are reset.

Counts appear to persist intentionally across reinitialization, and there is a test for that behavior. That is defensible, but some documentation describes them as counts “since initialization,” which is inconsistent. The cleanest resolution would be either a `reset_counts` argument or a separate `mize_reset()` operation, accompanied by precise lifetime-count documentation.

`mize_step()` should also fail immediately with an informative error when `opt$is_initialized` is false. Its documentation says initialization is required, but the function currently begins manipulating state without checking.

## 4. The full BFGS update is unnecessarily cubic-time

**Severity: high performance issue.**

The inverse-BFGS update is currently evaluated as two dense matrix–matrix multiplications:

```text
(I - ρsyᵀ) H (I - ρysᵀ) + ρssᵀ
```

The implementation constructs both dense correction matrices and multiplies them with `H`.

That makes each BFGS update (O(n^3)). The same update can be expanded into:

[
H^+ =
H
-\rho(s(Hy)^\top + (Hy)s^\top)
+(\rho+\rho^2 y^\top Hy)ss^\top,
]

which requires one matrix-vector multiplication and several outer products: (O(n^2)), matching BFGS’s expected dense-storage complexity.

Optimized BLAS may hide the problem at small dimensions, but it does not change the asymptotic cost. This is probably the clearest avoidable inefficiency in the repository. It is particularly worth fixing because the package describes itself as holding up well on higher-dimensional problems, even though L-BFGS is naturally the preferred method there.

A benchmark should accompany the change, covering a range such as (n=50,100,250,500), and verify numerical equivalence, symmetry, and descent behavior.

## 5. SR1 can reference nonexistent variables on its first iteration

**Severity: high but conditional.**

During an SR1 iteration, `hm`, `sm`, and `ym` are assigned only when an old gradient exists. The eventual non-descent fallback nevertheless calls:

```r
bfgs_update(hm, sm, ym, ...)
```

without checking whether a secant pair exists.

On the first iteration, a user-provided inverse Hessian that gives a zero, nonfinite, or ascent direction can therefore send execution into a block where those local variables have never been assigned.

The correct first-iteration fallback is simply steepest descent. A BFGS fallback only makes sense when an actual ((s,y)) pair exists.

The SR1 initial scaling also performs

[
\gamma = \frac{s^\top y}{y^\top y}
]

without the finite/curvature guard used by BFGS.  That scaling should receive equivalent protection.

## 6. TN can exceed `max_fg`, and its most useful extension is missing

**Severity: medium-high.**

The truncated-Newton inner loop checks only `max_gr` before finite-difference Hessian-vector evaluations. It does not check the combined `max_fg` budget. It can therefore consume gradients beyond the combined limit before control returns to the outer convergence checker.

There is a related accounting edge: `tn_inner_cg()` increments the gradient count after `bd_approx()` even when `bd_approx()` returns early without calling `fg$gr`, for example because its finite-difference step is invalid.

The inner limit is also hardcoded to 40 and not available through the public controls.

The obvious deterministic feature here is an exact Hessian-vector-product callback:

```r
fg$hvp(par, v)
```

with fallbacks in this order:

1. `fg$hvp(par, v)`
2. `fg$hs(par) %*% v`
3. Finite-difference gradient approximation

That would improve accuracy, reduce gradient counts, support operator-based large problems, and make TN substantially more compelling. It is a much higher-value addition than another line-search implementation.

## 7. Some line searches can turn budget exhaustion into an accepted, unverified step

**Severity: medium-high correctness/API issue.**

The common line-search interface returns a step and counts, but not a reliable `accepted` flag or failure reason. As a result, exhausting a budget can be treated much like successfully finding an Armijo or Wolfe step.

Two concrete examples:

### Bold Driver

Bold Driver evaluates an initial candidate before its `num_steps < max_fn` check. The counter tracks only subsequent reductions, so it can perform one more candidate evaluation than the computed allowance. If the budget is exhausted while the candidate still fails to reduce the objective, the last candidate remains the selected step and can be applied.

Its non-descent shortcut also uses the minimum positive step rather than a zero step, and its nonfinite fallback selects the minimum step without evaluating the objective there.

### Schmidt Armijo backtracking

The backtracking loop exits when the evaluation budget is exhausted, even if the current trial still fails Armijo. The current trial is then returned as the result.

The unified result should carry at least:

```text
accepted
reason
step
nfn
ngr
conditions_checked
```

Possible reasons include `armijo`, `wolfe`, `budget`, `nonfinite`, `step_too_small`, and `bracket_failed`. The outer optimizer can then make an explicit policy decision:

* Accept a condition-satisfying step.
* Accept a merely decreasing fallback, but label it.
* Return a zero step.
* Terminate with a line-search failure.

At present, those cases are too easily conflated.

## 8. Evaluation budgets do not have one consistent meaning

**Severity: medium.**

Documentation calls `max_fn`, `max_gr`, and `max_fg` maximum evaluations “allowed.”  In practice, they behave partly as hard caps and partly as thresholds checked after work has already occurred.

Examples include:

* The gradient required for the search direction is evaluated before the line-search budget check.
* Hager–Zhang step initialization can spend a function evaluation before remaining budgets are calculated.
* TN only considers `max_gr`.
* `mize()` always evaluates the returned objective if it is not cached, even after the main loop has exhausted the function budget.
* The outer loop uses a first-iteration heuristic to reserve a final function evaluation by reducing stored limits after observing algorithm behavior.

This needs a documented decision:

**Hard-cap interpretation:** never invoke a callback after a budget reaches zero. The final `f` may then be unavailable.

**Stopping-threshold interpretation:** stop at the first safe boundary after a threshold is reached, and document that reported counts may exceed it by a bounded amount.

Either is workable. The current mixture is difficult for callers who use function evaluations as a meaningful cost unit.

A useful test helper would wrap every callback in an independent counter and compare the package’s `nf` and `ng` against the actual invocations after every operation.

# API and maintainability

## 9. Input validation is much weaker than the numerical code deserves

**Severity: medium, with some high-impact malformed-input cases.**

The factory exposes a very large flat parameter surface, but many checks enforce only one inequality. They generally do not enforce:

* Scalarity
* Finiteness
* Integer-valued counts
* Nonmissing values
* Method-specific relevance
* Valid callback dimensions

Examples include fractional `memory` or iteration counts, zero convergence/logging cadences, invalid preconditioner strings that silently act like no preconditioner, and a constant line search with no numeric `step0`.

The callback boundary is even more important:

* `fn` should return one finite numeric value.
* `gr` should return a finite numeric vector of exactly `length(par)`.
* Combined `fg` should have the required fields and dimensions.
* `hs` and `hi` should return an (n)-vector or (n\times n) matrix.
* Matrix values should be finite; symmetry should be checked within tolerance where required.
* Diagonal Hessians should be checked for usable entries before division.

Without that, helpers such as `dot()` silently use ordinary R vector recycling for mismatched lengths.

The gradient checker already contains substantially better validation helpers: it verifies callback types, finite parameter values, scalar objective output, and correctly sized finite gradients. Those should become the basis of centralized optimizer validation rather than remaining isolated inside the diagnostic function.

For API organization, I would preserve the existing arguments for compatibility but add internal or public control constructors such as:

* convergence controls
* line-search controls
* momentum controls
* TN controls

That would make method-specific validation and documentation much more manageable.

## 10. The lifecycle/hook system is the clear overengineering hotspot

**Severity: medium maintainability risk.**

The package’s core iteration is represented through:

* Dynamically registered hooks
* Function attributes carrying string events
* Advice types and join points
* Stage and substage wrappers
* Dependencies resolved by constructing names such as `require_*` and calling `get0()`

The source itself warns, “You don't want to look too closely at any of this.”

A concrete sign of fragility is that the comments say a one-token event such as `"init"` is valid and should imply `"during"`, but the parser immediately assumes two tokens and reads `event_tok[2]`.  Dependencies are also resolved dynamically from the package namespace rather than represented as explicit functions or typed objects.

This machinery is not exposed as a public extension API. Consequently, the package pays much of the complexity cost without giving users the corresponding flexibility.

I would retain the valuable stateful API but gradually replace the internals with an explicit pipeline:

1. Ensure current objective/gradient requirements.
2. Calculate direction.
3. Search for or calculate a step.
4. Validate the candidate.
5. Commit or reject state.
6. Update method memory.
7. Notify optional observers.

Method objects can still supply `initialize`, `direction`, and `after_step` functions. Momentum can still be a sequence of explicit stages. What is unnecessary is encoding control flow in attributed strings and dynamically discovered function names.

This should happen only after the current behavior is protected by state-transition tests. A rewrite of this subsystem without such tests would be risky.

## 11. There are several dormant or stringly typed controls

These are not the top priorities, but they are useful cleanup targets:

* `cg_direction(ortho_check = TRUE)` calls `sub_stage$cg_restart(...)`, but no `cg_restart` member is stored in the substage. The default `FALSE` hides the problem.
* Arbitrary `preconditioner` strings are accepted; anything other than exact `"l-bfgs"` silently behaves as no preconditioner.
* Arbitrary character `tn_init` values fall through to the “previous direction” behavior rather than being validated.
* `ras_ls()` has a self-referential default, `verbose = verbose`, which should fail when the internal helper is called without explicitly supplying `verbose`.
* PHESS is selectable through the main method argument but is described as experimental; its useful refresh interval is not exposed. It should either become fully supported or be moved out of the ordinary stable method list.

# Convergence behavior and documentation

## 12. Convergence precedence and comparisons need an explicit policy

**Severity: medium.**

`check_mize_convergence()` checks evaluation budgets before step, gradient, and function tolerances. Therefore, a point that satisfies a convergence tolerance on exactly the evaluation that reaches a budget is reported as `budget_exhausted`, not `converged`.

That might be intentional, but it should be documented and tested. Many users would reasonably prefer successful convergence to take precedence when both happen simultaneously.

Other details:

* `max_iter` is tested with `==`, not `>=`, which is fragile in a manually driven stateful loop.
* Absolute and relative function tolerances use strict `<`; gradient tolerances use `<=`.
* Relative function change divides by `min(abs(f_old), abs(f_new))`. That is very strict and unstable around zero; it should at least be documented, and a scale such as `max(1, abs(f_old), abs(f_new))` would usually be easier to reason about.
* Tolerances should be validated as finite, nonnegative scalar numbers.

## 13. Some documentation currently describes different behavior from the code

These are professional-polish issues, but several are substantive:

* The main documentation says `abs_tol` is reached when the **absolute function value** falls below the threshold. The implementation tests the **absolute difference between consecutive function values**.
* `max_fg` is described in one place as a maximum number of gradient evaluations, although the code treats it as combined function-plus-gradient evaluations.
* The development NEWS attributes finite-difference gradient checking to `check_mize_convergence()` rather than `check_mize_gradient()`, and contains duplicated wording in the status-field entry.
* The stateful documentation and implementation disagree over whether counts are per initialization or over the lifetime of the optimizer.
* The `nest_q` descriptions are inconsistent: the user-facing text describes zero as strongly convex, while the implementation comments correctly associate zero with the largest momentum and one with zero momentum.
* Some older stateful-vignette language still refers to manually supplying an iteration number and to fields that are no longer part of the public state.

A documentation sweep should follow the state and convergence contract decisions, not precede them.

# Smaller inefficiencies and code-quality issues

These are suitable for a later cleanup PR:

* `update_progress()` repeatedly `rbind()`s a growing data frame, giving quadratic copying behavior for long traces. Accumulate rows as lists or preallocate, then bind once.
* `safe_chol()` accepts `eps` but hardcodes `1e-10`; it replaces only negative eigenvalues, leaving zero eigenvalues unchanged, and uses a general eigendecomposition rather than `symmetric = TRUE`.
* `normalize()` calculates the same norm twice.
* `opt_clear_cache()` retains old values and writes the string `"invalid"` into iteration fields. Clearing both value and tag to `NULL` would be more type-stable.
* `is_finite_numeric()` does not enforce scalarity or all-elements finiteness.
* `is_in_range()` has confusing “open”/“closed” terminology and uses function-valued `ifelse()` where direct comparisons would be clearer.
* There are parallel generations of line-search code with somewhat different counters, status conventions, and failure handling. Removing unused or superseded helpers after coverage analysis would lower the maintenance burden.

# Missing features worth adding

Within the deterministic, non–trust-region scope you specified, I would prioritize these:

## 1. Exact Hessian-vector products

This is the most obvious missing numerical capability. It directly upgrades TN and enables large structured problems without materializing a dense Hessian.

## 2. A general preconditioner callback

Rather than accepting only the string `"l-bfgs"`, allow an operator such as:

```r
preconditioner <- function(par, residual, state) ...
```

or a simpler `fg$precondition(par, residual)`. Validate the output and fall back safely when it is nonfinite or not positive in the residual inner product.

## 3. Parameter scaling

An equivalent of `optim()`’s `parscale`, or a diagonal scaling/operator abstraction, would benefit SD, CG, BFGS, TN, convergence tolerances, and finite-difference HVPs. Poorly scaled parameters are a more common practical problem than the absence of another optimizer.

## 4. Damped BFGS

The existing code comments acknowledge that standard references prefer damping to simply skipping bad curvature pairs.  Powell-style damping would be a useful robustness option for both full BFGS and possibly the L-BFGS memory update.

## 5. Observer and user-stop callbacks

The stateful API provides this manually, but standalone `mize()` would benefit from a callback receiving iteration, parameters, objective, gradient norms, counts, step status, and optimizer state. It could request continuation, clean termination, or checkpointing.

## 6. Hessian and HVP diagnostics

`check_mize_gradient()` is a good start. Natural companions would check:

* Hessian symmetry
* Hessian-vector consistency against finite differences
* `hi` consistency with `hs`
* Dimensions and finite values
* Directional second derivatives

## 7. Spectral gradient and nonmonotone Armijo, later

Barzilai–Borwein step initialization with a nonmonotone Armijo search would be a reasonable lightweight deterministic addition. I would not add it until the current line-search acceptance/status contract has been cleaned up.

I do **not** consider box constraints an “obvious omission” from this repository because it explicitly presents itself as an unconstrained optimizer. L-BFGS-B would be useful, but it changes package scope rather than simply filling a hole. I also would not add more legacy line-search ports before consolidating the ones already present.

# Tests I would add before architectural work

The most valuable next test group is small and targeted:

1. **Rejected first L-BFGS pair:** ordinary L-BFGS, CG preconditioner, and TN preconditioner, all beginning with empty history.
2. **SR1 first-iteration fallback:** zero, indefinite, malformed, and nonfinite user `hi` results.
3. **Stateful reinitialization:** after convergence, budget exhaustion, nonfinite termination, and adaptive restart; verify terminal state and function-change baselines are clean.
4. **Initialization override:** verify every convergence argument supplied to `mize_init()` actually wins.
5. **Exact budget accounting:** independently instrument `fn`, `gr`, and combined `fg`; test zero and one-evaluation boundaries, Hager–Zhang initialization, TN inner iterations, logging, summaries, and final return evaluation.
6. **Line-search failure status:** deliberately exhaust each search before its condition is satisfied and verify that an unaccepted step cannot masquerade as success.
7. **Malformed callback contracts:** wrong gradient length, vector-valued objective, malformed combined result, wrong Hessian dimensions, NA/NaN/Inf, and nonsymmetric matrices.
8. **Configuration fuzzing:** NA, NaN, Inf, length-two values, fractional iteration counts, zero cadence, unknown strings, and missing schedule values.
9. **Randomized numerical invariants:** random SPD quadratics, transformed/scaled versions of the same problem, and equality of L-BFGS and full-BFGS directions in low dimensions while sufficient history is retained.
10. **Performance regression:** benchmark the current BFGS update against its quadratic-time expansion and measure lifecycle overhead on a very cheap quadratic objective.

# Recommended cleanup order

1. **Correctness regressions first:** empty L-BFGS history, SR1 first step, stateful reinitialization, TN budget handling, and line-search budget failure.
2. **Replace the BFGS update:** retain numerical-equivalence tests and add a benchmark.
3. **Define contracts:** hard versus soft budgets, line-search acceptance, reinitialization, count lifetime, and convergence precedence.
4. **Centralize validation:** both configuration and callback outputs.
5. **Repair documentation:** regenerate reference material and update vignettes and NEWS from the finalized contracts.
6. **Simplify the lifecycle internals:** preserve the public stateful API, but move toward an explicit iteration pipeline.
7. **Add HVPs, general preconditioning, scaling, and damping.**

The key judgment is that `mize` already has enough optimization algorithms. Its next increase in quality will come from making failure, initialization, accounting, and callback behavior boringly predictable—not from adding another method.
