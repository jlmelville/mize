# Developer Scripts

This directory contains developer utilities that are useful when working from
the source repository. These scripts are not part of the package API and are
excluded from R package builds.

## Optimizer Benchmark Harness

`benchmark-optimizers.R` compares `mize` optimizers with `stats::optim()` on a
small set of smooth unconstrained optimization problems. It is a permanent
developer benchmark harness, not a pass/fail package validation command.
Optional `funconstrain` case discovery and starting-point resolution live in
the sourceable `funconstrain-mize-harness.R` helper; sourcing that helper only
defines its adapter functions and does not run a benchmark.

Use `testthat::test_local()` or `R CMD check` to check that the package works
as expected. Use this benchmark script when you want local evidence about
optimizer behavior, evaluation counts, wall time, termination modes, or
sensitivity to `mize` line-search and `step0` settings.

By default, mize receives each case's combined objective/gradient callback.
Use `--callbacks separate` or `--callbacks combined,separate` to measure the
separate-callback path under otherwise identical settings.

### Method Selection

Use `--methods` to select one or more explicit optimizer profiles:

- `stats-bfgs`, `stats-cg`, and `stats-l-bfgs-b`
- `mize-l-bfgs`, `mize-bfgs`, `mize-cg-pr+`, `mize-cg-hz+`, and `mize-tn`
- `mize-newton` for mize's exact-Hessian Newton method

The default is the original three `stats::optim()` and five mize profiles, so
existing commands retain their grids. Exact Newton is opt-in because it
requires a case with an `fg$hs` callback. Each selected mize profile is crossed
with the requested line searches, first-step settings, and callback modes;
`stats::optim()` profiles keep their ordinary separate-callback execution.

### Funconstrain Adapter Schema

The sourceable helper preserves the benchmark's `case$name`, `case$source`,
`case$par`, and `case$fg` fields. It maps funconstrain callbacks as follows:

- `fn` to `fg$fn`
- `gr` to `fg$gr`
- optional `fg` to `fg$fg`
- `he` to mize's `fg$hs` name

Additional developer metadata is kept outside `fg`:

- `start` records the resolved parameters, requested and actual dimensions,
  and whether a fixed, requested-dimension, default, or fallback start was
  selected.
- `factory` records the factory name, supplied arguments, and effective
  arguments after defaults.
- `provenance` records the installed funconstrain version. A source commit is
  marked unavailable unless it was supplied by the caller or recorded in the
  installed package metadata.
- `reference` retains `fmin` and `xmin` as non-callback metadata. Its
  applicability flags are `TRUE`, `FALSE`, or `NA` for applicable, inapplicable,
  or not yet encoded, with a corresponding basis string. The helper encodes the
  documented variable-dimension rules for Extended Rosenbrock, Variably
  Dimensioned, and Chebyquad; it does not infer applicability merely because a
  stored `xmin` has the requested length.

Each adapted case is validated once at its resolved starting parameters.
`fg$fn`, `fg$gr`, and `fg$hs` must be functions; the starting parameters must
be a non-empty numeric vector; and the callbacks must return a finite scalar
objective, an equally sized numeric gradient vector, and a finite numeric
Hessian of the corresponding square dimensions. The Hessian must also be
approximately symmetric using a scale-aware relative tolerance of
`sqrt(.Machine$double.eps)`. The optional combined `fg$fg` callback is not
required by the adapter.

This is a developer-boundary usability check, not a claim that the Hessian is
the derivative of the gradient or that it supplies a useful Newton direction.
Those numerical-integrity and direction checks remain separate benchmark
probes.

Each optimizer run gets fresh wrappers around the available `fg$fn`, `fg$gr`,
`fg$fg`, and `fg$hs` callbacks. Their `n_fn_call`, `n_gr_call`, `n_fg_call`,
and `n_hs_call` counts cover only the optimizer invocation. Adapter validation,
the reported initial objective, and final objective/gradient metrics use the
original callbacks and are excluded. A call to the combined `fg$fg` callback
increments only `n_fg_call`; it is not also recorded as separate objective and
gradient calls.

## Hessian Integrity Probe

`hessian-integrity-probe.R` is a deterministic developer check that is separate
from optimizer benchmarking. It compares the analytic Hessian-vector product
`H v` with the centered directional difference
`(g(x + h v) - g(x - h v)) / (2 h)` of the analytic gradient. The generic
`probe_hessian_integrity()` function accepts any finite numeric probe point;
the command-line entry point exercises only each case's resolved start unless
the extended point set is selected explicitly. Sourcing the script defines its
functions without running a probe.

The default directions come from at most three fixed dimensionless patterns:
uniform, alternating, and linear trend. Duplicate patterns are removed. Each
pattern is scaled coordinate-wise by `pmax(1, abs(x))` and then normalized to
unit Euclidean length. The recorded step is
`h = .Machine$double.eps^(1 / 3) * direction_scale`, so the actual coordinate
perturbations remain relative to the point's coordinate scales.

The CSV records Hessian and gradient shape and finiteness, Hessian symmetry,
direction and step scales, absolute and scaled directional residuals, every
tolerance contribution, the final thresholds, per-direction pass/fail, and an
overall point-level classification.
Vector norms are evaluated after scaling by their largest magnitude to avoid
overflow from squaring large finite entries. A direction cannot pass unless its
geometry, perturbed points, analytic and finite-difference products, residuals,
scales, and thresholds are all finite.
Centered differences also order subtraction and division according to scale:
small denominators use difference-then-divide, large denominators scale each
endpoint gradient first, and steps above half the largest double are divided by
the step and by two sequentially without forming an infinite `2 * h`.
The directional threshold is the sum of fixed absolute and relative terms,
both `100 * sqrt(.Machine$double.eps)`, plus the recorded cancellation allowance
`.Machine$double.eps * gradient_scale / h`. Symmetry uses the adapter's relative
threshold `sqrt(.Machine$double.eps) * max(1, max(abs(H)))`.

Run a bounded resolved-start probe over the exact built-in SPD quadratic,
fixed-dimension Rosenbrock and Brown badly scaled cases, and two requested
dimension-five funconstrain cases:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --out /private/tmp/mize-hessian-integrity-starts.csv
```

The SPD control always runs. If funconstrain is unavailable, the script reports
that fact and writes only the SPD rows; funconstrain remains an optional
developer dependency. A failed directional check is callback or adapter
evidence and must stop that problem from being interpreted as an exact-Newton
benchmark. The probe does not replace a failed Hessian numerically, classify
Newton directions, or establish behavior at points that were not selected.

Use `--point-set extended` to add deterministic optimizer-derived and reference
points. For each case, the script runs mize `NEWTON` from the resolved start
with separate callbacks, More-Thuente, its default `step0`, and `max_iter` set
to one and two. It probes the returned parameters and records both the requested
cap and the optimizer's actual iteration and termination. Point collection and
probing use the original callbacks, not the benchmark harness's counted
wrappers, so they do not alter or claim benchmark callback counts.

The extended point set also considers the exact SPD minimizer and each external
`xmin`. An external reference is eligible only when `xmin_applicable` is
exactly `TRUE` and its length equals the resolved dimension. `FALSE` and `NA`
are never promoted from shape alone. Exact duplicate vectors are probed once;
the point-selection manifest records excluded and duplicate candidates, their
basis, and the selected point to which a duplicate was mapped.

Extended runs require a separate manifest path so exclusions are not mixed
with directional probe rows. Active output paths are resolved during argument
parsing and must differ; an alias is rejected before point collection or file
output. Explicit output paths must be non-empty. Output-file symbolic links,
including dangling links, are also rejected rather than followed:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

The directional CSV repeats the optimizer configuration and point provenance
needed to interpret every selected point. The selection manifest has one row
per candidate, including point values, eligibility, requested/actual iteration,
termination, reference applicability and dimension match, selection status,
candidate dimension, and explicit deduplication target. This point-coverage
command does not add eigenspectrum classification unless `--spectrum-out` is
selected, and it never adds Newton/TN direction comparison, line-search
experiments, or performance conclusions.

### Hessian Eigenspectrum Classification

Use `--spectrum-out PATH` with the extended point set to write a separate
point-level eigenspectrum artifact. Only unique selected points whose existing
directional rows all have `probe_pass = TRUE` are classified. The directional
gate and both existing CSV schemas remain unchanged.

The integrity gate does not retain its analytic Hessian, so spectrum
classification reevaluates `fg$hs` exactly once at each eligible point. The
artifact records this as `reevaluated_after_integrity_gate` and records one
`spectrum_hs_calls`; the original callbacks are used directly, so this is not
mixed into the benchmark harness's external callback counts.

An exactly symmetric Hessian is sent to the eigensolver unchanged, preserving
representable subnormal entries. Otherwise, each transposed pair is averaged by
summing before halving when the sum is representable and by adding halves only
when the sum could overflow. For spectral scale
`s = max(abs(lambda))`, the absolute tolerance term is zero and the relative
term is `sqrt(.Machine$double.eps) * s`. An eigenvalue is positive above the
resulting tolerance, negative below its negative, exact zero only when it equals
zero, and numerically unresolved when it is nonzero but remains inside the
tolerance band. No external case changes this rule.

The exhaustive point classifications are `positive_definite`,
`positive_semidefinite`, `zero`, `negative_semidefinite`,
`negative_definite`, `indefinite`, `numerically_unresolved`, and
`calculation_failed`. Indefinite requires both resolved positive and negative
curvature; negative-definite and negative-semidefinite spectra are represented
directly. The CSV records algebraic and absolute eigenvalue extrema, all four
inertia counts, `has_resolved_positive_curvature` and
`has_resolved_negative_curvature`, scale and tolerance terms, and a separate
singularity state. A false resolved-curvature flag does not claim that an
unresolved eigenvalue has the opposite sign. Exact zero modes are singular even
when another mode is unresolved or the sign class is semidefinite or
indefinite.

`spectral_condition_estimate` is
`max(abs(lambda)) / min(abs(lambda))` for fully resolved nonzero spectra. It is
`Inf` whenever an exact zero mode is present, including when another mode is
unresolved. It is missing for unresolved spectra without an exact zero and for
failed spectra. This is an invertibility estimate based on eigenvalue
magnitudes, not a claim of positive definiteness.

All active file destinations are preflighted together before point collection,
including every selected spectrum, direction, globalization, manifest, and
integrity output. Paths must be non-empty, distinct after path resolution, not
live or dangling final symbolic links, and not existing hard links to the same
file. Existing-file identity uses the platform file identity when available
and otherwise requires the `fs` package; unused output paths remain
dependency-free. When prospective basenames differ only by case, preflight
creates and immediately removes one private temporary probe in their resolved
parent directory. Case-equivalent destinations are rejected on a
case-insensitive filesystem and remain distinct on a case-sensitive one; the
probe is removed before case collection. A bounded spectrum command is:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

The spectrum artifact is diagnostic evidence only. It does not construct a
Newton or TN direction, apply a Hessian safeguard, or support optimizer
performance conclusions.

### Raw Exact-Newton Direction Diagnostics

Use `--direction-out PATH` with `--point-set extended` and
`--spectrum-out PATH` to write one raw exact-Newton row for every unique
integrity-passing point, including points whose spectra are singular,
unresolved, or failed. This layer solves only the unglobalized system
`H p = -g`; it does not invoke a line search or construct mize's public
fallback, the internal safeguarded direction, or TN's approximate direction.

A solve is eligible when the existing spectrum evidence is fully resolved and
nonsingular, with a finite positive condition estimate no larger than
`1 / spectral_relative_tolerance`. Under the fixed spectrum rule that limit is
`1 / sqrt(.Machine$double.eps)`. Positive-definite, resolved indefinite, and
negative-definite spectra can therefore be eligible. Exact-singular and
numerically unresolved spectra retain rows with an explicit ineligibility basis
and no callback or solver work; positive definiteness is not silently treated
as a prerequisite for inspecting a raw Newton system.

At each eligible point, the original analytic `fg$gr` and `fg$hs` callbacks are
reevaluated once after the integrity and spectrum gates. The artifact records
those calls separately and marks them outside the benchmark callback counts.
The Hessian is represented with the same overflow-safe symmetrization rule used
for the spectrum. A finite returned direction counts as a successful solve only
when the residual satisfies
`||H p + g|| <= 100 * sqrt(.Machine$double.eps) *
max(1, ||H p||, ||g||)` and every residual term is finite. The CSV records the
condition limit, solve eligibility and basis, callback and solver status,
direction values and norm, residual scale and tolerance terms, `g' p`,
`p' H p`, and predicted decrease `-(g' p + 0.5 p' H p)`. Slope and predicted
decrease signs are exact floating-point sign classifications; no external case
tunes them.

Run the bounded raw-direction slice with a fourth distinct destination:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --direction-out /private/tmp/mize-raw-newton-direction.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

The direction artifact is diagnostic evidence only. It does not record an
accepted step or actual objective decrease and does not support performance or
safeguard-policy conclusions by itself.

### Current-Public Exact-Newton Direction Diagnostics

Use `--public-direction-out PATH` with the extended point set, spectrum output,
and raw `--direction-out` to write one row for mize's current public
exact-Newton policy at every selected integrity-passing point. Raw solve
eligibility does not gate public evaluation: singular and numerically
unresolved spectra still retain public-direction rows.

The probe loads mize, retrieves its actual internal `newton_direction()`
constructor, constructs `newton_direction(try_safe_chol = FALSE)`, and calls
the real direction object's `calculate()` method directly with no line search.
It evaluates the original analytic gradient once to seed the direction cache.
The public implementation then invokes one capturing `fg$hs` wrapper, which
retains the Hessian value used by the direction without a reporting
reevaluation. The artifact records these gradient and Hessian calls explicitly
and marks them outside benchmark callback accounting. No objective callback is
invoked.

The retained `direction_reason` identifies `hessian_solve`,
`cholesky_fallback`, or `direction_check_fallback`. A companion provenance
field spells out ordinary Cholesky solve, Cholesky failure followed by steepest
descent, or direction-check failure followed by steepest descent. Callback,
validation, invocation, returned-direction, and provenance failures receive
explicit non-success rows instead of aborting the artifact.

For a finite public direction `p`, the CSV records its values, overflow-safe
norm, `g' p`, descent status, `p' H p`, predicted decrease
`-(g' p + 0.5 p' H p)`, exact signs, and diagnostic completeness. Model
metrics use the Hessian captured from the real public call and the same
overflow-safe symmetrization convention as the spectrum and raw artifacts; no
extra Hessian callback is made.

When the raw row has `solve_success = TRUE`, the artifact compares the public
and raw vectors under one fixed rule. For
`d = ||p_public - p_raw||` and
`s = max(1, ||p_public||, ||p_raw||)`, a match requires
`d <= 100 * sqrt(.Machine$double.eps) * s`. The difference norm, scale,
tolerance terms, threshold, scaled difference, and match status are all
recorded. When the raw solve is unavailable, the public row remains and the
comparison is explicitly unavailable with the raw solve status as its basis.

Run the bounded current-public slice with a fifth distinct destination:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --direction-out /private/tmp/mize-raw-newton-direction.csv \
  --public-direction-out /private/tmp/mize-public-newton-direction.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

This artifact is current-policy direction evidence only. It does not invoke a
line search, construct `safe_chol()` or any safeguarded direction, evaluate TN,
record accepted-step behavior, or justify a runtime/default change.

### Safeguarded Exact-Newton Direction Diagnostics

Use `--safeguarded-direction-out PATH` with the extended point set and all
three predecessor outputs to write one row for the internal safeguarded
exact-Newton direction at every selected integrity-passing point. The option
requires `--spectrum-out`, raw `--direction-out`, and
`--public-direction-out`, because the new artifact compares with both raw and
current-public evidence. Numerically unresolved spectra and raw-ineligible
points remain evaluated and retain safeguarded rows.

The probe calls the same loaded package's actual
`newton_direction(try_safe_chol = TRUE)$calculate` implementation directly.
Its internal ordinary Cholesky attempt, fixed `1e-10` eigenvalue-floor repair,
solve, and final direction check therefore run without copying or tuning the
algorithm in the probe. As in the current-public layer, the original analytic
gradient seeds the cache once and a single capturing Hessian callback supplies
the implementation and model reporting. Both counts are explicit and outside
benchmark accounting; no objective callback or line search is invoked.

The returned `direction_reason` retains the implementation's branch. Because
that field does not distinguish an ordinary factor from a repaired factor, the
probe also replays only ordinary `chol()` on the captured validated Hessian,
without another callback or a separate `safe_chol()` call. Its result and the
actual branch distinguish ordinary solves, eigenvalue-floor repair solves,
ordinary candidates rejected by the final direction check, and a failed repair
whose `cholesky_fallback` reason survives. When ordinary Cholesky fails and the
returned reason is `direction_check_fallback`, the final check may have
overwritten either a repaired solve or a failed-repair steepest-descent
fallback. The probe does not infer which occurred: repair success is missing,
provenance and evaluation status are explicitly ambiguous, and the finite
returned direction and its metrics remain available. The artifact also records
the ordinary factorization error, whether repair was required and attempted,
and the fixed repair floor.

For each finite safeguarded direction, the CSV records values, overflow-safe
norm, `g' p`, descent status, `p' H p`, predicted decrease
`-(g' p + 0.5 p' H p)`, exact signs, and numerical completeness. Metrics reuse
the captured Hessian and existing overflow-safe symmetrization convention.
Callback, validation, invocation, direction-return, and provenance failures
retain non-success rows rather than aborting later points. Collapsed
post-check provenance retains its finite direction with an explicit ambiguous
status rather than being reported as a successful repair.

Raw/safeguarded comparison is available only for a successful finite raw
direction. Current-public/safeguarded comparison is available only when both
directions succeeded. Each applies the same fixed rule: for difference norm
`d` and scale `s = max(1, ||p_a||, ||p_b||)`, a match requires
`d <= 100 * sqrt(.Machine$double.eps) * s`. Both comparisons report scale,
tolerance terms, threshold, scaled difference, and match status. Unavailable
raw or current-public evidence retains its actual basis; safeguard behavior is
never used to infer a missing predecessor direction.

Run the bounded safeguarded slice with a sixth distinct destination:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --direction-out /private/tmp/mize-raw-newton-direction.csv \
  --public-direction-out /private/tmp/mize-public-newton-direction.csv \
  --safeguarded-direction-out /private/tmp/mize-safe-newton-direction.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

This is unglobalized internal-policy evidence only. It does not invoke a line
search or TN, record accepted steps or actual objective decrease, make
performance claims, change defaults or runtime behavior, or add a public API.

### Truncated-Newton Direction Diagnostics

Use `--tn-direction-out PATH` with the extended point set and all completed
Hessian/Newton outputs to write one row for mize's current default TN direction
at every selected integrity-passing point. The option requires
`--spectrum-out`, raw `--direction-out`, `--public-direction-out`, and
`--safeguarded-direction-out`, because the TN vector is compared independently
with all three predecessor direction artifacts. Numerically unresolved spectra
and raw-ineligible points remain evaluated and retain TN rows.

The probe loads the package's actual internal
`tn_direction(init = 0, exit_criterion = "curvature", preconditioner = "")`
constructor and invokes the real direction object's `calculate()` method at
iteration one. This is the default unpreconditioned, zero-initialized TN
direction policy. The call receives a preseeded analytic gradient and infinite
callback budgets; no objective callback, line search, or accepted step is
reachable.

TN does not call the analytic Hessian. Its `bd_approx()` helper obtains each
approximate Hessian-vector product from one forward analytic-gradient callback.
The artifact records the one gradient seed separately from every attempted HVP
gradient call, successful return, shape-valid return, finite return, and the
inner gradient count returned by mize. TN's analytic-Hessian count is therefore
zero. After a finite direction exists, one separately labeled analytic-Hessian
callback supplies only the common `p' H p` and predicted-decrease reporting
metrics. All of this work is explicitly outside benchmark callback accounting.
The returned inner count is usable only when it is a finite, nonnegative,
dimensionless scalar integer and exactly matches the observed shape-valid HVP
gradient returns.
Missing, malformed, or mismatched counts retain a finite direction but make
callback accounting and inner provenance explicitly incomplete.

The current inner solver returns its direction and gradient count but does not
retain an iteration count or stop reason. During each real HVP callback, the
probe observes the loaded `bd_approx()` and `tn_inner_cg()` call-frame state,
then replays only the decisive post-callback predicates on the captured state
and returned gradient. This classifies stationary or unusable-gradient exits,
nonpositive or unusable curvature, unusable step coefficients, residual
convergence or unusable residuals, and the fixed 40-iteration limit without
generating another iterate or callback. A separate field records whether
`tn_direction()`'s final descent check replaced the inner result with `-g`.
If the current internal state cannot be observed consistently, the valid
direction remains available but inner provenance is explicitly incomplete;
the probe does not infer a stop reason from callback counts or vector equality
alone.

For each finite TN direction `p`, the CSV records its values, overflow-safe
norm, `g' p`, descent status, analytic-model `p' H p`, predicted decrease
`-(g' p + 0.5 p' H p)`, exact signs, callback/provenance completeness, and
overall numerical completeness. Throwing or malformed seed/HVP callbacks,
direction invocation and return failures, missing or malformed returned inner
counts, observer inconsistency, and reporting-Hessian failures retain explicit
point-level statuses rather than aborting the artifact. Shape-valid nonfinite
HVP gradients follow and record the implementation's actual zero-product
fallback policy.

Raw/TN, current-public/TN, and safeguarded/TN comparisons are independently
available only when both relevant vectors exist and are finite. All three use
the same fixed rule as the earlier direction artifacts: for difference norm
`d` and scale `s = max(1, ||p_a||, ||p_b||)`, a match requires
`d <= 100 * sqrt(.Machine$double.eps) * s`. Each comparison reports its scale,
tolerance terms, threshold, scaled difference, and match status. A failed raw
solve stays `unavailable_raw_solve`; no relationship is inferred from the
spectrum or another direction.

Run the bounded TN slice with a seventh distinct destination:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --direction-out /private/tmp/mize-raw-newton-direction.csv \
  --public-direction-out /private/tmp/mize-public-newton-direction.csv \
  --safeguarded-direction-out /private/tmp/mize-safe-newton-direction.csv \
  --tn-direction-out /private/tmp/mize-tn-direction.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

This is unglobalized point-level direction evidence only. It does not invoke a
line search, accept a step, measure actual objective decrease or performance,
change defaults or runtime behavior, add public inner-solver diagnostics, or
add a public API.

### One-Step Globalization Diagnostics

Use `--globalization-out PATH` with the extended point set and all seven
predecessor outputs to write one row for each of the raw, current-public,
safeguarded, and TN directions at every selected integrity-passing point. The
artifact therefore has four rows per point. A raw direction rejected by the
earlier stable-solve gate still retains an explicit zero-callback unavailable
row, while the other three directions at that point are evaluated normally.

Every mathematically available direction is passed independently to the real
step-size object returned by mize's internal `more_thuente_ls()`. The fixed
policy is the current first-step exact-Hessian configuration: `c1 = 1e-4`,
`c2 = 0.9`, exact Armijo, strong curvature, Rasmussen first-alpha
initialization with `try_newton_step = TRUE`, quadratic later-alpha
initialization, 20 objective-limited trials, infinite gradient/combined and
alpha bounds, and cubic safeguarding disabled. The same policy is used for all
four direction families. This is a controlled comparison choice, not a new
default or a claim that the policy is best for every family.

Each available direction gets one isolated objective seed and one isolated
gradient seed at the selected point. Trial evaluations use transparent wrappers
around the original separate `fg$fn` and `fg$gr`; the combined callback is not
supplied. The CSV records seed, trial, and total calls; successful returns;
shape-valid and finite results; the step-size object's `ls_nf`/`ls_ng`; the
returned optimizer increments; and affirmative count matches. All calls are
outside benchmark accounting, and this slice makes no Hessian or HVP callback.
Callback accounting is complete only when every returned trial objective and
gradient is shape-valid and finite and every returned count matches the
observed attempts. A shape-valid nonfinite trial return remains an explicit
`trial_objective_value_nonfinite` or `trial_gradient_value_nonfinite`
non-success even when More-Thuente safely returns `nonfinite_recovery` and
selects `no_step`; the returned reason and outcome are retained separately.
The real zero-direction short-circuit has one seed of each kind and no trial
call. Its absent line-search reason/outcome remain absent, while the artifact
termination identifies `zero_direction_short_circuit` explicitly.

The accepted alpha is the alpha selected by the real step-size calculation. It
can be zero when the policy selects no step or short-circuits at a stationary
direction; `line_search_step_accepted` separately requires a positive alpha and
a `wolfe` or `improving_fallback` outcome. The artifact retains the actual
line-search reason and selected-point outcome, selected parameters, objective,
gradient and slope, actual decrease `f(x) - f(x + alpha p)`, and exact signs.
Selected objective/gradient data are reused from the returned line-search state
or the seed at a zero step, so reporting adds no callback.

The predecessor's full-step slope, curvature, and predicted decrease remain in
the row. The line-search seed gradient recomputes `g' p`; its agreement with the
predecessor slope uses
`100 * sqrt(.Machine$double.eps) * max(1, abs(s_seed),
abs(s_predecessor))`. With selected alpha `a`, the scaled model prediction is
`-(a * g' p + 0.5 * a^2 * p' H p)`, using the new seed slope and the
predecessor artifact's curvature without another Hessian callback. The
actual-to-predicted ratio is available only when the slope comparison matches,
the actual and predicted decreases are finite, and the predicted decrease is
strictly positive. Stationary or rejected zero steps therefore retain an
explicit `predicted_decrease_not_positive` basis instead of an invented ratio.
A slope mismatch makes numerical and overall diagnostics incomplete even when
the real line search returned a usable selected point.

Callback/validation errors, malformed line-search returns, count mismatches,
Wolfe acceptance, improving fallback, no-step selection, and the zero-direction
short-circuit all retain separate statuses. Direction availability,
line-search success, callback accounting, selected-point completeness,
numerical completeness, and overall diagnostic completeness are independent
fields, so one failed direction does not abort another row or point.

Run the bounded globalization slice with an eighth distinct destination:

```sh
R_LIBS=/private/tmp/mize-funconstrain-metadata-lib \
  Rscript scripts/hessian-integrity-probe.R \
  --dimension 5 \
  --funconstrain-cases rosen,brown_bs,var_dim,chebyquad \
  --point-set extended \
  --selection-out /private/tmp/mize-hessian-point-selection.csv \
  --spectrum-out /private/tmp/mize-hessian-spectrum.csv \
  --direction-out /private/tmp/mize-raw-newton-direction.csv \
  --public-direction-out /private/tmp/mize-public-newton-direction.csv \
  --safeguarded-direction-out /private/tmp/mize-safe-newton-direction.csv \
  --tn-direction-out /private/tmp/mize-tn-direction.csv \
  --globalization-out /private/tmp/mize-one-step-globalization.csv \
  --out /private/tmp/mize-hessian-integrity-extended.csv
```

This artifact is bounded one-step numerical evidence. It does not run a second
outer iteration, rank methods, compare performance, change a line-search or
direction policy, add runtime diagnostics, or change a public API.

### Dependencies

Built-in benchmark cases use base R plus packages that ship with R:

- `stats`
- `datasets`

When running from a source checkout, install `pkgload` so the script loads the
local source tree. This is preferred for the current-public, safeguarded, TN,
and globalization artifacts because they retrieve the loaded package's
internal `newton_direction()`, `tn_direction()`, or `more_thuente_ls()`
implementation:

```sh
Rscript -e 'install.packages("pkgload")'
```

If `pkgload` is unavailable, the script falls back to an installed `mize`
package.

The optional `funconstrain` cases require the separate `funconstrain` package.
Install it only if you want those external benchmark problems:

```sh
Rscript -e 'install.packages("pak"); pak::pak("jlmelville/funconstrain")'
```

If `funconstrain` is not installed, the script reports that and skips those
cases without failing. You can also skip them explicitly with
`--no-funconstrain`.

No benchmarking package such as `bench` or `microbenchmark` is required; the
harness uses base R timing.

### Smoke Run

Use a tiny smoke run to confirm that the harness itself is wired correctly:

```sh
Rscript scripts/benchmark-optimizers.R \
  --cases rosenbrock,spd-quadratic,mmds-eurodist \
  --spd-n 5 \
  --max-iter 2 \
  --reps 1 \
  --line-search More-Thuente \
  --step0 default \
  --out /tmp/mize-benchmark-smoke.csv
```

This command should produce a CSV file, but the numbers are not benchmark
evidence for optimizer quality. They are too small and too noisy for that.

### Output Policy

By default, results are written as CSV to standard output. Use `--out PATH` to
write a file. Prefer temporary locations such as `/tmp` or an intentionally
ignored local results directory for exploratory runs.

Do not commit large generated benchmark output by default. Commit summarized
results only when their review purpose, command, environment, and interpretation
are documented deliberately.

The CSV columns include:

- problem case and source
- actual problem dimension and initial objective value
- optimizer and method
- `mize` line search and `step0` setting
- callback mode (`combined` or `separate`)
- final objective
- final gradient norm
- `nf` and `ng`: mize's existing reported counts for mize rows, and the
  existing external objective/gradient counts for `stats::optim()` rows
- external `n_fn_call`, `n_gr_call`, `n_fg_call`, and `n_hs_call` callback
  invocation counts for every row
- elapsed wall time
- iteration count where available
- termination and failure mode
- warnings and errors

### Examples

Run the default built-in problem grid and skip optional external problems:

```sh
Rscript scripts/benchmark-optimizers.R \
  --no-funconstrain \
  --out /tmp/mize-benchmark.csv
```

Run a small line-search and first-step sensitivity check on Rosenbrock:

```sh
Rscript scripts/benchmark-optimizers.R \
  --cases rosenbrock \
  --max-iter 100 \
  --reps 3 \
  --line-search More-Thuente,Hager-Zhang \
  --step0 default,1 \
  --out /tmp/mize-rosenbrock-benchmark.csv
```

Run optional `funconstrain` cases after installing `funconstrain`:

```sh
Rscript scripts/benchmark-optimizers.R \
  --cases rosenbrock \
  --funconstrain-cases rosen,chebyquad \
  --out /tmp/mize-funconstrain-benchmark.csv
```

Run a bounded exact-Newton wiring check across a built-in SPD control, fixed-
dimension external problems, requested variable dimensions, and one rejected
dimension that falls back to its factory default:

```sh
Rscript scripts/benchmark-optimizers.R \
  --cases spd-quadratic \
  --funconstrain-cases rosen,brown_bs,meyer,var_dim,chebyquad,ex_powell \
  --spd-n 5 \
  --methods mize-newton \
  --max-iter 1 \
  --line-search More-Thuente \
  --step0 default \
  --callbacks combined,separate \
  --out /tmp/mize-newton-smoke.csv
```

This is a harness smoke test, not derivative-integrity or Newton-direction
evidence. `n_hs_call` confirms that the optimizer consumed the mapped Hessian.
The CSV records the actual dimension, while the adapter's `start` metadata and
the fallback message retain the requested-versus-resolved dimension decision.

The `--spd-n` value is also requested from variable-dimension `funconstrain`
factories. If a factory rejects that dimension, the harness reports the
rejection and uses the factory's documented default starting-point dimension;
the actual dimension is recorded in the CSV.

Compare combined and separate callbacks on the same external cases:

```sh
Rscript scripts/benchmark-optimizers.R \
  --cases rosenbrock \
  --funconstrain-cases rosen,chebyquad \
  --callbacks combined,separate \
  --line-search More-Thuente,Hager-Zhang \
  --step0 default \
  --out /tmp/mize-callback-benchmark.csv
```
