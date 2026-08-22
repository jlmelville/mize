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

All active file destinations are preflighted together before point collection:
`--out`, `--selection-out`, and `--spectrum-out` must be non-empty, distinct
after path resolution, not live or dangling final symbolic links, and not
existing hard links to the same file. Existing-file identity uses the platform
file identity when available and otherwise requires the `fs` package; unused
output paths remain dependency-free. When prospective basenames differ only by
case, preflight creates and immediately removes one private temporary probe in
their resolved parent directory. Case-equivalent destinations are rejected on
a case-insensitive filesystem and remain distinct on a case-sensitive one; the
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

### Dependencies

Built-in benchmark cases use base R plus packages that ship with R:

- `stats`
- `datasets`

When running from a source checkout, install `pkgload` so the script benchmarks
the local source tree:

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
