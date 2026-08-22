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
