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

### Funconstrain Adapter Schema

The sourceable helper preserves the benchmark's `case$name`, `case$source`,
`case$par`, and `case$fg` fields. It maps funconstrain callbacks as follows:

- `fn` to `fg$fn`
- `gr` to `fg$gr`
- optional `fg` to `fg$fg`
- optional `he` to mize's `fg$hs` name

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
- function and gradient evaluation counts
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
