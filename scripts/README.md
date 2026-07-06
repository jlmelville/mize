# Developer Scripts

This directory contains developer utilities that are useful when working from
the source repository. These scripts are not part of the package API and are
excluded from R package builds.

## Optimizer Benchmark Harness

`benchmark-optimizers.R` compares `mize` optimizers with `stats::optim()` on a
small set of smooth unconstrained optimization problems. It is a permanent
developer benchmark harness, not a pass/fail package validation command.

Use `testthat::test_local()` or `R CMD check` to check that the package works
as expected. Use this benchmark script when you want local evidence about
optimizer behavior, evaluation counts, wall time, termination modes, or
sensitivity to `mize` line-search and `step0` settings.

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
- optimizer and method
- `mize` line search and `step0` setting
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
