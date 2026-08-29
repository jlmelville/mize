# Reverse-dependency check results

Last updated: 2026-08-28

## Result

The CRAN reverse-dependency comparison found no regressions.

| Reverse dependency | Version | CRAN mize | Development mize | New problems |
|:-------------------|:--------|:----------|:-----------------|:-------------|
| ctsem | 3.11.1 | OK | OK | None |

- Direct CRAN reverse dependencies checked: 1
- New problems: 0
- Packages that failed to check: 0

Both `ctsem` check logs ended with `Status: OK` (0 errors, 0 warnings, and
0 notes).

## Comparison details

- Baseline: CRAN `mize` 0.2.5.
- Candidate: development `mize` 0.2.5.9002 at commit `372d96a`.
- Reverse dependency: CRAN `ctsem` 3.11.1 (source MD5
  `062c3b8ce38aaf8d16359fc68d025676`).
- Platform: Ubuntu 26.04.1 LTS, x86_64, R 4.5.2.
- Check date: 2026-08-28.

The subsequent release-preparation commit, `d223d36`, changed only the package
version, release notes and comments, and two vignette sources. It did not change
R implementation code or tests, so this comparison covers the runtime and API
changes in mize 0.3.0. The exact 0.3.0 candidate also passed a local
`R CMD check --as-cran` and the package's GitHub Actions checks.

## Setup cost

Although there was only one direct reverse dependency, preparing its dependency
closure required 62 packages; 33 were marked as needing compilation. The
isolated runner was seeded with those packages plus mize. The disposable cache,
libraries, sources, and results occupied approximately 1.9 GB, and the old and
new `ctsem` installation logs were approximately 24 MB each.

The generated libraries, caches, and check trees are intentionally not tracked
in this repository. This file records the durable result. A new comparison is
worth running if mize's implementation or public API changes again, or if its
CRAN reverse-dependency set changes.
