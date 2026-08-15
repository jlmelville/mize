Audit: jlmelville/mize

One correction to the premise up front: mize is pure R — there's no C++ anywhere in the repo (no src/, no compiled code). Its only hard dependency is methods. Everything below is an R-code audit.

What the repo is

mize is a CRAN package (v0.2.5 released; repo is at dev version 0.2.5.9000) implementing unconstrained gradient-based optimization: SD, BFGS, SR1, L-BFGS, CG (10 update variants), truncated Newton, Newton, NAG, momentum, and DBD, combined with five line searches (More-Thuente, Rasmussen, Schmidt/minFunc, Hager-Zhang, backtracking/bold-driver). It offers both a black-box mize() call and a stateful make_mize()/mize_step() API. ~10,700 lines of R across 26 files, ~7,300 lines of tests. The architecture is an unusual but coherent hook/lifecycle ("aspect-oriented") event system in life_cycle.R that wires directions, step sizes, momentum, validation, and caching together.

Health snapshot (verified locally): full testthat suite passes cleanly (all 24 test files, thousands of assertions); R CMD check passes with no package-related warnings or notes; test coverage is 86.9% overall; CI is strong (R CMD check across 6 OS/R-version configs with pinned action SHAs, plus lint, coverage, rhub, and pkgdown workflows). Docs (man/) are in sync with roxygen sources. Recent history shows a large, active cleanup burst (new status fields, gradient checker, many new tests, safeguards) — and, as it turns out, one of the new safeguards introduced a regression.

________________________________

A. Confirmed bugs (reproduced against the current head)

A1. Crash: L-BFGS memory can be empty when lbfgs_solve() is called — regression from commit 9fa0c6f ("skip bad BFGS updates", Jul 6 2026, unreleased). lbfgs_memory_update() now silently skips updates failing the curvature test, but lbfgs_solve() still assumes at least one stored update: for (i in length(rhos):1) becomes 0:1 and indexes rhos[[0]]. If the first update is skipped (zero/negative curvature — e.g. a locally linear region, or any non-Wolfe line search), the optimizer errors. Reproduced through the public API in three method configurations:

lin <- list(fn = function(x) sum(2*x), gr = function(x) c(2, 2))
mize(c(1,1), lin, method = "L-BFGS", line_search = "constant", step0 = 0.1, ...)
# Error: attempt to select less than one element in integerOneIndex

The same crash fires for method = "TN", preconditioner = "l-bfgs" and method = "CG", preconditioner = "l-bfgs", since both feed lbfgs_memory_update() results straight into lbfgs_solve(). Fix is a one-line guard: fall back to the identity/steepest-descent guess when the memory is empty.

A2. Crash: SR1 references out-of-scope variables on a non-descent first direction. In sr1_direction$calculate, hm, sm, ym are only defined inside the if (!is.null(gm_old)) block, but the try_bfgs fallback that uses them runs whenever descent >= 0 — including iteration 1 with a user-supplied non-positive-definite fg$hi:

fg <- list(fn = function(x) sum(x^2), gr = function(x) 2*x, hi = function(x) diag(-1, 2))
mize(c(1,1), fg, method = "SR1")   # Error: object 'ym' not found

A3. cubic_extrapolate() silently returns NULL with default arguments. The function body is a single if (ignoreWarnings) { ... } with no else branch, so ignoreWarnings = FALSE (the default) returns NULL. The only current caller passes TRUE, so this is latent rather than live, but it's a plainly broken code path and the single test in test_polation.R only exercises the TRUE path.

A4. Gradient evaluations can go uncounted in mize_step_summary(). When a gradient norm is needed and the gradient isn't cached, fg$gr(par) is invoked, but counts$gr is only incremented when grad_is_first_stage(opt) is true. For Nesterov-momentum-first optimizers the evaluation really happens yet is never counted (or cached), so ng under-reports and max_gr/max_fg budgets can be silently exceeded. Skipping the cache there is defensible (the gradient position differs under Nesterov); skipping the count is not.

B. Correctness risks and design smells

B1. Unvalidated cache reads (stale-value hazard). The cache design pairs each value (fn_curr, gr_curr) with an _iter validity marker, and most code checks it via has_fn_curr()/has_gr_curr(). Three places read cached values without checking validity: check_mize_convergence() reads opt$cache$fn_curr directly for abs/rel tolerance checks; check_gr_conv() reads opt$cache$gr_curr (so a stale gradient can trigger grad_tol/gr_inf, especially via the loop API where users control call order); and the best-result restoration in opt_loop() compares best_fn != opt$cache$fn_curr / best_grn != norm_inf(opt$cache$gr_curr) — a stale cached value can make the "is best already current?" test wrong in either direction. These paths mostly work in the canonical mize() flow because summaries refresh the cache first, but the invariant is enforced by call order rather than by the accessors, which is fragile for the stateful API.

B2. SR1 first-step scaling lacks the curvature guard BFGS received. Commit 9fa0c6f added has_bfgs_curvature() to BFGS's iteration-2 gamma scaling, but the identical scaling in sr1_direction still computes gamma <- dot(sm, ym) / dot(ym) unguarded — a non-positive sm·ym yields a negative/degenerate scaling of the inverse-Hessian estimate.

B3. Operator-precedence ambiguity in register_hook(). The condition join_point == "stage" || join_point == "gradient_descent" || join_point == "momentum" && is.null(stage_type) parses as a || b || (c && d) — the is.null(stage_type) guard applies only to "momentum". If that asymmetry is intentional it deserves parentheses and a comment; if not, it's a latent bug in hook routing.

B4. max_iter detection in check_mize_convergence() uses equality. if (opt$iter == convergence$max_iter) never fires when check_conv_every > 1 and max_iter isn't a multiple of it. opt_loop() compensates with its own loop bound, but users of the documented stateful API (mize_step() + check_mize_convergence()) can iterate past max_iter without ever seeing the termination flag. >= would be safer.

B5. Minor: update_progress() grows a data frame with rbind per check — O(n²) for long, frequently-logged runs; in the restore-best branch of opt_loop(), opt is not reassigned from mize_step_summary()'s return (harmless today because counts are read from step_info, but inconsistent with every other call site); mize() validates max_iter/max_fn/etc. after the full optimizer has been constructed; normalize() computes the norm twice.

C. Documentation defects

C1. NEWS.md dev entry misattributes the new function. It announces "check_mize_convergence(), which will compare your an analytic gradient with a finite-different approximation" — that describes the genuinely new check_mize_gradient(); check_mize_convergence() has existed for years. The same entry contains garbled text ("your an analytic", "finite-different", "returns additive status fields more status fields").

C2. abs_tol prose contradicts the code and its own param doc. The Convergence section of ?mize says abs_tol triggers when "the absolute value of the function falls below this threshold" (and only for functions with minimum zero), but check_fn_conv() tests abs(fn_old - fn_new) < abs_tol — a successive-difference criterion, exactly as the @param text ("comparing two function evaluations") says. One of the two descriptions is wrong; the section text should be corrected.

C3. Smaller doc issues: max_fg is described in the Convergence section as "Maximum number of gradient evaluations" (copy-paste from max_gr); dbd_weight "must be an integer between 0 and 1" (should be a real number; code validates a range, not integrality); the comment above is_in_range() calls the inclusive default an "open range" (terminology inverted); assorted typos ("strng", "currents step size", "appears to be not be making progress").

D. Testing

The suite is genuinely strong — behavioral invariants (test_optimizer_invariants.R), regression tests, line-search condition tests, and loop-API tests — and everything passes. Gaps line up with the findings above: no tests for empty L-BFGS memory / skipped-first-update paths (A1), SR1 with user-supplied hi (A2), cubic_extrapolate's default path (A3), or gradient counting under Nesterov-first stages (A4). Coverage is weakest in util.R (51%, mostly logging helpers) and the vendored line searches (hager_zhang.R 70%, rasmussen.R 78%, more_thuente.R 82%) — precisely the numerically delicate code where untested branches are failure-recovery paths.

E. Provenance, licensing, housekeeping

Translated third-party algorithms vs. BSD-2 license. schmidt.R is explicitly "a translation of Mark Schmidt's minFunc line search code," rasmussen.R derives from Carl Rasmussen's minimize.m, and more_thuente.R is based on Dianne O'Leary's Matlab translation of MINPACK. The package license is BSD 2-clause with a single copyright holder, and none of the upstream authors appear in Authors@R (e.g. as ctb/cph). I'm not asserting incompatibility — attribution in comments is present, and the package has been on CRAN for years — but the license terms of the originals (minFunc in particular) are worth verifying and documenting, and CRAN policy generally expects derived-code authors to be credited in Authors@R.
docs/ contains four stray files (nesterov.Rmd/html, qhm.html/md) that duplicate/shadow vignettes/articles/; pkgdown now builds to pkgdown-site/, so these look like stale leftovers.
.lintr excludes R/RcppExports.R / src/RcppExports.cpp, which don't exist — harmless template residue, but consistent with the "no C++" point.

________________________________

Overall assessment

This is a mature, well-engineered package in unusually good shape for a solo-maintained numerical library: excellent CI, a serious test suite, clean checks, synced docs, and thoughtful numerical safeguards throughout the line searches. The hook/lifecycle architecture is exotic and hard to statically reason about, but it's consistently applied and well tested at the behavioral level.

The most important takeaway is that the unreleased dev head contains a real regression: the July 2026 "skip bad BFGS updates" safeguard made L-BFGS (and the L-BFGS preconditioner used by CG and TN) crash whenever the first curvature update is rejected (A1). Together with the SR1 scope bug (A2), both crashers live in first-iteration/degenerate-curvature edge paths that the otherwise-thorough test suite doesn't reach — worth fixing and pinning with regression tests before the next release. The remaining findings are secondary: cache-validity discipline for the stateful API (B1, B4), the gradient-count leak (A4), and a round of documentation corrections led by the NEWS misattribution (C1) and the abs_tol contradiction (C2).

Happy to turn any of these into minimal reproducible issue reports or suggested patches if you'd like.