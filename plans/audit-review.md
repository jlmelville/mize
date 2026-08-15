## Bottom line

The repo has a bounded correctness backlog, not a roadmap problem.

I would accept the crashers, state-reset defects, stale-cache behavior, budget overruns, rejected-line-search steps, malformed callback handling, the BFGS update improvement, exact documentation corrections, and the licensing gap. I would reject every proposed new optimizer capability, public diagnostic layer, control-object hierarchy, lifecycle rewrite, provenance system, and generalized quality campaign.

No implementation files were changed during this review.

## Assessment of the reviewers

The Sol audit was impressively accurate despite lacking R. Its strongest findings—the L-BFGS, SR1, initialization, budget, cache, and line-search problems—survive inspection. Its weakness is classification: it repeatedly puts correctness defects, performance improvements, architectural preferences, and speculative features on the same “fix-first” continuum.

The Fable audit is more restrained and, overall, the better execution review. It reproduced the main crashers and found several additional real defects: `cubic_extrapolate()`, gradient undercounting, stale cache reads, and documentation problems. Its operator-precedence observation is not established as a bug, and some minor findings do not justify work.

The audit baseline was `92fffdf`. Since the Sol snapshot at `d4a3901`, no R implementation or tests had changed—only CI, ignore files, and `DESCRIPTION`. Consequently, the algorithmic findings were current when reviewed.

At that baseline, the suite passed all 2,359 assertions. That did not rebut either audit: the failures occurred on inputs the suite did not exercise.

## Decisive reproductions

These were reproduced against the audit baseline. The first two rows were
resolved in phase 1; the remaining rows are still queued below.

| Behavior | Baseline result |
|---|---|
| Rejected first L-BFGS pair | Public L-BFGS, CG preconditioning, and TN preconditioning all crash |
| SR1 with indefinite initial `hi` | Crashes with `object 'ym' not found` |
| Explicit `mize_init(max_iter = 30)` after factory default 100 | Retains 100 |
| Reinitialize a terminated optimizer | Terminal fields survive and the next step does nothing |
| Call `mize_step()` before initialization | Silently returns `numeric(0)` parameters |
| `cubic_extrapolate()` with its default argument | Returns `NULL` |
| Bold Driver with `max_fn = 1` on \(x^4\) | Uses two function calls and moves from \(f=1\) to \(f=81\) |
| Schmidt Armijo with one line-search evaluation | Applies the known failing step, also moving to \(f=81\) |
| Best-result restoration after a worsening constant step | Returns \(f=9\) although the initial \(f=1\) was recorded |
| TN with `max_fg = 1` | Performs one function and two gradient calls |
| Invalid one-element gradient for two parameters | Silently recycles it and updates both parameters |
| Vector-valued objective or wrong-sized inverse Hessian | Fails later with opaque base-R errors |

The BFGS algebraic expansion agreed with the current implementation to about \(10^{-14}\). It was roughly 4–6.5 times faster at dimensions 100–500 in a small same-process comparison. That makes it a justified performance fix, though not a correctness or release-blocking defect.

## Clarifications after owner review

### Licensing correction: resolve the source chain or rewrite from the paper

`mize` can remain BSD-2, but the earlier conclusion about the More–Thuente
source was too optimistic. The permissive notice in Netlib's MINPACK directory
does not establish permission for the More–Thuente line search. This is a
repository-maintenance assessment, not legal advice.

- Mark Schmidt distributes the software on the `minFunc` page under a
  two-clause FreeBSD-style license. Those terms are compatible with `mize`'s
  BSD-2 license. They require the complete Schmidt copyright notice,
  conditions, and disclaimer to accompany source and binary redistributions.
- The contemporaneous GPML 4.0 archive from October 2016 also uses the
  two-clause FreeBSD license. Its notice names Carl Edward Rasmussen and Hannes
  Nickisch; `util/minimize.m` separately identifies Carl Edward Rasmussen as
  its author and copyright holder. This is compatible with BSD-2, but the
  upstream notice must be retained.
- Netlib's permissively licensed MINPACK directory contains the nonlinear
  equation and least-squares routines, but not `cvsrch`, `cstep`, `dcsrch`, or
  `dcstep`. The line search is in the separate MINPACK-2 `csrch` distribution.
- The authoritative MINPACK-2 `csrch` README permits copying and modification
  for internal research use, expressly prohibits commercial use without
  permission, and requires its notice to remain. Those terms are not a basis
  for putting a direct translation into a BSD-2 package.
- Dianne O'Leary's official `cvsrch.m` and `cstep.m` sources identify her 1991
  MATLAB translation and the MINPACK originals, but neither the files nor the
  official download page supplies a separate permissive licence. The current
  R file explicitly says it translates those MATLAB files and retains a large
  amount of source-shaped structure and commentary.

"Translate the Fortran and then R-ify it" would not cure this. Translation and
cosmetic restructuring can still copy the source's protected expression. A
genuinely new R implementation based on Moré and Thuente's 1994 paper remains
the definitive fallback: copyright does not protect an algorithm, procedure,
or method as such, although it does protect a paper's and a program's
particular expression.

Before paying that implementation cost, make one bounded attempt to establish
a usable licence chain through O'Leary. Her permission must cover her own
copyrightable contribution to `cvsrch.m` and `cstep.m`; it must also establish
whether she had authority to distribute or license the underlying MINPACK
material. Permission for her MATLAB translation alone cannot license code
owned by Moré, Thuente, Argonne, or another upstream holder.

The revised repair is therefore:

1. Email O'Leary now, while the unrelated correctness phases proceed.
2. Treat the response as sufficient only if it explicitly permits use,
   modification, translation, and source/binary redistribution under terms
   compatible with BSD-2 and resolves the scope of the underlying MINPACK
   rights.
3. If that complete chain is established, retain the implementation and add
   the exact required notices and contributor attribution. If it is not
   established by the time the preceding implementation phases finish,
   proceed with the independent paper-derived rewrite; an unanswered email is
   not a reason to defer the release decision indefinitely.
4. Do not release another package tarball containing the current translation
   before either branch is complete.

Both branches preserve More–Thuente as the default and avoid a disruptive API
removal. Neither requires a new dependency, compiled code, provenance
manifest, or recurring compliance machinery.

There is no licensing reason to remove the Schmidt or Rasmussen line searches:
their upstream terms are explicit and compatible. If the replacement branch is
needed, replacing the current file makes future source package distributions
clean, but the old translation will remain in Git history. Do not turn history
rewriting into an automatic part of this code change; whether to seek
retrospective permission or alter the public repository history is a separate
owner/legal decision.

Primary materials checked:

- CRAN, *Writing R Extensions*, sections 1.1.1–1.1.2:
  <https://cran.r-project.org/doc/manuals/r-release/R-exts.html>
- Mark Schmidt's software copyright terms:
  <https://www.cs.ubc.ca/~schmidtm/Software/copyright.html>
- GPML copyright notice:
  <https://gaussianprocess.org/gpml/code/matlab/Copyright>
- Netlib MINPACK directory and copyright notice (showing both the permissive
  terms and the absence of the line-search routines):
  <https://www.netlib.org/minpack/>
- MINPACK-2 `csrch` notice:
  <https://ftp.mcs.anl.gov/pub/MINPACK-2/csrch/README>
- Dianne O'Leary's official software page and translation sources:
  <https://www.cs.umd.edu/users/oleary/software/>
- Moré and Thuente (1994), *Line Search Algorithms with Guaranteed Sufficient
  Decrease*, DOI 10.1145/192115.192132.
- 17 U.S.C. § 102(b) and U.S. Copyright Office Circular 33 on the distinction
  between methods and their particular expression.

#### O'Leary contact and permission request

O'Leary's [current UMD homepage](https://www.cs.umd.edu/~oleary/) says to email
her last name, without the apostrophe, at `umd.edu`: `oleary@umd.edu`. Current
[UMD Computer Science](https://www.cs.umd.edu/people/oleary) and
[QuICS](https://quics.umd.edu/people/dianne-oleary) profiles also list
`oleary@cs.umd.edu`. Use the former first and the latter only as a fallback
rather than sending duplicate messages.

Suggested subject: **Permission to redistribute your More–Thuente MATLAB
translation in the R package mize**

> Dear Professor O'Leary,
>
> I maintain `mize`, an open-source R package for numerical optimization:
> <https://github.com/jlmelville/mize>. The package is distributed under the
> BSD 2-Clause licence.
>
> Its More–Thuente line search in `R/more_thuente.R` was translated from the
> `cvsrch.m` and `cstep.m` files on your University of Maryland software page.
> During a licensing review I noticed that those files identify you as the
> MATLAB translator but do not state terms for modification and
> redistribution.
>
> Would you be willing to license your copyrightable contributions to
> `cvsrch.m` and `cstep.m` under the BSD 2-Clause licence, including permission
> to use, reproduce, modify, translate, and redistribute them and derivative R
> code in source and binary form? I would retain your preferred copyright and
> attribution notice.
>
> The files also identify the underlying routines as MINPACK work by Jorge J.
> Moré and David J. Thuente. Could you please also tell me whether you had
> permission to distribute the translated routines and whether you have
> authority to license the files as a whole? If your permission covers only
> your MATLAB contribution, any pointer to the relevant upstream permission or
> rights holder would be very helpful.
>
> If you are happy to grant permission, wording such as the following would
> give me a clear record:
>
> "I grant James Melville and recipients of the `mize` package permission to
> use, reproduce, modify, translate, and redistribute my copyrightable
> contributions to `cvsrch.m` and `cstep.m`, including derivative R code, in
> source and binary form under the BSD 2-Clause licence. [This grant covers the
> files as a whole / this grant is limited to my own contribution.] Please use
> the following copyright and attribution notice: [your preferred wording]."
>
> Thank you for making the MATLAB implementation available, and for considering
> this request.
>
> Best wishes,
> James Melville

### BFGS: the quadratic update is the standard implementation

The proposed \(O(n^2)\) implementation is not a new quasi-Newton method and is
not an LLM invention. Full BFGS is a rank-two update of an \(n \times n\)
Hessian or inverse-Hessian approximation, so the expected update cost is
quadratic. The current `mize` expression is itself a standard, factored form of
the inverse-BFGS identity:

\[
H_+ = (I - \rho s y^T) H (I - \rho y s^T) + \rho s s^T.
\]

The problem is only that `mize` asks generic dense matrix multiplication to
evaluate that identity literally. Multiplying the three dense \(n \times n\)
matrices makes this particular implementation cubic. Compute \(Hy\) once and
expand the same identity:

\[
H_+ = H - \rho\{(Hy)s^T + s(Hy)^T\}
      + \{\rho + \rho^2 y^T H y\}ss^T.
\]

That requires one matrix–vector product, dot products, and rank-one/rank-two
outer-product updates, all \(O(n^2)\). SciPy's maintained BFGS implementation
uses this same expanded inverse-Hessian formula and explicitly describes it as
the more efficient implementation of the standard formula; its BLAS `syr` and
`syr2` calls are optimized forms of the same outer-product operations. The
classic Dennis–Moré review gives the factored inverse-BFGS identity, and the
original 1970 BFGS literature establishes the rank-two method.

This changes evaluation order, not mathematics. Floating-point results can
differ in their last few bits, so implementation must retain the existing
curvature guard and `eps` convention and be protected by algebraic-equivalence,
symmetry, descent-direction, and public optimizer tests. The earlier benchmark
is supporting evidence, not the justification for inventing a variant.

Technical references checked:

- Dennis and Moré (1977), *Quasi-Newton Methods, Motivation and Theory*,
  especially the inverse-BFGS identity in equation 7.25.
- Fletcher (1970), *A new approach to variable metric algorithms*,
  DOI 10.1093/comjnl/13.3.317.
- Shanno (1970), *Conditioning of quasi-Newton methods for function
  minimization*, DOI 10.1090/S0025-5718-1970-0274029-X.
- SciPy's `BFGS._update_inverse_hessian()` implementation:
  <https://github.com/scipy/scipy/blob/main/scipy/optimize/_hessian_update_strategy.py>

## Accepted execution plan

Correctness determines the first four phases. Within each phase, fixes should remain independently testable and commit-ready.

### Execution status

Phase 1 is complete, independently reviewed, and accepted. The implementation:

- treats empty L-BFGS history as a normal state and uses the documented initial
  inverse-Hessian guess, fixing public L-BFGS plus the CG and TN L-BFGS
  preconditioners without caller-specific branches;
- makes SR1 use steepest descent when no secant pair exists, applies the BFGS
  curvature guard before initial scaling, and keeps the later BFGS rescue only
  for iterations with a real pair; and
- validates `fg$hi` narrowly in the BFGS, SR1, and L-BFGS paths. Valid finite
  zero or indefinite approximations may fall back to steepest descent;
  malformed or nonfinite results fail with an informative error. Phase 4
  should extend or reuse this boundary rather than introduce a competing
  quasi-Newton validator.

Public regressions cover rejected first pairs in all three L-BFGS consumers,
identity and user-supplied empty-history guesses, SR1 vector and matrix
fallbacks, guarded scaling, and invalid inverse-Hessian results. Focused tests,
the full 2,381-assertion suite, Air, lintr, `git diff --check`, and independent
manual BFGS/L-BFGS/SR1 `hi` probes passed. No generated documentation, public
API, dependency, or later-phase implementation was changed.

Two implementation details mattered: the rejected-pair state itself was valid
and belonged in `lbfgs_solve()`, while SR1's first-iteration failure came from
trying to use `hm`, `sm`, and `ym` before a secant pair existed. Central fixes
covered the public configurations; neither required architectural work.

The audits and this file are the durable record. Future phase handoffs should
be supplied in chat and should update this execution-status section when a
phase completes; do not create separate phase handoff files. Phase 2 is next.

### 1. Fix the direct crashers — complete

- Make empty-history `lbfgs_solve()` return the normal initial inverse-Hessian guess: user `hi` when supplied, otherwise identity behavior.
- Add public regressions for L-BFGS, CG with L-BFGS preconditioning, and TN with L-BFGS preconditioning.
- Make the first-iteration SR1 fallback steepest descent when no secant pair exists.
- Apply the BFGS curvature guard to SR1’s initial scaling.
- Test finite zero/ascent `hi` results through exported APIs and require a
  steepest-descent fallback. Reject nonfinite or malformed `hi` results with an
  informative error; keep that validation local to the inverse-Hessian paths
  needed by this phase, leaving generalized callback validation to phase 4.

These are unequivocal correctness defects and should precede all cleanup.

### 2. Repair the stateful optimizer contract and cache discipline

- Explicit arguments supplied to `mize_init()` must override factory convergence values.
- Reinitialization must clear terminal/transient state, the function-convergence baseline, and adaptive-restart timing, then rerun algorithm initialization.
- Preserve evaluation counts across reinitialization. Existing tests show that this is intentional; document them as optimizer-lifetime counts.
- Preserve explicitly registered custom hooks rather than rebuilding an allegedly “clean” registry that silently removes them.
- Make uninitialized `mize_step()` fail immediately with an informative error.
- Change the stateful `max_iter` comparison to `>=`.
- Require cache-validity predicates everywhere current values are consumed.
- Fix best-result restoration and retain the `opt` returned by `mize_step_summary()` so counts and cache state remain consistent.
- Count every actual summary gradient evaluation, even when the value cannot safely be cached.

The stale-cache item is not merely theoretical: it currently causes `mize()` to return a worse point as its “best” result.

### 3. Establish hard evaluation budgets and safe line-search failure

The documentation and existing tests already promise conservative maximums, so hard caps are the correct interpretation.

- Centralize callback accounting sufficiently that `fn`, `gr`, combined `fg`, TN finite differences, convergence summaries, logging, initialization, and final reporting cannot invoke a callback with no remaining budget.
- TN must consider `max_fg`, not just `max_gr`, and increment counts only when `fg$gr` was actually called.
- Define the zero-budget behavior explicitly rather than returning an `f` obtained by violating the budget.
- Prefer genuine convergence or a detected nonfinite failure over “budget exhausted” when both are established by the final allowed evaluation.
- Bold Driver must never apply a nondecreasing, nonfinite, unevaluated, or over-budget candidate.
- Schmidt Armijo must return a zero step or a verified decreasing fallback when Armijo remains false at exhaustion.
- Apply the same safety rule to other line searches: condition-satisfying steps are accepted; a documented decreasing fallback may be accepted; otherwise return zero.
- Instrument callbacks independently in regression tests and compare actual invocations with `nf` and `ng`.

I reject the proposed public `accepted/reason/conditions_checked` result schema. A small internal success/failure indication is acceptable if needed to enforce safe behavior, but users have not presented a use case for a new public diagnostics contract.

### 4. Add narrow validation at boundaries mize directly consumes

Accept validation for:

- Scalarity, missingness, finiteness, and integer-valued count/cadence controls where applicable.
- Positive integer `memory`.
- Valid `preconditioner` and `tn_init` strings.
- A numeric `step0` for a constant line search.
- Nonempty finite numeric `par`.
- Scalar numeric objective results.
- Gradient vectors whose length exactly matches `par`.
- Combined `fg` shape.
- Dimensions, finiteness, and required symmetry of Hessian and inverse-Hessian
  results.

Nonfinite objective or gradient values encountered during optimization should continue to use the package’s existing failure paths where meaningful; validation should not destroy intentional line-search handling of `Inf`.

I reject public control constructors, rejection of irrelevant controls, and validation of values merely passed through to another package. No dependency behavior should become mize’s responsibility.

### 5. Use the standard quadratic inverse-BFGS update

- Replace the literal dense matrix–matrix evaluation with the standard
  matrix-vector/rank-two expansion recorded above.
- Preserve the existing curvature guard and `eps` denominator convention.
- Confirm or enforce the symmetric inverse-Hessian contract before using the
  one-`Hy` symmetric expansion.
- Add numerical-equivalence, symmetry, and public optimizer regressions.
- Record a one-time before/after benchmark at representative dimensions.
- Do not add a benchmark dependency, recurring performance CI, or permanent performance-reporting system.

This restores the expected \(O(n^2)\) dense-BFGS update without changing the API.

### 6. Small, bounded internal correctness cleanup

- Fix `cubic_extrapolate()`’s default branch.
- Replace `ras_ls(verbose = verbose)` with a valid default.
- Correct `safe_chol()` so its `eps` argument is honored, zero eigenvalues are repaired, and the symmetric eigensolver is used.
- Remove the unreachable private CG `ortho_check` branch rather than expanding the public CG API around it.
- Validate `tn_init` as zero or the documented `"previous"` value; remove the hidden `"l-bfgs"` initialization branch rather than documenting another control.
- Remove the false lifecycle comment claiming one-token events are supported. There is no current caller or public contract requiring that feature.

### 7. Correct documentation, attribution, and stale repository artifacts

The documentation work is a correction sweep, not an expansion:

- Fix NEWS’s `check_mize_gradient()` attribution and garbled status wording.
- Correct `abs_tol`, `max_fg`, `dbd_weight`, `nest_q`, gradient-tolerance, and count-lifetime descriptions.
- Update the stateful vignette’s obsolete iteration-passing text, `opt$error` reference, misspelled `mize_check_convergence`, duplicate `par_old`, and accidentally state-dependent example.
- Document the exact existing relative-tolerance formula rather than changing it.
- Regenerate Rd files from roxygen and render affected vignettes.
- Delete the four obsolete tracked files under `docs/`; pkgdown now builds from the vignette sources to `pkgdown-site/`.
- Remove nonexistent Rcpp paths from `.lintr`.
- `^plans$` was added to `.Rbuildignore` when the source audits and this plan
  became tracked, avoiding the package-check top-level-file note.

Apply the Schmidt and GPML notices and attribution recorded above. Leave the
More–Thuente-specific attribution and licensing edits to the final phase so
they match the permission-or-rewrite branch actually taken.

No SHAs, checksums, provenance manifests, or recurring verification process should be added.

### 8. Resolve More–Thuente last

Do not begin the licensing-driven replacement or More–Thuente attribution work
while phases 1–7 are in progress. Phase 3 may still apply its generic budget
and safe-failure contract at the line-search boundary, but it should not
restructure `R/more_thuente.R`. The permission request can proceed
asynchronously; this final phase has a definite two-branch decision rather than
an open-ended deferral.

If O'Leary supplies an explicit BSD-2-compatible grant covering the files as a
whole, or supplies evidence that completes the upstream MINPACK permission
chain:

- Retain the current implementation and public behavior.
- Add the exact copyright, licence, and attribution text required by the grant.
- Add O'Leary and any actual upstream code contributors to `Authors@R` in the
  appropriate contributor role; do not infer contributor roles merely from
  algorithm authorship.
- Run the dedicated More–Thuente tests, line-search condition tests, public
  optimizer regressions, and the full package suite.

If the response covers only O'Leary's contribution, is ambiguous about the
underlying code, declines permission, or no adequate response has arrived by
the time phases 1–7 finish:

- Treat Moré and Thuente (1994) as the implementation specification. Write new
  R code around the paper's interval, stage, acceptance, and safeguarded
  interpolation concepts; do not translate or paraphrase the MINPACK-2 or
  O'Leary source line by line.
- Preserve the exported behavior and existing controls: initial and maximum
  step bounds, step-growth limit, strong versus weak curvature, approximate
  Armijo, cubic safeguarding, and evaluation budgets. Clearly identify the
  weak-curvature, approximate-Armijo, and Xie–Schlick safeguards as `mize`
  extensions rather than features claimed from the 1994 algorithm.
- Use the paper's Tables 1–6 and their published test functions as the primary
  numerical oracle. Correct existing expectations where they deliberately
  reproduce O'Leary-specific differences from the paper instead of preserving
  those differences as compatibility requirements.
- Add property-level tests that every successful result satisfies the selected
  decrease and curvature conditions, remains within step bounds, respects hard
  callback budgets, and handles nonfinite evaluations and exhausted searches
  through the safe-failure contract established in phase 3.
- Retain public optimizer regressions, but do not require bit-for-bit internal
  trajectories, evaluation counts, private helper names, or legacy status
  codes when those are artifacts of the old translation rather than documented
  API behavior.
- Delete the translated comment corpus and source-shaped private helpers once
  replacement coverage is in place. Do not retain the old implementation as a
  fallback or second code path.
- Cite Moré and Thuente as algorithm authors rather than code contributors. Do
  not add the Netlib MINPACK notice to imply a licence chain the new code does
  not use.
- Run the full package suite and targeted comparisons against the published
  tables. No permanent cross-implementation harness or new dependency is
  warranted.

Completion of either branch is a release prerequisite. Reordering it last does
not permit another source or binary release containing code with an unresolved
licence chain.

## Rejected work

These recommendations do not solve a demonstrated current problem and are rejected:

- Rewrite the lifecycle/hook system into a new pipeline. It is unusual but coherent, covered by focused state-transition tests, and not itself causing the confirmed defects.
- Add exact HVP callbacks, arbitrary preconditioner callbacks, parameter scaling, damped BFGS, observer/user-stop callbacks, Hessian diagnostics, spectral gradients, or nonmonotone Armijo.
- Expose TN’s inner iteration limit or PHESS refresh interval.
- Remove or “graduate” PHESS. It is already public, documented, and behaviorally tested; neither audit showed incorrect results.
- Add `mize_reset()` or `reset_counts`. The existing lifetime-count behavior is defensible once documented.
- Add public convergence-control object constructors.
- Change the relative-function denominator to `max(1, ...)`. The current definition is strict but not shown to be wrong.
- Normalize every tolerance comparison to `<` or `<=`; exact equality semantics are not causing a demonstrated defect.
- Add configuration fuzzing, generalized randomized property campaigns, or tests merely to raise coverage in vendored recovery branches.
- Benchmark lifecycle overhead.
- Optimize `update_progress()` or the duplicate norm calculation in `normalize()`. Neither is shown to matter in a real workload.
- Purge parallel generations of line-search code wholesale. Remove a helper only when it is proved unused and the deletion is independently reviewable.
- Add method-comparison tutorials, broader documentation, or “publication-ready” material.

## Compact audit crosswalk

For the Sol audit:

- Accept 1, 2, 4, 5, 8, and the correctness portions of 3, 6, 7, 9, 11, 12, and 13.
- Reject the new reset API in 3; HVP and public TN controls in 6; public line-search diagnostics in 7; control constructors and over-validation in 9; lifecycle rewrite in 10; PHESS expansion/removal in 11; and the proposed relative-tolerance formula in 12.
- Reject all seven “missing features.”
- Accept only tests directly protecting accepted contracts, plus the one-time BFGS benchmark.

For the Fable audit:

- Accept A1–A4, although A4 is a narrow accounting invariant rather than a common public path.
- Accept B1, B2, and B4.
- Reject B3 as an unproven ambiguity; no current route demonstrates incorrect hook registration.
- From B5, accept the lost `opt` update and earlier validation as part of the relevant fixes; reject progress accumulation and `normalize()` as standalone work.
- Accept C1–C3.
- Accept targeted missing regressions, not coverage-driven testing.
- Accept the bounded licensing repair, stale `docs/` removal, and `.lintr` cleanup.
