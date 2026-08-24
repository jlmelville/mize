# Metric MDS

Metric multidimensional scaling (metric MDS) finds a low-dimensional
configuration whose pairwise distances resemble supplied
dissimilarities. It is also a useful example of an optimization problem
where a closure keeps fixed data out of the parameter vector.

We will reconstruct a two-dimensional configuration for the 21 European
cities in
[`datasets::eurodist`](https://rdrr.io/r/datasets/eurodist.html). The
input contains road distances in kilometres, so a two-dimensional
Euclidean configuration can only approximate them.

## Stress and its gradient

Let \\R = (r\_{ij})\\ be the target distances and let \\D = (d\_{ij})\\
be the Euclidean distances among candidate coordinates \\Y\\. We
minimize raw stress:

\\ C(Y) = \sum\_{i \< j} (r\_{ij} - d\_{ij})^2. \\

The condition \\i \< j\\ counts each unordered pair once. The matrix
implementation below instead sums both symmetric triangles, so it
multiplies that sum by one half. With \\n\\ objects, the corresponding
distance RMSE is \\\sqrt{C / {n \choose 2}}\\ and has the same units as
the input distances.

For distinct points, the derivative with respect to the coordinates of
object \\i\\ is

\\ \frac{\partial C}{\partial \mathbf{y}\_i} = -2 \sum\_{j \ne i}
\frac{r\_{ij} - d\_{ij}}{d\_{ij}} (\mathbf{y}\_i - \mathbf{y}\_j). \\

At coincident points, \\d\_{ij} = 0\\ and Euclidean distance is not
differentiable. The implementation uses a named `distance_epsilon` in
the denominator to keep the calculation finite. This numerical
convention slightly perturbs the displayed derivative. The derivative
remains undefined at a coincident configuration, so we validate the
implementation at a seeded, nondegenerate starting point.

The callback factory below captures the fixed target matrix in a
closure. Its `prepare()` helper performs the shared conversion from a
parameter vector to coordinates and pairwise distances. The counter will
let us compare shared work without relying on a noisy timing
measurement.

``` r

distance_epsilon <- 1e-10

make_mmds_callbacks <- function(
  distances,
  combined = FALSE,
  epsilon = distance_epsilon
) {
  target <- as.matrix(distances)
  preparation_count <- 0L

  prepare <- function(par) {
    preparation_count <<- preparation_count + 1L
    coordinates <- matrix(par, ncol = 2, byrow = TRUE)
    list(
      coordinates = coordinates,
      distances = as.matrix(stats::dist(coordinates))
    )
  }

  stress <- function(prepared) {
    0.5 * sum((target - prepared$distances)^2)
  }

  gradient <- function(prepared) {
    weights <- (target - prepared$distances) /
      (prepared$distances + epsilon)
    coordinates <- prepared$coordinates
    result <- matrix(
      nrow = nrow(coordinates),
      ncol = ncol(coordinates)
    )

    for (i in seq_len(nrow(coordinates))) {
      differences <- sweep(-coordinates, 2, -coordinates[i, ])
      result[i, ] <- colSums(differences * weights[, i])
    }

    as.vector(t(result)) * -2
  }

  callbacks <- list(
    fn = function(par) stress(prepare(par)),
    gr = function(par) gradient(prepare(par))
  )
  if (combined) {
    callbacks$fg <- function(par) {
      prepared <- prepare(par)
      list(
        fn = stress(prepared),
        gr = gradient(prepared)
      )
    }
  }

  list(
    fg = callbacks,
    preparations = function() preparation_count,
    reset_preparations = function() {
      preparation_count <<- 0L
      invisible(NULL)
    }
  )
}
```

The separate variant supplies `fn` and `gr`. The combined variant also
supplies `fg`, which lets `mize` request an objective and gradient from
one prepared configuration.

``` r

target_distances <- datasets::eurodist
set.seed(42)
initial_par <- rnorm(attr(target_distances, "Size") * 2)

separate_callbacks <- make_mmds_callbacks(target_distances)
combined_callbacks <- make_mmds_callbacks(
  target_distances,
  combined = TRUE
)
```

## Validate the gradient

[`check_mize_gradient()`](https://jlmelville.github.io/mize/reference/check_mize_gradient.md)
compares the analytic gradient with central finite differences. Because
the combined callback is present, it also checks that `fg` agrees with
the separate `fn` and `gr` callbacks.

``` r

gradient_check <- check_mize_gradient(
  combined_callbacks$fg,
  initial_par
)

gradient_summary <- data.frame(
  max_abs_error = gradient_check$max_abs_error,
  max_rel_error = gradient_check$max_rel_error,
  fg_fn_abs_error = gradient_check$fg_consistency$fn$abs_error,
  fg_gr_max_abs_error =
    gradient_check$fg_consistency$gr$max_abs_error
)
signif(gradient_summary, digits = 4)
#>   max_abs_error max_rel_error fg_fn_abs_error fg_gr_max_abs_error
#> 1      0.008302     1.748e-06               0                   0
```

The relative error is the more useful scale-free summary here. The exact
zero consistency errors show that the separate and combined paths
implement the same objective and gradient at the checked point.

## Optimize with separate and combined callbacks

This case study uses one explicit convergence criterion. Relative change
in stress is natural for this scale-dependent objective, so unrelated
tolerances are disabled. See the
[Convergence](https://jlmelville.github.io/mize/articles/convergence.md)
article when choosing a criterion for another problem.

``` r

mmds_controls <- list(
  method = "L-BFGS",
  max_iter = 200,
  abs_tol = NULL,
  rel_tol = 1e-8,
  grad_tol = NULL,
  step_tol = NULL
)

separate_callbacks$reset_preparations()
combined_callbacks$reset_preparations()

separate_result <- do.call(
  mize,
  c(
    list(par = initial_par, fg = separate_callbacks$fg),
    mmds_controls
  )
)
combined_result <- do.call(
  mize,
  c(
    list(par = initial_par, fg = combined_callbacks$fg),
    mmds_controls
  )
)

results <- list(
  separate = separate_result,
  combined = combined_result
)
initial_stress <- gradient_check$fn
pair_count <- choose(attr(target_distances, "Size"), 2)

termination_table <- data.frame(
  callbacks = names(results),
  status = vapply(results, `[[`, character(1), "status"),
  termination = vapply(
    results,
    function(result) result$terminate$what,
    character(1)
  ),
  iterations = vapply(results, `[[`, numeric(1), "iter"),
  nf = vapply(results, `[[`, numeric(1), "nf"),
  ng = vapply(results, `[[`, numeric(1), "ng")
)
termination_table
#>          callbacks    status termination iterations nf ng
#> separate  separate converged     rel_tol         37 57 57
#> combined  combined converged     rel_tol         37 57 57
```

``` r

fit_table <- data.frame(
  callbacks = names(results),
  initial_stress = initial_stress,
  final_stress = vapply(results, `[[`, numeric(1), "f"),
  initial_rmse_km = sqrt(initial_stress / pair_count),
  final_rmse_km = sqrt(
    vapply(results, `[[`, numeric(1), "f") / pair_count
  )
)
fit_display <- fit_table
fit_display[-1] <- lapply(fit_display[-1], signif, digits = 7)
fit_display
#>          callbacks initial_stress final_stress initial_rmse_km final_rmse_km
#> separate  separate      643185300      3356497        1750.082      126.4252
#> combined  combined      643185300      3356497        1750.082      126.4252
```

Here, “final stress” is stress at the returned parameters;
[`mize()`](https://jlmelville.github.io/mize/reference/mize.md) returns
the best observed parameters. Both variants report status `converged`
and terminate on `rel_tol`; they reduce the distance RMSE from about
1750 km to 126 km. The tables show the observed iteration and callback
counts without turning an exact count into a lasting promise.

The logical function and gradient counts are identical because both runs
ask for the same information. The underlying preparation counts expose
the work that `fg` shares:

``` r

preparation_counts <- c(
  separate = separate_callbacks$preparations(),
  combined = combined_callbacks$preparations()
)
data.frame(
  callbacks = names(preparation_counts),
  shared_preparations = unname(preparation_counts)
)
#>   callbacks shared_preparations
#> 1  separate                 114
#> 2  combined                  58
```

This deterministic count shows that a combined callback avoids repeated
work when the objective and gradient share an expensive intermediate.
The tiny example provides no useful wall-clock estimate.

## Interpret the configuration

Distances leave absolute position and orientation unidentified. Adding
one vector to every point leaves all distances unchanged, as do
rotations and reflections. Centering removes the arbitrary translation
for display. The plot shows distance geometry; its axes and handedness
carry no compass meaning.

``` r

center_coordinates <- function(par) {
  coordinates <- matrix(par, ncol = 2, byrow = TRUE)
  sweep(coordinates, 2, colMeans(coordinates))
}

plot_mmds <- function(coordinates, distances, ...) {
  horizontal_cutoff <- 0.08 * diff(range(coordinates[, 1]))
  label_positions <- ifelse(
    coordinates[, 1] < -horizontal_cutoff,
    2,
    ifelse(
      coordinates[, 1] > horizontal_cutoff,
      4,
      rep(c(1, 3), length.out = nrow(coordinates))
    )
  )
  label_positions[which.min(coordinates[, 2])] <- 3
  label_positions[which.max(coordinates[, 2])] <- 1
  city_names <- labels(distances)
  label_positions[city_names == "Geneva"] <- 2
  label_positions[city_names == "Lyons"] <- 4
  graphics::plot(
    coordinates,
    type = "n",
    asp = 1,
    xlab = "Centered dimension 1",
    ylab = "Centered dimension 2",
    main = "Metric MDS configuration of European road distances"
  )
  graphics::points(coordinates, pch = 19, cex = 0.45)
  graphics::text(
    coordinates[, 1],
    coordinates[, 2],
    labels = city_names,
    pos = label_positions,
    offset = 0.25,
    ...
  )
}

centered_coordinates <- center_coordinates(combined_result$par)
plot_mmds(centered_coordinates, target_distances, cex = 0.55)
```

![](mmds_files/figure-html/plot-results-1.png)

The invariance holds at the objective level. The hidden check below
translates, rotates, and reflects the returned configuration and
verifies that stress is unchanged to numerical tolerance.

The same closure pattern applies whenever callbacks need fixed data or
share an intermediate calculation. For the full callback contract, see
[`?mize`](https://jlmelville.github.io/mize/reference/mize.md).

## See also

- [Getting started](https://jlmelville.github.io/mize/articles/mize.md)
  introduces callback lists and gradient checks.
- [Convergence](https://jlmelville.github.io/mize/articles/convergence.md)
  explains the relative tolerance used here.
- [Stateful
  optimization](https://jlmelville.github.io/mize/articles/stateful.md)
  shows step-by-step integration with retained optimizer state.
