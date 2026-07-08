# Convergence

There are a variety of ways that the optimization can terminate, running
the gamut from good (you have reached the minimum and further work is
pointless) to bad (the solution diverged and an infinity or NaN turned
up in a calculation).

We’ll use the 2D Rosenbrock function for the examples, which has a
minimum at `c(1, 1)`, where the function equals `0`.

``` r

rb_fg <- list(
   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
                          200 *        (x[2] - x[1] * x[1])) })
rb0 <- c(-1.2, 1)
```

## Iteration tolerance

An obvious way for the optimization to terminate is if you run out of
iterations:

``` r

res <- mize(rb0, rb_fg, max_iter = 10)
res$terminate
#> $what
#> [1] "max_iter"
#> 
#> $val
#> [1] 10
res$f
#> [1] 2.552598
res$par
#> [1] -0.5690157  0.2936478
```

When comparing different methods, the number of iterations is obviously
less important than the amount of actual CPU time you spent. Comparing
results with a fixed number of iterations is not a very good idea,
because different methods may do a lot more work within an iteration
than others. See the section on function and gradient tolerance below.

## Function tolerance

There are two ways to specify a function tolerance, based on comparing
the difference between consecutive function values. `abs_tol` measures
absolute tolerance.

``` r

res <- mize(rb0, rb_fg, max_iter = 100, abs_tol = 1e-8)
res$terminate
#> $what
#> [1] "abs_tol"
#> 
#> $val
#> [1] 4.061025e-10
res$f
#> [1] 4.606476e-11
res$par
#> [1] 0.9999967 0.9999928
```

However, relative tolerance is often preferred, because it measures the
change in value relative to the size of the values themselves.

``` r

res <- mize(rb0, rb_fg, max_iter = 100, rel_tol = 1e-3)
# hit relative tolerance
res$terminate
#> $what
#> [1] "abs_tol"
#> 
#> $val
#> [1] 1.06609e-08

# but stopped too early!
res$iter
#> [1] 37
res$f
#> [1] 4.521673e-10
res$par
#> [1] 0.9999972 0.9999923
```

In this example we stopped way too early. Even efficient methods like
L-BFGS may make little progress on some iterations, so don’t be too
aggressive with relative tolerance.

## Gradient tolerance

Gradient tolerances measure the difference between the size of the
gradient on consecutive step. `grad_tol` uses the 2-norm (sometimes
referred to as the Euclidean norm) of the gradient to measure
convergence.

``` r

res <- mize(rb0, rb_fg, abs_tol = 0, grad_tol = 1e-3)

res$terminate
#> $what
#> [1] "grad_tol"
#> 
#> $val
#> [1] 0.0009378802
res$f
#> [1] 4.521673e-10
res$par
#> [1] 0.9999972 0.9999923
```

This seems like a good stopping criterion because it is always zero at a
minimum, even if the function isn’t. It is also used to compare
different methods in Nocedal and Wright’s book. However, it has also
been recognized that it is not always reliable, see for instance this
[paper by Nocedal and
co-workers](https://doi.org/10.1023/A:1014897230089).

Other workers suggest using the infinity norm (the maximum absolute
component) of the gradient vector, particularly for larger problems. For
example, see this [conjugate gradient paper by Hager and
Zhang](https://doi.org/10.1137/030601880). To use the infinity norm, set
the `ginf_norm` parameter.

``` r

res <- mize(rb0, rb_fg, rel_tol = NULL, abs_tol = NULL, ginf_tol = 1e-3)

res$terminate
#> $what
#> [1] "ginf_tol"
#> 
#> $val
#> [1] 0.0008377527
res$f
#> [1] 4.521673e-10
res$par
#> [1] 0.9999972 0.9999923
```

While the gradient norms aren’t as reliable for checking convergence,
they almost never incur any overhead for checking, because the gradient
that’s calculated at the end of the iteration for this purpose can
nearly always be re-used for the gradient descent calculation at the
beginning of the next iteration, whereas the function-based convergence
requires the function to be calculated at the end of the iteration and
this is not always reused, although for many line search methods it is.

## Step tolerance

You can also look out for the change in `par` itself getting too small:

``` r

# set abs_tol to zero to stop it from triggering instead of step_tol
res <- mize(rb0, rb_fg, abs_tol = 0, step_tol = .Machine$double.eps)
res$terminate
#> $what
#> [1] "max_iter"
#> 
#> $val
#> [1] 100
res$iter
#> [1] 100
res$f
#> [1] 2.383932e-23
res$par
#> [1] 1 1
```

In most cases, the step tolerance should be a reasonable way to spot
convergence. Some optimization methods may allow for a step size of zero
for some iterations, preferring to commence the next iteration using the
same initial value of `par`, but with different optimization settings.
The step tolerance criterion knows when this sort of “restart” is being
attempted, and does not triggered under these conditions.

## Function and gradient count tolerance

For most problems, the time spent calculating the function and gradient
values will drown out any of the house-keeping that individual methods
do, so the number of function and gradient evaluations is the usual
determinant of how long an optimization takes. You can therefore decide
to terminate based on the number of function evaluations:

``` r

res <- mize(rb0, rb_fg, max_fn = 10)

res$terminate
#> $what
#> [1] "max_fn"
#> 
#> $val
#> [1] 10
res$nf
#> [1] 10
res$f
#> [1] 4.097612
res$par
#> [1] -1.015021  1.049581
```

Number of gradient evaluations:

``` r

res <- mize(rb0, rb_fg, max_gr = 10)

res$terminate
#> $what
#> [1] "max_gr"
#> 
#> $val
#> [1] 10
res$ng
#> [1] 10
res$f
#> [1] 4.097612
res$par
#> [1] -1.015021  1.049581
```

or both:

``` r

res <- mize(rb0, rb_fg, max_fg = 10)

res$terminate
#> $what
#> [1] "max_fg"
#> 
#> $val
#> [1] 10
res$nf
#> [1] 5
res$ng
#> [1] 5
res$f
#> [1] 4.280634
res$par
#> [1] -1.048507  1.070341
```

The function and gradient termination criteria are checked both between
iterations and during line search. On the assumption that if you specify
a maximum number of evaluations, that means these calculations are
expensive, `mize` errs on the side of caution and will sometimes
calculate fewer evaluations than you ask for, because it thinks that
attempting another iteration will exceed the limit.

## A minor complication with convergence checking

By default, convergence is checked at every iteration. For `abs_tol` and
`rel_tol`, this means that the function needs to have been evaluated at
the current value of `par`. A lot of optimization methods do this as
part of their normal working, so it doesn’t cost very much to do the
convergence check. However, not all optimization methods do. If you
specify a non-`NULL` value for `rel_tol` and `abs_tol` and the function
value isn’t available, it will be calculated. This could, for some
methods, add a lot of overhead.

If this is important, then using a gradient-based tolerance will be a
better choice.

`mize` internally uses the function value as a way to keep track of the
best `par` found during optimization. If this isn’t available, it will
use a gradient norm if that is being calculated. This is less reliable
than using function values, but better than nothing. If you turn off all
function and gradient tolerances then `mize` will be unable to return
the best set of parameters found over the course of the optimization.
Instead, you’ll get the last set of parameters it used.

If convergence checking at every iteration is too much of a burden, you
can reduce the frequency with which it is carried out with the
`check_conv_every` parameter:

``` r

res <- mize(rb0, rb_fg, grad_tol = 1e-3, check_conv_every = 5, verbose = TRUE)
#> 07:20:34 iter 0 f = 24.2 g2 = 232.9 nf = 1 ng = 1 step = 0 alpha = 0
#> 07:20:34 iter 5 f = 4.139 g2 = 1.773 nf = 7 ng = 7 step = 0.002565 alpha = 1
#> 07:20:34 iter 10 f = 2.553 g2 = 11.67 nf = 18 ng = 18 step = 0.04657 alpha = 0.05068
#> 07:20:34 iter 15 f = 1.37 g2 = 6.954 nf = 25 ng = 25 step = 0.0922 alpha = 0.405
#> 07:20:34 iter 20 f = 0.5142 g2 = 3.001 nf = 31 ng = 31 step = 0.1319 alpha = 1
#> 07:20:34 iter 25 f = 0.1203 g2 = 2.398 nf = 37 ng = 37 step = 0.03943 alpha = 0.9408
#> 07:20:34 iter 30 f = 0.009862 g2 = 3.333 nf = 42 ng = 42 step = 0.03554 alpha = 0.1706
#> 07:20:34 iter 35 f = 2.304e-06 g2 = 0.01537 nf = 47 ng = 47 step = 0.01136 alpha = 1
#> 07:20:34 iter 40 f = 4.147e-18 g2 = 4.386e-08 nf = 52 ng = 52 step = 6.386e-07 alpha = 1
```

This also has the side effect of producing less output to the console
when `verbose = TRUE`, because `log_every` is set to the same value of
`check_conv_every` by default. If you set them to different values,
`log_every` must be an integer multiple of `check_conv_every`. If it’s
not, it will be silently set to be equal to `check_conv_every`.

In many cases, however, convergence checking every iteration imposes no
overhead, so this is a non-issue. The vignette that runs through the
methods available in `mize` mentions where it might be an issue.
