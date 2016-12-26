# mize

Unconstrained Numerical Optimization Algorithms.

`mize` can be used as a standalone function like the `stats::optim` function, 
or can be integrated into other packages by creating a stateful optimizer and 
handling the iterations, convergence, logging and so on externally. 

`mize` knows how to do Broyden-Fletcher-Goldfarb-Shanno (BFGS), 
the limited-memory BFGS (L-BFGS), various flavors of Conjugate Gradient (CG), 
Nesterov Accelerated Gradient (NAG) and momentum-based methods, among others.

## Installing

```R
# install.packages("devtools")
devtools::install_github("jlmelville/mize")
```

## Documentation

```R
?mize
```

A vignette is on its way.

## Examples

```R
# Make a list containing the function and gradient:
rosenbrock_fg <- list(
   fn = function(x) { 100 * (x[2] - x[1] * x[1]) ^ 2 + (1 - x[1]) ^ 2  },
   gr = function(x) { c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
                          200 *        (x[2] - x[1] * x[1])) })
# Starting point:
rb0 <- c(-1.2, 1)

# Minimize using L-BFGS
res <- mize(rb0, rosenbrock_fg, method = "L-BFGS")
# Optimized parameters are in res$par

# Or create an optimizer and then loop manually
opt <- make_mize(method = "L-BFGS")
opt <- mize_init(opt, rb0, rosebrock_fg)

par <- rb0
done <- FALSE
iter <- 0
while (!done) {
  iter <- iter + 1
  res <- mize_step(opt, par, rosenbrock_fg, iter)
  par <- res$par
  opt <- res$opt
  # Look at res$f for current function value
  # you get to (i.e. have to) decide when to stop
}
```

## See also

The Wolfe line search functionality uses a conversion of Carl Edward Rasmussen's
[Matlab code](http://learning.eng.cam.ac.uk/carl/code/minimize/) and Dianne 
O'Leary's Matlab translation of the 
[More'-Thuente line search](https://www.cs.umd.edu/users/oleary/software/)
algorithm from [MINPACK](http://www.netlib.org/minpack/).

## License

[BSD 2-Clause](https://opensource.org/licenses/BSD-2-Clause).
