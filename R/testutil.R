# Functions used only for testing


# Rosenbrock ---------------------------------------------------------------

out0 <- c(-1.2, 1)
# taken from the optim man page
rosenbrock_fg <- list(
  fn = function(x) {
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1) ^ 2 + (1 - x1) ^ 2
  },
  gr = function(x) {
    x1 <- x[1]
    x2 <- x[2]
    c(
      -400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
      200 *      (x2 - x1 * x1))
  },
  hs = function(x) {
    xx <- 1200 * x[1] * x[1] - 400 * x[2] + 2
    xy <- x[1] * -400
    yy <- 200
    matrix(c(xx, xy, xy, yy), nrow = 2)
  },
  fg <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    a <- (x2 - x1 * x1)
    b <- 1 - x1
    list(
      f = 100 * a * a + b * b,
      g = c(
        -400 * x1 * a - 2 * b,
        200 * a
      )
    )
  },
  n = 2
)

tricky_fg <- function() {
  res <- list(
    fn = function(par) {
      x1 <- par[1]
      x2 <- par[2]
      (1.5 - x1 + x1 * x2)^2 +
        (2.25 - x1 + x1 * x2 * x2)^2 +
        (2.625 - x1 + x1 * x2 * x2 * x2)^2
    }
  )
  res$gr <- make_gfd(res$fn)
  res$hs <- make_hfd(res$fn)

  res
}

wrap_fg <- function(fg) {
  nc <- fg$n
  nr <- 1

  list(
    nc = nc,
    nr = nr,
    method = list(
      cost = list(
        fn = function(inp, out, method) {
          fg$fn(out$ym)
        }
      ),
      eps = .Machine$double.eps
    ),
    grad_func = function(inp, out, method, mat_name) {
      list(gm = matrix(fg$gr(out[[mat_name]]), nrow = nr, ncol = nc))
    }
  )
}
rosenbrock <- wrap_fg(rosenbrock_fg)


# Create Initial Step Value
#
# Given a set of start parameters and a search direction, initializes the
# step data. Utility function for testing.
make_step0 <- function(fg, x, pv, f = fg$fn(x), df = fg$gr(x)) {
  list(
    x = x,
    alpha = 0,
    f = f,
    df = df,
    d = dot(pv, df)
  )
}

# More'-Thuente test functions --------------------------------------------


# Test Function 1 ---------------------------------------------------------

# 1 phi(a) = -(a) / (a^2 + b) phi'(a) = (a^2 - b) / (a^2 + b)^2 b = 2
fn1 <- function(alpha, beta = 2) {
  -alpha / (alpha ^ 2 + beta)
}
gr1 <- function(alpha, beta = 2) {
  (alpha ^ 2 - beta) / ((alpha ^ 2 + beta) ^ 2)
}
fcn1 <- function(x) {
  list(f = fn1(x), g = gr1(x))
}


# Test Function 2 ---------------------------------------------------------

# 2 phi(a) = (a + b)^5 - 2(a + b)^4  phi'(a) = (a+b)^3*(5a+5b-8)  b = 0.004
fn2 <- function(alpha, beta = 0.004) {
  (alpha + beta) ^ 5 - 2 * (alpha + beta) ^ 4
}
gr2 <- function(alpha, beta = 0.004) {
  (alpha + beta) ^ 3 * (5 * alpha + 5 * beta - 8)
}
fcn2 <- function(x) {
  list(f = fn2(x), g = gr2(x))
}


# Test Function 3 ---------------------------------------------------------

# 3 phi(a) = phi_0(a) + [2(1-b)/(l*pi)]sin(l*pi*alpha*0.5)
#   phi'(a) =  phi_0'(a) + (1-b)cos(l*pi*alpha*0.5)
#   phi_0(a) = 1 - a if a <= 1 - b  phi_0'(a) = -1
#            = a - 1 if a >= 1 + b  phi_0'(a) = 1
#            = (a-1)^2/2b + b/2 if a in [1-b,1+b] phi_0'(a) =  a - 1 /b
#     b = 0.01 l = 39
fn3 <- function(alpha, beta = 0.01, l = 39) {
  if (alpha <= 1 - beta) {
    phi_0 <- 1 - alpha
  } else if (alpha >= 1 + beta) {
    phi_0 <- alpha - 1
  } else {
    phi_0 <- (alpha - 1) ^ 2 / (2 * beta) + (beta / 2)
  }
  phi_0 + (2 * (1 - beta) * sin(l * pi * alpha * 0.5)) / (l * pi)
}
#   phi'(a) =  phi_0'(a) + (1-b)cos(l*pi*alpha*0.5)

gr3 <- function(alpha, beta = 0.01, l = 39) {
  if (alpha <= 1 - beta) {
    dphi_0 <- -1
  } else if (alpha >= 1 + beta) {
    dphi_0 <- 1
  } else {
    dphi_0 <- (alpha - 1) / beta
  }

  dphi_0 + (1 - beta) * cos(l * pi * alpha * 0.5)
}
fcn3 <- function(x) {
  list(f = fn3(x), g = gr3(x))
}


# Utils for Test Functions 4-6 --------------------------------------------

# gamma(b) = (1+b^2)^1/2 - b
yanaig <- function(beta) {
  sqrt(1 + beta ^ 2) - beta
}

yanai1 <- function(alpha, beta) {
  sqrt((1 - alpha) ^ 2 + beta ^ 2)
}
gryanai1 <- function(alpha, beta) {
  (alpha - 1) / yanai1(alpha, beta)
}

yanai2 <- function(alpha, beta) {
  sqrt(alpha ^ 2 + beta ^ 2)
}
gryanai2 <- function(alpha, beta) {
  alpha / yanai2(alpha, beta)
}

# phi(a) = gamma(b_1)[(1-a)^2 + b_2^2]^1/2 + gamma(b_2)[a^2 + b_1^2]^1/2
yanai <- function(alpha, beta1, beta2) {
  yanaig(beta1) * yanai1(alpha, beta2) +
    yanaig(beta2) * yanai2(alpha, beta1)
}

#   phi'(a) = -[gamma(b_1)(1-a)]/sqrt[(1-a)^2 + b_2^2] + gamma(b_2)a/sqrt([a^2 + b_1^2])
gryanai <- function(alpha, beta1, beta2) {
  (yanaig(beta1) * gryanai1(alpha, beta2)) + (yanaig(beta2) * gryanai2(alpha, beta1))
}

# Test Function 4 ---------------------------------------------------------

fn4 <- function(alpha) {
  yanai(alpha, beta1 = 0.001, beta2 = 0.001)
}
gr4 <- function(alpha) {
  gryanai(alpha, beta1 = 0.001, beta2 = 0.001)
}
fcn4 <- function(x) {
  list(f = fn4(x), g = gr4(x))
}


# Test Function 5 ---------------------------------------------------------

fn5 <- function(alpha) {
  yanai(alpha, beta1 = 0.01, beta2 = 0.001)
}
gr5 <- function(alpha) {
  gryanai(alpha, beta1 = 0.01, beta2 = 0.001)
}
fcn5 <- function(x) {
  list(f = fn5(x), g = gr5(x))
}


# Test Function 6 ---------------------------------------------------------


fn6 <- function(alpha) {
  yanai(alpha, beta1 = 0.001, beta2 = 0.01)
}
gr6 <- function(alpha) {
  gryanai(alpha, beta1 = 0.001, beta2 = 0.01)
}
fcn6 <- function(x) {
  list(f = fn6(x), g = gr6(x))
}

f1 <- list(fn = fn1, gr = gr1)
f2 <- list(fn = fn2, gr = gr2)
f3 <- list(fn = fn3, gr = gr3)
f4 <- list(fn = fn4, gr = gr4)
f5 <- list(fn = fn5, gr = gr5)
f6 <- list(fn = fn6, gr = gr6)


# Euro Cities -------------------------------------------------------------

make_mmds_fg <- function(dist_mat) {
  dxm <- as.matrix(dist_mat)

  par_to_ym <- function(par) {
    matrix(par, ncol = 2, byrow = TRUE)
  }

  par_to_dym <- function(par) {
    as.matrix(dist(par_to_ym(par)))
  }

  f <- function(par) {
    dym <- par_to_dym(par)
    diff2 <- (dxm - dym) ^ 2
    sum(diff2)
  }

  g <- function(par) {
    dym <- par_to_dym(par)
    km <- (dxm - dym) / (dym + 1.e-10)

    ym <- par_to_ym(par)
    gm <- matrix(nrow = nrow(ym), ncol = ncol(ym))
    for (i in 1:nrow(ym)) {
      dyij <- sweep(-ym, 2, -ym[i, ])
      gm[i, ] <- apply(dyij * km[, i], 2, sum)
    }

    - 4 * as.vector(t(gm))
  }

  h <- make_hfd(f)

  list(fn = f,
       gr = g,
       hs = h)
}

eurodist_fg <- function() { make_mmds_fg(eurodist) }
ed0 <- as.vector(t(-cmdscale(eurodist, add = TRUE)$points))

us_fg <- function() { make_mmds_fg(UScitiesD)}
us0 <- as.vector(t(-cmdscale(UScitiesD, add = TRUE)$points))


# Rotates qm onto pm
kabsch <- function(pm, qm) {
  # center the points
  pm <- scale(qm, center = TRUE, scale = FALSE)
  qm <- scale(qm, center = TRUE, scale = FALSE)

  am <- t(pm) %*% qm

  svd_res <- svd(am)
  # use the sign of the determinant to ensure a right-hand coordinate system
  d <- determinant(svd_res$v %*% t(svd_res$u))$sign
  dm <- diag(c(1, d))

  # rotation matrix
  um <- svd_res$u %*% dm %*% t(svd_res$v)

  t(um %*% t(qm)) + attr(pm, "scaled:center")
}

plot_mmds <- function(coords, dist, ...) {
  if (class(coords) == "numeric") {
    coords <- matrix(coords, ncol = 2, byrow = TRUE)
  }
  plot(coords, type = 'n')
  text(coords[, 1], coords[, 2], labels = labels(dist), ...)
}

rotate_mmds_results <- function(par, dist) {
  pm <- -cmdscale(dist, add = TRUE)$points
  qm <- matrix(par, ncol = 2, byrow = TRUE)

  kabsch(pm, qm)
}

# Finite Difference -------------------------------------------------------

gfd <- function(par, fn, eps =  1.e-3) {
  g <- rep(0, length(par))
  for (i in 1:length(par)) {
    oldx <- par[i]

    par[i] <- oldx + eps
    fplus <- fn(par)

    par[i] <- oldx - eps
    fminus <- fn(par)
    par[i] <- oldx

    g[i] <- (fplus - fminus) / (2 * eps)
  }
  g
}

make_gfd <- function(fn, eps = 1.e-3) {
  function(par) {
    gfd(par, fn, eps)
  }
}

hfd <- function(par, fn, eps =  1.e-3) {
  hs <- matrix(0, nrow = length(par), ncol = length(par))
  for (i in 1:length(par)) {
    for (j in i:length(par)) {
      if (i != j) {
        oldxi <- par[i]
        oldxj <- par[j]

        par[i] <- par[i] + eps
        par[j] <- par[j] + eps
        fpp <- fn(par)

        par[j] <- oldxj - eps
        fpm <- fn(par)

        par[i] <- oldxi - eps
        par[j] <- oldxj + eps
        fmp <- fn(par)

        par[j] <- oldxj - eps
        fmm <- fn(par)

        par[i] <- oldxi
        par[j] <- oldxj

        val <- (fpp - fpm - fmp + fmm) / (4 * eps * eps)

        hs[i, j] <- val
        hs[j, i] <- val
      }
      else {
        f <- fn(par)
        oldxi <- par[i]

        par[i] <- oldxi + 2 * eps
        fpp <- fn(par)

        par[i] <- oldxi + eps
        fp <- fn(par)

        par[i] <- oldxi - 2 * eps
        fmm <- fn(par)

        par[i] <- oldxi - eps
        fm <- fn(par)

        par[i] <- oldxi

        hs[i, i] <- (-fpp + 16 * fp - 30 * f + 16 * fm - fmm) / (12 * eps * eps)
      }
    }
  }
  hs
}

make_hfd <- function(fn, eps = 1.e-3) {
  function(par) {
    hfd(par, fn, eps)
  }
}



