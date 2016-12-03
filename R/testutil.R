# Functions used only for testing

# Create Initial Step Value
#
# Given a set of start parameters and a search direction, initializes the
# step data. Utility function for testing.
make_step0 <- function(fn, gr, x, pv, f = fn(x), df = gr(x)) {
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


