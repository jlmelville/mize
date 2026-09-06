test_that("scalar inverse Hessians match one by one matrices", {
  for (method in c("BFGS", "SR1")) {
    for (inverse in c(0.5, 2, 3)) {
      fg <- list(fn = function(x) x^2, gr = function(x) 2 * x)
      fg$hi <- function(x) inverse
      scalar <- mize(1, fg, method = method, max_iter = 2)
      fg$hi <- function(x) matrix(inverse, 1, 1)
      full <- mize(1, fg, method = method, max_iter = 2)
      expect_equal(scalar, full, info = paste(method, inverse))
    }
  }
})
