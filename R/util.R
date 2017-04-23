# partial function application
partial <- function(f, ...) {
  args <- list(...)
  function(...) {
    do.call(f, c(args, list(...)))
  }
}

# Square of the Euclidean norm of a vector
sqnorm2 <- function(v) {
  dot(v, v)
}

# l1 norm of a vector
norm1 <- function(v) {
  sum(abs(v))
}

# l2 (Euclidean) norm of a vector
norm2 <- function(v) {
  sqrt(dot(v, v))
}

# Infinity norm of a vector
norm_inf <- function(v) {
  max(abs(v))
}

# normalize a vector to length 1
normalize <- function(v) {
  l <- norm2(v)
  if (l < .Machine$double.eps) {
    v
  }
  else {
    v / norm2(v)
  }
}

# dot product of a and b
dot <- function(a, b) {
  sum(a * b)
}

clamp <- function(x, min_val = .Machine$double.eps, max_val = NULL) {
  x[x < min_val] <- min_val
  if (!is.null(max_val)) {
    x[x > max_val] <- max_val
  }
  x
}

sclamp <- function(x, min, max) {
  base::min(base::max(x, min), max)
}

vec_formatC <- function(v) {
  paste(Map(function(x) { formatC(x) }, v), collapse = ", ")
}

# convert a list to a strng
format_list <- function(ll) {
  Reduce(function(acc, elem) {
    paste0(acc,
           ifelse(nchar(acc) == 0, "", " "),
           elem,
           " = ",
           ifelse(length(ll[[elem]]) == 1,
                  formatC(ll[[elem]]), vec_formatC(ll[[elem]])))
  },
  names(ll), "")
}

# returns TRUE if x is in the range (left, right). By default, this is
# an open range, i.e. x == left and x == right is in the range.
# lopen if FALSE then range is [left, right) i.e. x = left is not in the range
# ropen if FALSE then range is (left, right] i.e. x = right is not in the range
is_in_range <- function(x, left, right, lopen = TRUE, ropen = TRUE) {
  `%lop%` <- ifelse(lopen, `<=`, `<`)
  `%rop%` <- ifelse(ropen, `<=`, `<`)

  left %lop% x && x %rop% right
}

# Checks if nullable x is finite
is_finite_numeric <- function(x) {
  is.numeric(x) && is.finite(x)
}

# Logging Hooks -----------------------------------------------------------


require_log_vals <- function(opt, stage, par, fg, iter) {
  message(iter, " ", substr(stage$type, 1, 2)
          ," par = ", vec_formatC(par)
          ," p = ", vec_formatC(stage$direction$value)
          , " a = ", formatC(stage$step_size$value)
          , " ap = ", vec_formatC(stage$result)
          , " f = ", formatC(fg$fn(par)))
  list(opt = opt)
}
attr(require_log_vals, 'name') <- 'log_vals'
attr(require_log_vals, 'event') <- 'after stage'

require_keep_stage_fs <- function(opt, stage, par, fg, iter) {
  if (is.null(opt$all_fs)) { opt$all_fs <- c() }
  f <- fg$fn(par)
  opt$all_fs <- c(opt$all_fs, f)
  list(opt = opt)
}
attr(require_keep_stage_fs, 'name') <- 'require_keep_stage_fs'
attr(require_keep_stage_fs, 'event') <- 'after stage'

