partial <- function(f, ...) {
  args <- list(...)
  function(...) {
    do.call(f, c(args, list(...)))
  }
}

norm2 <- function(v) {
  sqrt(sum(v * v))
}

normalize <- function(v) {
  l <- norm2(v)
  if (l < .Machine$double.eps) {
    v
  }
  else {
    v / norm2(v)
  }
}

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

lreplace <- function(l, ...) {
  varargs <- list(...)
  for (i in names(varargs)) {
    l[[i]] <- varargs[[i]]
  }
  l
}

vec_formatC <- function(v) {
  paste(Map(function(x) { formatC(x) }, v), collapse = ", ")
}

format_vec <- function(vec) {
  paste(formatC(vec), collapse = ' ')
}


