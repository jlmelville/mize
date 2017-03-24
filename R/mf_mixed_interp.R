mixedInterp <-
  function(bracket, bracketFval, bracketDval,
           Tpos,
           oldLOval, oldLOFval, oldLODval,
           debug) {
    # Mixed Case
    nonTpos <- -Tpos + 3

    gtdT <- bracketDval[Tpos]
    gtdNonT <- bracketDval[nonTpos]
    oldLOgtd <- oldLODval
    if (bracketFval[Tpos] > oldLOFval) {
      alpha_c <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)))
      alpha_q <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, NA)))
      if (abs(alpha_c - oldLOval) < abs(alpha_q - oldLOval)) {
        if (debug) {
          message('Cubic Interpolation')
        }
        res <- alpha_c
      }
      else {
        if (debug) {
          message('Mixed Quad/Cubic Interpolation')
        }
        res <- (alpha_q + alpha_c) / 2
      }
    }
    else if (dot(gtdT, oldLOgtd) < 0) {
      alpha_c <- polyinterp(point_matrix(
          c(oldLOval, bracket[Tpos]),
          c(oldLOFval, bracketFval[Tpos]),
          c(oldLOgtd, gtdT)))
      alpha_s <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, NA),
        c(oldLOgtd, gtdT)))
      if (abs(alpha_c - bracket[Tpos]) >= abs(alpha_s - bracket[Tpos])) {
        if (debug) {
          message('Cubic Interpolation')
        }
        res <- alpha_c
      }
      else {
        if (debug) {
          message('Quad Interpolation')
        }
        res <- alpha_s
      }
    }
    else if (abs(gtdT) <= abs(oldLOgtd)) {
      alpha_c <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(oldLOFval, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)), min(bracket), max(bracket))
      alpha_s <- polyinterp(point_matrix(
        c(oldLOval, bracket[Tpos]),
        c(NA, bracketFval[Tpos]),
        c(oldLOgtd, gtdT)), min(bracket), max(bracket))

      if (alpha_c > min(bracket) && alpha_c < max(bracket)) {
        if (abs(alpha_c - bracket[Tpos]) < abs(alpha_s - bracket[Tpos])) {
          if (debug) {
            message('Bounded Cubic Extrapolation')
          }
          res <- alpha_c
        }
        else {
          if (debug) {
            message('Bounded Secant Extrapolation')
          }
          res <- alpha_s
        }
      }
      else {
        if (debug) {
          message('Bounded Secant Extrapolation')
        }
        res <- alpha_s
      }

      if (bracket[Tpos] > oldLOval) {
        res <- min(bracket[Tpos] + 0.66 * (bracket[nonTpos] - bracket[Tpos]),
                   res)
      }
      else {
        res <- max(bracket[Tpos] + 0.66 * (bracket[nonTpos] - bracket[Tpos]),
                   res)
      }
    }
    else {
      res <- polyinterp(point_matrix(
        c(bracket[nonTpos], bracket[Tpos]),
        c(bracketFval[nonTpos], bracketFval[Tpos]),
        c(gtdNonT, gtdT)))
    }
    res
  }
