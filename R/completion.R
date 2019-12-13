#' @importFrom utils .DollarNames
#' @export
.DollarNames.rxEt <- function(x, pattern) {
  grep(pattern, .Call(`_RxODE_etDollarNames`, x), value = TRUE)
}

#' @export
.DollarNames.rxSolve <- function(x, pattern) {
  grep(pattern, .Call(`_RxODE_rxSolveDollarNames`, x), value = TRUE)
}

#' @export
str.rxSolve <- function(object, ..., nchar.max=128) {
  if (rxIs(object, "rxSolve")) {
    .dn <- .Call(`_RxODE_rxSolveDollarNames`, object)
    .max <- max(sapply(.dn, nchar))
    cat(sprintf("Classes 'rxSolve' and 'data.frame':\t%s rows of  %s variables:\n",
                object$.check.nrow, object$.check.ncol))
    if (any(names(object) == "sim.id")) {
      cat(paste0(" $ sim.id",
                 paste(rep(" ", .max - 6), collapse = ""),
                 ":"));
      str(object$sim.id)
    }
    if (any(names(object) == "id")) {
      cat(paste0(" $ id",
                 paste(rep(" ", .max - 2), collapse = ""),
                 ":"));
      str(object$id)
    }
    if (any(names(object) == "evid")) {
      cat(paste0(" $ evid",
                 paste(rep(" ", .max - 4), collapse = ""),
                 ":"));
      str(object$evid)
    }
    if (any(names(object) == "cmt")) {
      cat(paste0(" $ cmt",
                 paste(rep(" ", .max - 3), collapse = ""),
                 ":"));
      str(object$cmt)
    }
    if (any(names(object) == "ss")) {
      cat(paste0(" $ ss",
                 paste(rep(" ", .max - 2), collapse = ""),
                 ":"));
      str(object$ss)
    }
    if (any(names(object) == "amt")) {
      cat(paste0(" $ amt",
                 paste(rep(" ", .max - 3), collapse = ""),
                 ":"));
      str(object$amt)
    }
    if (any(names(object) == "rate")) {
      cat(paste0(" $ rate",
                 paste(rep(" ", .max - 4), collapse = ""),
                 ":"));
      str(object$rate)
    }
    if (any(names(object) == "dur")) {
      cat(paste0(" $ dur",
                 paste(rep(" ", .max - 3), collapse = ""),
                 ":"));
      str(object$dur)
    }
    if (any(names(object) == "ii")) {
      cat(paste0(" $ ii",
                 paste(rep(" ", .max - 2), collapse = ""),
                 ":"));
      str(object$ii)
    }
    cat(paste0(" $ time", paste(rep(" ", .max - 4), collapse = ""), ":"));
    str(object$time, nchar.max = nchar.max - .max - 8);
    .mv <- rxModelVars(object);
    if (length(.mv$lhs) > 0) {
      cat("Left Handed values ($lhs):\n")
      for (.l in .mv$lhs) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l)), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
    }
    if (length(.mv$state) > 0) {
      cat("State values ($state):\n")
      for (.l in .mv$state) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l)), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
      cat("State ini values:\n")
      for (.l in .mv$state) {
        cat(paste0(sprintf(" $ %s0", .l),
                   paste(rep(" ", .max - nchar(.l) - 1), collapse = ""),
                   ":"));
        str(object[[paste0(.l, "0")]],
            nchar.max = nchar.max - .max - nchar(.l) - 5);
      }
    }
    if (length(.mv$params) > 0) {
      cat("Parameter values ($params):\n")
      for (.l in .mv$params) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l) - 1), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
    }
    .vars <- c("sim.id", "id", "evid", "cmt", "ss", "amt", "rate",
               "dur", "ii", "time", .mv$lhs, .mv$state)
    .n2 <- names(object)[!(names(object) %in% .vars)]
    if (length(.n2) > 0) {
      cat("Other Variables:\n")
      for (.l in .n2) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l)), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
    }
    .dn <- .dn[!(.dn %in% c(names(object), paste0(.mv$state, "0"),
                            "t", "params", "inits",
                            "model", #fixme
                            .mv$params))]
    .fns <- sapply(.dn, function(x) {
      inherits(`$.rxSolve`(object, x), "function")
    })
    .fns <- names(.fns[.fns])
    if (length(.fns) > 0) {
      cat("Functions:\n")
      for (.l in .fns) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l)), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
    }
    .dn <- .dn[!(.dn %in% .fns)]
    if (length(.dn) > 0) {
      cat("Other:\n")
      for (.l in .dn) {
        cat(paste0(sprintf(" $ %s", .l),
                   paste(rep(" ", .max - nchar(.l)), collapse = ""),
                   ":"));
        str(object[[.l]], nchar.max = nchar.max - .max - nchar(.l) - 4);
      }
    }
  } else {
    NextMethod()
  }
}
