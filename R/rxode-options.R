.ggplot2Fix <- function() {
  .ggplot2 <- loadNamespace("ggplot2")
  if (any(ls(.ggplot2) == "guide_none")) {
    assignInMyNamespace("guide_none", .ggplot2$guide_none)
  }
}
.onLoad <- function(libname, pkgname) { ## nocov start
  .s3register("pillar::type_sum", "rxEvid")
  .s3register("pillar::type_sum", "rxRateDur")
  .s3register("pillar::pillar_shaft", "rxEvid")
  .s3register("pillar::pillar_shaft", "rxRateDur")
  .s3register("pillar::type_sum", "units")
  .s3register("pillar::type_sum", "mixed_units")
  .s3register("pillar::pillar_shaft", "units")
  .s3register("pillar::pillar_shaft", "mixed_units")
  .s3register("pillar::format_type_sum", "type_sum_units")
  .s3register("tibble::as_tibble", "rxEt")
  .s3register("data.table::as.data.table", "rxEt")

  ## .s3Register("tibble::as_tibble", "rxEt")
  ## .s3Register("as.data.table", "rxEt")

  backports::import(pkgname)
  ## Setup RxODE.prefer.tbl
  .Call(`_RxODE_setRstudio`, Sys.getenv("RSTUDIO") == "1")
  rxPermissive(respect = TRUE) ## need to call respect on the first time
  suppressMessages(.rxWinRtoolsPath(retry = NA))
  rxTempDir()
  if (!interactive()) {
    setProgSupported(0)
  }
  .getDTEnv()
  .ggplot2Fix()
} ## nocov end

.onAttach <- function(libname, pkgname) {
  ## For some strange reason, mvnfast needs to be loaded before RxODE to work correctly
  .Call(`_RxODE_setRstudio`, Sys.getenv("RSTUDIO") == "1")
  rxPermissive(respect = TRUE) ## need to call respect on the first time
  if (!.rxWinRtoolsPath(retry = NA)) {
    ## nocov start
    packageStartupMessage("Rtools is not set up correctly!\n\nYou need a working Rtools installation for RxODE to work.\nYou can set up Rtools using the command 'rxWinSetup()'.\n")
    ## nocov end
  }
  if (!interactive()) {
    setProgSupported(0)
  }
  rxTempDir()
  .getDTEnv()
  .ggplot2Fix()
  v <- utils::packageVersion("RxODE")
  packageStartupMessage("RxODE ", v, " using ", getRxThreads(verbose=FALSE),
                        " threads (see ?getRxThreads)")
  if (!.Call(`_rxHasOpenMp`)) {
    packageStartupMessage("========================================\n",
        "RxODE has not detected OpenMP support and will run in single-threaded mode\n",
        if (Sys.info()["sysname"]=="Darwin")
          "This is a Mac. Please read https://mac.r-project.org/openmp/"
        else
          paste0("The system is ", Sys.info()["sysname"], "; To get best performance enable OpenMP"),
        "\n========================================\n")
  }
}

.onUnload <- function(libpath) {
  ## nocov start
  rxUnloadAll()
  gc() # Force garbage collection finalization
  library.dynam.unload("RxODE", libpath)
  ## nocov end
}

.rxTempDir0 <- NULL
.cacheDefault <- NULL
##' Get the RxODE temporary directory
##'
##' @return RxODE temporary directory.
##' @export
rxTempDir <- function() {
  if (is.null(getFromNamespace(".rxTempDir0", "RxODE"))) {
    .tmp <- Sys.getenv("rxTempDir")
    if (.tmp == "") {
      if (is.null(.cacheDefault)) {
        assignInMyNamespace(".cacheDefault", R_user_dir("RxODE", "cache"))
      }
      if (getOption("RxODE.cache.directory", .cacheDefault) != ".") {
        .tmp <- getOption("RxODE.cache.directory", .cacheDefault)
      } else {
        .tmp <- R_user_dir("RxODE", "cache")
      }
    }
    if (!file.exists(.tmp)) {
      dir.create(.tmp, recursive = TRUE)
    }
    .tmp <- .normalizePath(.tmp)
    Sys.setenv(rxTempDir = .tmp)
    utils::assignInMyNamespace(".rxTempDir0", .tmp)
    utils::assignInMyNamespace("RxODE.cache.directory", .tmp)
    return(.tmp)
  } else {
    .tmp <- getFromNamespace(".rxTempDir0", "RxODE")
    if (!file.exists(.tmp)) {
      dir.create(.tmp, recursive = TRUE)
    }
    utils::assignInMyNamespace("RxODE.cache.directory", .tmp)
    return(.tmp)
  }
}


##' Clear memoise cache for RxODE
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxForget <- function() {
  for (fn in ls(envir = getNamespace("RxODE"))) {
    if (memoise::is.memoised(getFromNamespace(fn, "RxODE"))) {
      memoise::forget(getFromNamespace(fn, "RxODE"))
    }
  }
}

## strict/permissive
rxOpt <- list(
  RxODE.prefer.tbl = c(FALSE, FALSE),
  RxODE.warn.on.assign = c(TRUE, TRUE),
  RxODE.syntax.assign = c(FALSE, TRUE),
  RxODE.syntax.star.pow = c(FALSE, TRUE),
  RxODE.syntax.require.semicolon = c(TRUE, FALSE),
  RxODE.syntax.allow.dots = c(FALSE, TRUE),
  RxODE.syntax.allow.ini0 = c(FALSE, TRUE),
  RxODE.syntax.allow.ini = c(FALSE, TRUE),
  RxODE.calculate.jacobian = c(FALSE, FALSE),
  RxODE.calculate.sensitivity = c(FALSE, FALSE),
  RxODE.verbose = c(TRUE, TRUE),
  RxODE.suppress.syntax.info = c(FALSE, FALSE),
  RxODE.sympy.engine = c("", ""),
  RxODE.cache.directory = c(.cacheDefault, .cacheDefault),
  RxODE.syntax.assign.state = c(FALSE, FALSE),
  RxODE.tempfiles = c(TRUE, TRUE),
  RxODE.sympy.run.internal = c(FALSE, FALSE),
  RxODE.syntax.require.ode.first = c(TRUE, TRUE),
  RxODE.compile.O = c("3", "3"),
  RxODE.unload.unused = c(FALSE, FALSE)
)
RxODE.prefer.tbl <- NULL
RxODE.warn.on.assign <- NULL
RxODE.syntax.assign <- NULL
RxODE.syntax.star.pow <- NULL
RxODE.syntax.require.semicolon <- NULL
RxODE.syntax.allow.dots <- NULL
RxODE.syntax.allow.ini0 <- NULL
RxODE.syntax.allow.ini <- NULL
RxODE.calculate.jacobian <- NULL
RxODE.calculate.sensitivity <- NULL
RxODE.verbose <- NULL
RxODE.suppress.syntax.info <- NULL
RxODE.sympy.engine <- NULL
RxODE.cache.directory <- NULL
RxODE.delete.unnamed <- NULL
RxODE.syntax.assign.state <- NULL
RxODE.tempfiles <- NULL
RxODE.sympy.run.internal <- NULL
RxODE.syntax.require.ode.first <- NULL
RxODE.compile.O <- NULL
RxODE.unload.unused <- NULL

.isTestthat <- function() {
  return(regexpr("/tests/testthat/", getwd(), fixed = TRUE) != -1) # nolint
}

##' Permissive or Strict RxODE syntax options
##'
##' This sets the RxODE syntax to be permissive or strict
##'
##' @param expr Expression to evaluate in the permissive/strict
##'     environment.  If unspecified, set the options for the current
##'     environment.
##' @param respect when TRUE, respect any options that are specified.
##'     This is called at startup, but really should not be called
##'     elsewhere, otherwise the options are not changed.
##' @param silent when true, also silence the syntax errors and
##'     interactive output (useful in testing).
##' @param test When specified as a string, the enclosed test is
##'     skipped unless the environmental variable "NOT_CRAN" equals this
##'     value. Special values are "cran" which means test on CRAN
##' @author Matthew L. Fidler
##' @export
rxPermissive <- function(expr, silent = .isTestthat(),
                         respect = FALSE,
                         test = "cran") {
  args <- as.list(match.call())[-1]
  args$op.rx <- 2
  do.call(getFromNamespace("rxOptions", "RxODE"), args, envir = parent.frame(1))
}
##' @rdname rxPermissive
##' @export
rxStrict <- function(expr, silent = .isTestthat(), respect = FALSE) {
  ## nocov start
  args <- as.list(match.call())[-1]
  args$op.rx <- 1
  do.call(getFromNamespace("rxOptions", "RxODE"), args, envir = parent.frame(1))
  ## nocov end
}
##' Options for RxODE
##'
##' This is a backend for \code{rxPermissive} (with
##' \code{op.rx} = \code{2}) and \code{rxStrict} (with
##' \code{op.rx} =\code{1})
##'
##' When \code{expr} is missing and \code{op.rx} is NULL, this
##' displays the current RxODE options.
##'
##' @inheritParams rxPermissive
##' @param op.rx A numeric for strict (1) or permissive (2) syntax.
##' @author Matthew L. Fidler
##' @export
rxOptions <- function(expr, op.rx = NULL, silent = .isTestthat(), respect = FALSE,
                      test = "cran") {
  if (!respect) rxUnloadAll()
  rxSetSilentErr(1L)
  do.it <- TRUE
  .test <- .test0 <- Sys.getenv("NOT_CRAN")
  if (Sys.getenv("rxCran") != "") {
    .test <- .test0 <- Sys.getenv("rxCran")
  }
  ## (identical to cran testing)
  on.exit({
    rxSetSilentErr(0L)
    if (!respect) rxUnloadAll()
  })
  if (any(.test == c("false", "", "cran"))) {
    ## devtools sets NOT_CRAN=false for revdep testing
    ## I don't believe CRAN sets this at all.
    ## Sys.unsetenv("NOT_CRAN")
    if (any(test == c("false", "", "cran"))) {
      do.it <- TRUE
    } else {
      do.it <- FALSE
    }
  } else {
    ## devtools sets NOT_CRAN=true for testing everything
    if (test == .test || .test == "true") {
      do.it <- TRUE
    } else {
      do.it <- FALSE
    }
  }
  if (do.it) {
    .lastCran <- Sys.getenv("NOT_CRAN")
    Sys.setenv("NOT_CRAN" = "true")
    on.exit(
      {
        Sys.setenv("NOT_CRAN" = .lastCran)
      },
      add = TRUE
    )
    if (missing(expr) && is.null(op.rx)) {
      op <- options()
      op <- op[regexpr(rex::rex("RxODE."), names(op)) != -1]
      op <- op[order(names(op))]
      sapply(names(op), function(n) {
        message(sprintf("%s: %s", n, op[[n]]))
      })
      return(invisible(op))
    } else {
      if (is(op.rx, "character")) {
        if (op.rx == "strict") {
          op.rx <- 1
        } else {
          op.rx <- 2
        }
      }
      if (is(op.rx, "numeric")) {
        if (op.rx <= 2) {
          x <- op.rx
          op.rx <- list()
          for (v in names(rxOpt)) {
            op.rx[[v]] <- rxOpt[[v]][x]
          }
        }
      }
      if (!missing(silent)) {
        op.rx$RxODE.verbose <- !silent
        op.rx$RxODE.suppress.syntax.info <- silent
      }
      if (!missing(expr)) {
        opOld <- options()
        .oldProg <- getProgSupported()
        if (silent) {
          setProgSupported(-1)
        }
        on.exit(
          {
            options(opOld)
            setProgSupported(.oldProg)
            rxSyncOptions()
          },
          add = TRUE
        )
      }
      if (respect) {
        op <- options()
        w <- !(names(op.rx) %in% names(op))
        if (any(w)) options(op.rx[w])
        rxSyncOptions()
      } else if (length(op.rx) > 0) {
        options(op.rx)
        rxSyncOptions()
      }
      if (is(substitute(expr), "{")) {
        if (silent) {
          return(suppressMessages(eval(substitute(expr), envir = parent.frame(1))))
        } else {
          return(eval(substitute(expr), envir = parent.frame(1)))
        }
      }
    }
  }
}

##' Sync options with RxODE variables
##'
##' Accessing RxODE options via getOption slows down solving.  This
##' allows the options to be synced with variables.
##'
##' @author Matthew L. Fidler
##' @export
rxSyncOptions <- function() {
  for (var in names(rxOpt)) {
    assignInMyNamespace(var, getOption(var, rxOpt[[var]][1]))
  }
}

rxSkipValidate <- function() {
  if (!identical(Sys.getenv("RxODE_VALIDATION_FULL"), "true")) {
    testthat::skip("Only run on full validation")
  }
}

rxSyncOptions()
