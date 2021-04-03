.ggplot2Fix <- function() {
  .ggplot2 <- loadNamespace("ggplot2")
  if (any(ls(.ggplot2) == "guide_none")) {
    assignInMyNamespace("guide_none", .ggplot2$guide_none)
  }
}
.hasUnits <- FALSE
.PreciseSumsVersion <- utils::packageVersion("PreciseSums")
.dparserVersion <- utils::packageVersion("dparser")
.onLoad <- function(libname, pkgname) { ## nocov start
  if (!identical(.dparserVersion, utils::packageVersion("dparser"))) {
    stop("RxODE compiled with dparser '", as.character(.dparserVersion),
         "' but dparser '", as.character(utils::packageVersion("dparser")),
         "' is loaded\nRecompile RxODE with the this version of dparser",
         call.=FALSE)
  }
  if (!identical(.PreciseSumsVersion, utils::packageVersion("PreciseSums"))) {
    stop("RxODE compiled with PreciseSums '", as.character(.PreciseSumsVersion),
         "' but PreciseSums '", as.character(utils::packageVersion("PreciseSums")),
         "' is loaded\nRecompile RxODE with the this version of PreciseSums",
         call.=FALSE)
  }
  if (requireNamespace("pillar", quietly = TRUE)) {
    .s3register("pillar::type_sum", "rxEvid")
    .s3register("pillar::type_sum", "rxRateDur")
    .s3register("pillar::pillar_shaft", "rxEvid")
    .s3register("pillar::pillar_shaft", "rxRateDur")
    .s3register("pillar::type_sum", "units")
    .s3register("pillar::type_sum", "mixed_units")
    .s3register("pillar::pillar_shaft", "units")
    .s3register("pillar::pillar_shaft", "mixed_units")
    .s3register("pillar::format_type_sum", "type_sum_units")
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    .s3register("tibble::as_tibble", "rxEt")
  }
  if (requireNamespace("data.table", quietly = TRUE)) {
    .s3register("data.table::as.data.table", "rxEt")
  }
  if (requireNamespace("units", quietly=TRUE)) {
    .s3register("units::set_units", "rxEt")
    .s3register("units::set_units", "rxRateDur")
    .s3register("units::drop_units", "rxEt")
    .s3register("units::drop_units", "rxSolve")
    .s3register("units::units<-","rxEvid")
    assignInMyNamespace(".hasUnits", TRUE)
    .units <- loadNamespace("units")
    assignInMyNamespace("type_sum.units", .units$type_sum.units)
    assignInMyNamespace("format_type_sum.type_sum_units", .units$format_type_sum.type_sum_units)
    assignInMyNamespace("pillar_shaft.units", .units$pillar_shaft.units)
    assignInMyNamespace("type_sum.mixed_units", .units$type_sum.mixed_units)
    assignInMyNamespace("pillar_shaft.mixed_units", .units$pillar_shaft.mixed_units)
  } else {
    assignInMyNamespace(".hasUnits", FALSE)
  }
  backports::import(pkgname)
  ## Setup RxODE.prefer.tbl
  .Call(`_RxODE_setRstudio`, Sys.getenv("RSTUDIO") == "1")
  rxSyncOptions("permissive")
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
  rxSyncOptions("permissive")
  if (!.rxWinRtoolsPath(retry = NA)) {
    ## nocov start
    packageStartupMessage("Rtools is not set up correctly!\n\nYou need a working Rtools installation for RxODE to compile models\n")
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
                        " threads (see ?getRxThreads)",
                        ifelse(.cacheIsTemp, "\n  no cache: create with `rxCreateCache()`", ""))
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

.mkCache <- function(.tmp) {
  if (!file.exists(.tmp)) {
    dir.create(.tmp, recursive = TRUE)
  } else if (!file.exists(file.path(.tmp, paste0(RxODE.md5, ".md5")))) {
    if (!.cacheIsTemp) packageStartupMessage("detected new version of RxODE, cleaning cache")
    unlink(.tmp, recursive=TRUE, force=TRUE)
    dir.create(.tmp, recursive = TRUE)
    writeLines("RxODE", file.path(.tmp, paste0(RxODE.md5, ".md5")))
  }
}

.cacheIsTemp <- TRUE
.rxTempDir0 <- NULL
.cacheDefault <- NULL
#' Get the RxODE temporary directory
#'
#' @return RxODE temporary directory.
#' @export
rxTempDir <- function() {
  if (is.null(getFromNamespace(".rxTempDir0", "RxODE"))) {
    .tmp <- Sys.getenv("rxTempDir")
    .rxUserDir <- R_user_dir("RxODE", "cache")
    assignInMyNamespace(".cacheIsTemp", FALSE)
    if (!file.exists(.rxUserDir)) {
      .rxUserDir <- file.path(tempdir(),"RxODE")
      assignInMyNamespace(".cacheIsTemp", TRUE)
    }
    if (.tmp == "") {
      if (is.null(.cacheDefault)) {
        assignInMyNamespace(".cacheDefault", .rxUserDir)
      }
      if (getOption("RxODE.cache.directory", .cacheDefault) != ".") {
        .tmp <- getOption("RxODE.cache.directory", .cacheDefault)
      } else {
        .tmp <- .rxUserDir
      }
    }
    .mkCache(.tmp)
    .tmp <- .normalizePath(.tmp)
    Sys.setenv(rxTempDir = .tmp)
    utils::assignInMyNamespace(".rxTempDir0", .tmp)
    utils::assignInMyNamespace("RxODE.cache.directory", .tmp)
    return(.tmp)
  } else {
    .tmp <- getFromNamespace(".rxTempDir0", "RxODE")
    .mkCache(.tmp)
    utils::assignInMyNamespace("RxODE.cache.directory", .tmp)
    return(.tmp)
  }
}
#' This will create the cache directory for RxODE to save between sessions
#'
#' When run, if the `R_user_dir` for RxODE's cache isn't present,
#' create the cache
#'
#' @return nothing
#'
#' @author Matthew Fidler
#'
#' @export
rxCreateCache <- function() {
  .tmp <- R_user_dir("RxODE", "cache")
  assignInMyNamespace(".cacheDefault", R_user_dir("RxODE", "cache"))
  .mkCache(.tmp)
  .tmp <- .normalizePath(.tmp)
  Sys.setenv(rxTempDir = .tmp)
  utils::assignInMyNamespace(".rxTempDir0", .tmp)
  utils::assignInMyNamespace("RxODE.cache.directory", .tmp)
  invisible()
}


#' Clear memoise cache for RxODE
#'
#' @author Matthew L. Fidler
#' @return nothing; called for side effects
#' @keywords internal
#' @export
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
#' Respect suppress messages
#'
#' This turns on the silent REprintf in C when `suppressMessages()` is
#' turned on. This makes the `REprintf` act like `messages` in R,
#' they can be suppressed with `suppressMessages()`
#'
#' @return Nothing
#' @author Matthew Fidler
#' @export
#' @examples
#'
#' # rxSupressMsg() is called with RxODE()
#'
#' # Note the errors are output to the console
#'
#' try(RxODE("d/dt(matt)=/3"),silent=TRUE)
#'
#' # When using suppressMessages, the output is suppressed
#'
#' suppressMessages(try(RxODE("d/dt(matt)=/3"),silent=TRUE))
#'
#' # In RxODE, we use REprintf so that interrupted threads do not crash R
#' # if there is a user interrupt. This isn't captured by R's messages, but
#' # This interface allows the `suppressMessages()` to suppress the C printing
#' # as well
#'
#' # If you  want to suppress messages from RxODE in other packages, you can use
#' # this function
rxSuppressMsg <- function() {
  if (requireNamespace("knitr", quietly = TRUE)) {
    if (!is.null(knitr::opts_knit$get('rmarkdown.pandoc.to'))) {
      return(invisible(NULL))
    } else {
      rxSetSilentErr(as.integer(length(capture.output(message(" "),type="message")) == 0L))
    }
  } else {
    rxSetSilentErr(as.integer(length(capture.output(message(" "),type="message")) == 0L))
  }
  invisible(NULL)
}

#' Sync options with RxODE variables
#'
#' Accessing RxODE options via getOption slows down solving.  This
#' allows the options to be synced with variables.
#'
#' @param setDefaults This will setup RxODE's default solving options with the following options:
#'
#' - `"none"` leave the options alone
#' - `"permissive"` This is a permissive option set similar to R language specifications.
#' - `"strict"` This is a strict option set similar to the original
#'    RxODE(). It requires semicolons at the end of lines and equals for
#'    assignment
#'
#' @author Matthew L. Fidler
#' @return nothing; called for side effects
#' @export
rxSyncOptions <- function(setDefaults=c("none", "permissive", "strict")) {
  x <- c("none" = 0L, "permissive" = 2L,
         "strict" = 1L)[match.arg(setDefaults)]
  if (x > 0) {
    op.rx <- list()
    for (v in names(rxOpt)) {
      op.rx[[v]] <- rxOpt[[v]][x]
    }
    options(op.rx) # nolint
  }
  for (var in names(rxOpt)) {
    assignInMyNamespace(var, getOption(var, rxOpt[[var]][1]))
  }
}
