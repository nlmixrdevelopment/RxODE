##'@export
print.rxEtTran <- function(x, ...) {
    print(as.data.frame(x));
    .cls <- class(x);
    .lst <- attr(.cls, ".RxODE.lst")
    cat("\nCovariates (non time-varying):\n")
    print(.lst$cov1)
    cat("\nCompartment translation:\n");
    print(data.frame("Compartment Name" = .lst$cmtInfo,
                     "Compartment Number" = seq_along(.lst$cmtInfo),
                     check.names = FALSE))
}

##'@export
print.rxHidden <- function(x,...) {
    cat("\r");
}

##'@rdname rxEvid
##' @export
print.rxEvid <- function(x, ...) {
    cat(paste(.colorFmt.rxEvid(x),collapse="\n"),"\n")
    return(invisible(x))
}

##'@export
print.rxRateDur <- function(x, ...) {
    cat(paste(.colorFmt.rxRateDur(x),collapse="\n"),"\n")
    return(invisible(x))
}

##'@export
print.rxEt <- function(x,...) {
  if (rxIs(x, "rxEt")) {
    bound <- .getBound(x, parent.frame(2));
    .et1 <- paste0("EventTable with ",x$nobs+x$ndose, " records")
    .et2 <- NULL
    .units <- x$.units;
    .maxId <- length(x$IDs)
    if (.maxId !=1) {
      .et2 <- sprintf("   %s individuals", .maxId)
    }
    .et3 <-sprintf("   %s dosing records (see %s$%s(); add with %s or %s)",
                   x$ndose, bound, "get.dosing", "add.dosing", "et")
    .et4 <- sprintf("   %s observation times (see %s$%s(); add with %s or %s)",
                    x$nobs, bound, "get.sampling", "add.sampling", 
                    "et")
    .et5 <- NULL
    if (x$show["addl"]) {
      .et5 <- sprintf("   multiple doses in `addl` columns, expand with %s$%s(); or %s(%s)",
                  bound, "expand", "etExpand", bound)
    }
    .et <- c(.et2, .et3, .et4, .et5)
    .df <- data.frame(et=.et,stringsAsFactors=FALSE)
    names(.df) <- .et1
    class(.df) <- c(sprintf("EventTable Info: %s", bound),
                    "paged_df","data.frame")
    .out <- utils::capture.output({print(.df)})
    .nb <- TRUE
    if (length(.out) > 0){
      .nb <- FALSE
      cat(cli::cli_format_method({
        cli::cli_rule(center=crayon::bold(.et1))
      }), sep="\n")
      cat(paste0(.et2,"\n"))
      cat(sprintf("   %s dosing records (see %s$%s(); add with %s or %s)\n",
                  x$ndose, crayon::yellow(bound), crayon::blue("get.dosing"),
                  crayon::blue("add.dosing"), crayon::blue("et")))
      cat(sprintf("   %s observation times (see %s$%s(); add with %s or %s)\n",
                  x$nobs, crayon::yellow(bound), crayon::blue("get.sampling"),
                  crayon::blue("add.sampling"), crayon::blue("et")))
      if (x$show["addl"]) {
        cat(sprintf("   multiple doses in `addl` columns, expand with %s$%s(); or %s(%s)\n",
                    crayon::yellow(bound), crayon::blue("expand"),
                    crayon::blue("etExpand"), crayon::yellow(bound)))
      }
    }
    if (x$nobs!=0 | x$ndose!=0) {
      if (!.nb) {
        cat(cli::cli_format_method({
          cli::cli_rule(crayon::bold(paste0("First part of ",crayon::yellow(bound),":")))
        }), sep="\n")
      }
      print(tibble::as_tibble(data.frame(.etAddCls(x))));
    }
    invisible(x)
  } else {
    print.data.frame(x)
  }
}

##' Print information about the RxODE object.
##'
##' This prints the model name and its status for being able to be solved
##'
##' @param x An rxode object
##' @param ... Ignored parameters
##' @author Matthew L.Fidler
##' @export
print.RxODE <- function(x, ...){
    rxModelVars(x);
    x  <- .getReal(x);
    .bound <- .getBound(x, parent.frame(2));
    .valid <- x$isValid()
    .msg2 <- "";
    if (!.valid){
        .msg <- crayon::red$bold("invalid");
        .msg2 <- paste0(" re-create with ", crayon::blue("RxODE::"), crayon::yellow("RxODE"))
        .ico <- crayon::red(cli::symbol$cross)
    } else {
        .loaded <- x$isLoaded();
        if (.loaded){
            .msg <- crayon::green$bold("ready")
            .ico <- crayon::green(cli::symbol$tick)
        } else{
            .msg <- crayon::yellow$bold("unloaded");
            .ico <- crayon::yellow(cli::symbol$warning);
            .msg2 <- paste0(" reload with ", crayon::blue("RxODE::"), crayon::yellow("rxLoad"));
        }
    }
    if (.useUtf()){
        .ico <- paste0(.ico, " ");
    } else {
        .ico <- ""
    }
    .dll <- getOption("RxODE.basename.print",basename(RxODE::rxDll(x)));
    .env <- attr(x, ".env")
    .pkg  <- "";
    .new  <- ""
    if (!is.null(x$package)){
        .dll <- .bound
        .pkg  <- crayon::blue(paste0(x$package, "::"))
        if (regexpr("_new",.dll)!=-1){
            .dll  <- gsub("_new","",.dll);
            .new  <- " (updated)";
        }
    } else {
        .dll <- substr(.dll, 1, nchar(.dll) - nchar(.Platform$dynlib.ext) - 1)
    }
    cat(paste0(crayon::bold("RxODE "), as.vector(RxODE::rxVersion()["version"]), " model named ",.pkg,
               crayon::yellow$bold(.dll), .new, " model (", .ico, .msg,
               .msg2, ")."), "\n")
    .indLin <- rxModelVars(x)$indLin;
    if (length(.indLin) > 0){
        cat(crayon::bold("  indLin: "));
        if (.indLin[["fullIndLin"]]){
            if (is.null(.indLin[["f"]])){
                cat(paste0(crayon::yellow("homogenous matrix exponential + Inductive Linearization")))
            } else {
                cat(paste0(crayon::red("in"), crayon::yellow("homogenous matrix exponential + Inductive Linearization")))
            }
        } else {
            if (is.null(.indLin[["f"]])){
                cat(crayon::yellow("homogenous matrix exponential"));
            } else {
                cat(paste0(crayon::red("in"), crayon::yellow("homogenous matrix exponential")))
            }
        }
        cat("\n");
    }
    if (!any(names(list(...)) == "rxSuppress") && .valid){
        .cur <- RxODE::rxState(x);
        if (length(.cur) > 0)
            cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$state"), ": ", paste(.cur, collapse=", "), "\n"))
        .cur <- x$stateExtra
        if (length(.cur) > 0)
            cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$stateExtra"), ": ", paste(.cur, collapse=", "), "\n"))
        .cur <- RxODE::rxParams(x);
        if (length(.cur) > 0)
            cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$params"), ": ", paste(.cur, collapse=", "), "\n"))
        .cur <- RxODE::rxLhs(x);
        if (length(.cur) > 0)
            cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$lhs"), ": ", paste(.cur, collapse=", "), "\n"))
    }
    invisible(x)
}

##'@export
print.rxModelVars <- function(x, ...) {
    .bound <- .getBound(x, parent.frame(2));
    cat("RxODE model variables (see str to see all variables)");
    .cur <- x$state;
    if (length(.cur) > 0)
        cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$state"), ": ", paste(.cur, collapse=", "), "\n"))
    .cur <- x$stateExtra;
    if (length(.cur) > 0)
        cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$stateExtra"), ": ", paste(.cur, collapse=", "), "\n"))
    .cur <- x$params;
    if (length(.cur) > 0)
        cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$params"), ": ", paste(.cur, collapse=", "), "\n"))
    .cur <- x$lhs;
    if (length(.cur) > 0)
        cat(paste0(crayon::yellow(.bound), crayon::blue$bold("$lhs"), ": ", paste(.cur, collapse=", "), "\n"))
    invisible(x)
}

##' Print the rxCoef object
##'
##' This prints out the user supplied arguments for rxCoef object
##'
##' @param x rxCoef object
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
print.rxCoef <- function(x, ...){
  .rxDllObj <- x$RxODE;
  if (length(rxParams(.rxDllObj)) > 0){
    cat(cli::cli_format_method({
      cli::cli_rule(left="User supplied parameters:")
    }), "\n");
    print(RxODE::rxInits(.rxDllObj, c(), RxODE::rxParams(.rxDllObj), NA, TRUE))
    cat(cli::cli_format_method({
      cli::cli_rule(left="User initial conditions:")
    }), "\n")
    .tmp <- RxODE::rxInits(.rxDllObj, c(), RxODE::rxState(.rxDllObj), 0, TRUE);
    if (length(x$sens) > 0){
      .tmp <- .tmp[regexpr(getFromNamespace("regSens", "RxODE"), names(.tmp)) == -1];
    }
    .tmp <- .tmp[!(names(.tmp) %in% x$fn.ini)];
    print(.tmp);
  }
  if (length(x$fn.ini) > 0){
    cat(cli::cli_format_method({
      cli::cli_rule(left="Parameter-based initial conditions:")
    }), "\n")
    print(x$fn.ini);
  }
  cat(cli::cli_format_method({
    cli::cli_rule(left="Compartments:")
  }), "\n")
  .tmp <- RxODE::rxState(.rxDllObj);
  if (length(.tmp) > 0){
    names(.tmp) <- paste0("cmt=", seq_along(.tmp));
    if (length(x$sens) > 0){
      .tmp1 <- .tmp[regexpr(getFromNamespace("regSens", "RxODE"), .tmp) == -1];
      print(.tmp1);
      cat(cli::cli_format_method({
        cli::cli_rule(left="Sensitivities:")
      }), "\n")
      .tmp2 <- gsub(getFromNamespace("regSens", "RxODE"), "d/dt(d(\\1)/d(\\2))",
                    .tmp[regexpr(regSens, .tmp) != -1]);
      print(.tmp2);
    } else {
      print(.tmp);
    }
  } else {
    cat("No ODEs in this DLL.\n");
  }
  return(invisible());
}

##' Print the rxCoefSolve object
##'
##' This prints out the user supplied arguments for the rxCoef object
##'
##' @param x rxCoefSolve object
##' @param ... Other (ignored) parameters.
##' @keywords Internal
##' @author Matthew L.Fidler
##' @export
print.rxCoefSolve <- function(x, ...){
    cat("\nUser supplied parameters ($params):\n");
    print(x$params);
    cat("\nUser initial conditions ($state):\n");
    print(x$state);
    rxDllObj <- x$RxODE;
    if (length(RxODE::rxInits(rxDllObj)) > 0){
        cat("\nDefault parameter values:\n")
        print(RxODE::rxInits(rxDllObj));
    }
    cat("\nCompartments:\n");
    tmp <- RxODE::rxState(rxDllObj);
    names(tmp) <- paste0("cmt=", seq_along(tmp));
    sens <- tmp[regexpr(getFromNamespace("regSens", "RxODE"), tmp) != -1];
    if (length(sens) > 0){
        tmp <- tmp[regexpr(getFromNamespace("regSens", "RxODE"), tmp) == -1]
    }
    print(tmp);
    if (length(sens) > 0){
        sens <- gsub(getFromNamespace("regSens", "RxODE"), "d/dt(d(\\1)/d(\\2))", sens);
        cat("\nSensitivities:\n");
        print(sens)
    }
    return(invisible());
}

##' @export
print.rxC <- function(x, ...){
    cat(sprintf("C file: %s  ('summary' for code)\n", x));
}


##' Print rxDll object
##'
##' This tells if the rxDll is loaded, ready and/or deleted.
##'
##' @keywords internal
##' @author Matthew L.Fidler
##' @export
print.rxDll <- function(x, ...){
    if (file.exists(x$dll)){
        cat(sprintf("RxODE DLL named \"%s\"", getOption("RxODE.basename.print", basename(x$dll))));
        if (rxDllLoaded(x)){
            cat(" is loaded and ready to use.\n");
        } else {
            cat(" is not loaded now.\n");
        }
    } else {
        cat(sprintf("RxODE DLL named \"%s\" has been deleted.\n", getOption("RxODE.basename.print", basename(x$dll))));
    }
    invisible(x);
}

##'@export
print.rxSymInvChol <- function(x, ...){
    d <- dim(x$fmat)[1]
    cat(sprintf("Object to create Omega and Omega^-1 & derivitaves for a %sx%s matrix:\n", d, d))
    print(x$fmat);
    cat("Use `rxSymInvChol' for the matrix.\n");
}

##' @export
print.boundCovs <- function(x, ...){
  cat(format(x, ...), sep="\n");
}


##' @export
print.rxSolveCovs <- function(x, ...){
  .args <- as.list(match.call(expand.dots = TRUE));
  if (any(names(.args) == "bound")) {
    .bound <- .args$bound;
  } else {
    .bound <- .getBound(x, parent.frame(2));
  }
  .df <- x$covs;
  if (!is.null(.df)) {
    if (rxIs(.df, "data.frame")){
      print(structure(.bound, class="boundCovs"), ...)
      print(tibble::as_tibble(.df), ...);
    }
  }
  NextMethod();
}

##'@export
print.boundInits <- function(x, ...){
  cat(format(x, ...), sep="\n");
}

##' @export
print.rxSolveInits <- function(x, ...){
  .args <- as.list(match.call(expand.dots = TRUE));
  if (any(names(.args) == "bound")) {
    .bound <- .args$bound;
  } else {
    .bound <- .getBound(x, parent.frame(2));
  }
  print(structure(.bound, class="boundInits"), ...)
  .df <- x$inits;
  print(.df)
  NextMethod();
}

##'@export
print.rxSolveSimType <- function(x, ...){
  if (any(names(x) == "sim.id")){
    cat(format(x, ...), sep="\n")
  }
}


##'@export
print.rxSolve <- function(x, ...){
  if (rxIs(x, "rxSolve")) {
    .nb <- TRUE
    .args <- as.list(match.call(expand.dots = TRUE));
    if (any(names(.args) == "n")) {
      .n <- .args$n;
    } else {
      .n <- 6L;
    }
    if (any(names(.args) == "width")) {
      .width <- .args$width;
    } else {
      .width <- NULL;
    }
    if (any(names(.args) == "bound")) {
      .bound <- .args$bound;
    } else {
      .bound <- .getBound(x, parent.frame(2));
    }
    if (.nb){
      .df <- x$pars
      if (rxIs(.df, "data.frame")) {
        .cls <- c(paste0("Parameters ", .bound, "$params"),
                "paged_df", "data.frame");
        class(.df) <- .cls
        .out <- utils::capture.output({print(.df)})
        if (length(.out) > 0) .nb <- FALSE
      }
    }
    if (.nb){
      .df <- x$covs;
      if (!is.null(.df)) {
        if (rxIs(.df, "data.frame")){
          .cls <- c(paste0("Covariates ", .bound, "$covs"),
                    "paged_df", "data.frame");
          class(.df) <- .cls;
          .out <- utils::capture.output({print(.df)})
          if (length(.out) > 0) .nb <- FALSE
        }
      }
    }
    if (.nb) {
      .df <- x$inits;
      .df <- as.data.frame(t(x$inits))
      .cls <- c(paste0("Initial\u00A0State ", .bound, "$inits"),
                "paged_df", "data.frame");
      class(.df) <- .cls
      .out <- utils::capture.output({print(.df)})
      if (length(.out) > 0) .nb <- FALSE
    }
    if (.nb){
      print.rxSolveSimType(x);
      .df <- x;
      .cls <- c(paste0("Solved\u00A0Data: ", .bound),
                "paged_df", "data.frame");
      class(.df) <- .cls
      print(.df)
      return(invisible(x));
    } else {
      .summary <- any(names(.args) == ".summary")
      if (!.summary){
        cat(cli::cli_format_method({
          d <- cli::cli_div(theme = list(rule = list(
            "line-type" = "bar2")))
          cli::cli_rule(center=crayon::bold("Solved RxODE object"))
          cli::cli_end(d);
        }), sep="\n")
      }
      NextMethod()
      if (.summary) {
        cat(cli::cli_format_method({
          cli::cli_rule(left=crayon::bold("Summary of data (object):"))
        }), sep="\n")
        print(summary.data.frame(x))
        cat(cli::cli_format_method({
          d <- cli::cli_div(theme = list(rule = list(
            "line-type" = "bar2")))
          cli::cli_rule()
          cli::cli_end(d)
        }), sep="\n")
      } else {
        cat(cli::cli_format_method({
          cli::cli_rule(left=crayon::bold("First part of data (object):"))
        }), sep="\n")
        .isDplyr <- requireNamespace("tibble", quietly = TRUE) && RxODE.display.tbl;
        if (!.isDplyr) {
          print(head(as.matrix(x), n = .n));
        } else {
          print(tibble::as_tibble(x), n = .n, width = .width);
        }
        cat(cli::cli_format_method({
          d <- cli::cli_div(theme = list(rule = list(
            "line-type" = "bar2")))
          cli::cli_rule()
          cli::cli_end(d)
        }), sep="\n")
      }
    }
  } else {
    print.data.frame(x)
  }
}

##'@export
print.rxModelCode <- function(x, ...){
  htmltools::code(x)
}


##'@export
print.rxModelText <- function(x, ...) {
  .args <- as.list(match.call(expand.dots = TRUE));
  .summary <- any(names(.args) == ".summary")
  if (any(names(.args) == "bound")) {
    .bound <- .args$bound;
  } else {
    .bound <- .getBound(x, parent.frame(2));
  }
  .code <- deparse(body(eval(parse(text=paste("function() {",as.vector(x),"}")))))
  .code[1]  <- "RxODE({"
  .code[length(.code)]  <- "})";
  if (.summary) {
    cat(cli::cli_format_method({
      cli::cli_rule(left=.fmt3("Model", .bound, "model"));
    }), sep="\n")
  } else {
    cat(cli::cli_format_method({
      d <- cli::cli_div(theme = list(rule = list(
        "line-type" = "bar2")))
      cli::cli_rule(center=crayon::bold("RxODE Model Syntax"))
      cli::cli_end(d);
    }), sep="\n")
  }
  cat(paste(.code,collapse="\n"), "\n");
  if (!.summary) {
    cat(cli::cli_format_method({
      d <- cli::cli_div(theme = list(rule = list(
        "line-type" = "bar2")))
      cli::cli_rule()
      cli::cli_end(d);
    }), sep="\n")
  }
}

##' @export
print.boundParams <- function(x, ...){
  cat(format(x, ...), sep="\n");
}

##' @export
print.rxSolveParams <- function(x, ...){
  .args <- as.list(match.call(expand.dots = TRUE));
  if (any(names(.args) == "bound")) {
    .bound <- .args$bound;
  } else {
    .bound <- .getBound(x, parent.frame(2));
  }
  class(.bound) <- "boundParams"
  .df <- x$.params.single
  if (length(.df) > 0) {
      print(.bound, ...)
      print(.df, ...)
  } else {
    .df <- x$pars
    if (rxIs(.df, "data.frame")) {
      print(.bound, ...)
      print(tibble::as_tibble(.df), ...);
    }
  }
  NextMethod();
}

##'@export
print.rxSymInvCholEnv <- function(x, ...) {
    if (is.null(x$theta)) {
        cat(sprintf("Uninitialized $theta, please assign (requires %s arguments)!\n", x$ntheta))
    } else {
        cat(sprintf("$theta=c(%s) for:\n\n", paste(x$theta, collapse=", ")))
        print(x$invobj$fmat)
        cat("\nThis allows accessing $omegaInv, $omega, etc. For a full list see str(.)\n");
    }
}
