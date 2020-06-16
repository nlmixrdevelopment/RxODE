.addExpr <- function(.ret) {
  .new <- .rxOptEnv$.rep[[.ret]];
  if (!is.null(.new)) {
    if (length(.rxOptEnv$.exclude) != 1) .rxOptEnv$.exclude <- ""
    if (.new == .rxOptEnv$.exclude) {
      return(.ret)
    } else {
      .rxOptEnv$.new <- c(.rxOptEnv$.new, .new);
      return(.new)
    }
  } else {
    if (is.null(.rxOptEnv$.list[[.ret]])) {
      .rxOptEnv$.list[[.ret]] <- 1L;
    } else {
      .rxOptEnv$.list[[.ret]] <- .rxOptEnv$.list[[.ret]] + 1L;
    }
  }
  return(.ret)
}

.rxOptFn <- function(fn) {
  force(fn);
  function(...) {
    .ret <- paste0(fn, "(", paste(unlist(list(...)), collapse = ", "), ")")
    return(.addExpr(.ret))
  }
}
.rxOptBin <- function(sep) {
  force(sep)
  function(e1, e2) {
    .add <- TRUE
    if (missing(e2)) {
      if (sep == "+") {
        .ret <- paste0(e1)
      } else {
        .ret <- paste0(gsub(" ", "", sep), e1)
      }
      if (regexpr(rex::rex(start, any_spaces, regNum, any_spaces, end),
                  .ret, perl=TRUE) != -1) {
        .add <- FALSE
      }
    } else {
      if (sep == "^" && isTRUE(checkmate::checkIntegerish(suppressWarnings(as.numeric(e2)), lower=2,
                                                          any.missing=FALSE))) {
        .ret <- paste0("(", paste(rep(paste0("(", e1, ")"), as.numeric(e2)), collapse="*"), ")")
      } else {
        if ((regexpr(rex::rex(start, any_spaces, regNum, any_spaces, end),
                     paste0(e1), perl=TRUE) != -1) &&
              (regexpr(rex::rex(start, any_spaces, regNum, any_spaces, end),
                       paste0(e2), perl=TRUE) != -1)) {
          .add <- FALSE
        }
        .ret <- paste0(e1, sep, e2)
      }
    }
    if (.add) {
      return(.addExpr(.ret))
    } else {
      return(.ret)
    }
  }
}

.rxOptEnv <- new.env(parent = emptyenv())
.rxOptEnv[["^"]] <- .rxOptBin("^")
.rxOptEnv[["**"]] <- .rxOptBin("^")

.rxOptEnv[["*"]] <- .rxOptBin("*");

.rxOptEnv[["/"]] <- .rxOptBin("/");

.rxOptEnv[["+"]] <- .rxOptBin("+");

.rxOptEnv[["-"]] <- .rxOptBin("-");

.rxOptEnv[["&&"]] <- .rxOptBin("&&");
.rxOptEnv[["||"]] <- .rxOptBin("||");
.rxOptEnv[["|"]] <- .rxOptBin("|");
.rxOptEnv[["&"]] <- .rxOptBin("&");
.rxOptEnv[["=="]] <- .rxOptBin("==");
.rxOptEnv[["<="]] <- .rxOptBin("<=");
.rxOptEnv[[">="]] <- .rxOptBin(">=");
.rxOptEnv[["<"]] <- .rxOptBin("<");
.rxOptEnv[[">"]] <- .rxOptBin(">");
.rxOptEnv[["!="]] <- .rxOptBin("!=");

.rxOptEnv[["["]] <- function(name, val) {
  .n <- toupper(name)
  .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
  if (any(.n == c("THETA", "ETA")) && is.numeric(val)) {
    if (round(val) == val && val > 0) {
      return(sprintf("%s[%s]", .n, val));
    } else {
      stop(.err);
    }
  } else {
    stop(.err)
  }
}
.rxOptEnv[["{"]] <- function(...) {
  return(sprintf("{\n%s\n}", paste(unlist(list(...)), collapse = "\n")))
}
.rxOptEnv[["["]] <- function(name, val) {
  .n <- toupper(name)
  .err <- "RxODE only supports THETA[#] and ETA[#] numbers."
  if (any(.n == c("THETA", "ETA")) && is.numeric(val)) {
    if (round(val) == val && val > 0) {
      return(sprintf("%s[%s]", .n, val));
    } else {
      stop(err);
    }
  } else {
    stop(err)
  }
}

.rxOptEnv[["("]] <- unaryOp("(", ")")

.rxOptEnv$.list <- list();
.rxOptEnv$.rep <- list();
.rxOptEnv$.exclude <- "";
.rxOptEnv$.new <- c();
.rxOptEnv$.added <- c();

.rxOptGetEnv <- function(expr) {
  ## Known functions
  .calls <- allCalls(expr)
  .callList <- setNames(lapply(.calls, .rxOptFn), .calls)
  .callEnv <- list2env(.callList);
  .currEnv <- cloneEnv(.rxOptEnv, .callEnv);
  .names <- allNames(expr)
  .n1 <- .names;
  .n2 <- .names;
  .symbolList <- setNames(as.list(.n2), .n1);
  .symbolEnv <- list2env(.symbolList, parent = .currEnv);
  return(.symbolEnv)
}

.rxOptExpr <- function(x) {
  .ret <- eval(x, .rxOptGetEnv(x))
  return(..rxOpt(eval(parse(text=paste0("quote(", .ret, ")")))))
}

..rxOpt <- function(x, progress=FALSE) {
  if (is.atomic(x)) {
    return(as.character(x))
  } else if (is.name(x)) {
    return(as.character(x))
  } else if (is.call(x)) {
    .x2 <- x[-1];
    if (identical(x[[1]], quote(`{`))) {
      ## Probably can implement a progress bar here.
      if (progress) {
        rxProgress(length(.x2));
        on.exit({
          rxProgressAbort("Stopped optimizing duplicate expressions")
        });
        .ret <- unlist(lapply(.x2, function(x) {
          rxTick()
          ..rxOpt(x)
        }))
        rxProgressStop();
        return(.ret)
      } else {
        return(unlist(lapply(.x2, ..rxOpt)))
      }
    } else if (identical(x[[1]], quote(`(`))){
      .x2 <- x[[2]]
      .test <- TRUE
      while (.test){
        if (length(.x2) == 2){
          if (identical(.x2[[1]], quote(`(`))) {
            .x2 <- .x2[[2]]
          } else {
            .test <- FALSE
          }
        } else {
          .test <- FALSE
        }
      }
      return(paste0("(", ..rxOpt(.x2), ")"))
    } else if (identical(x[[1]], quote(`*`)) ||
                 identical(x[[1]], quote(`^`)) ||
                 identical(x[[1]], quote(`+`)) ||
                 identical(x[[1]], quote(`-`)) ||
                 identical(x[[1]], quote(`/`)) ||
                 identical(x[[1]], quote(`==`)) ||
                 identical(x[[1]], quote(`>=`)) ||
                 identical(x[[1]], quote(`<=`)) ||
                 identical(x[[1]], quote(`>`)) ||
                 identical(x[[1]], quote(`<`)) ||
                 identical(x[[1]], quote(`!=`)) ||
                 identical(x[[1]], quote(`&&`)) ||
                 identical(x[[1]], quote(`||`)) ||
                 identical(x[[1]], quote(`&`)) ||
                 identical(x[[1]], quote(`|`))) {
      if (length(x) == 3) {
        if (length(x[[2]]) == 2) {
          if (identical(x[[2]][[1]], quote(`-`)) &&
                is.atomic(x[[2]][[2]])) {
            if (is.atomic(x[[3]])){
              if (identical(x[[1]], quote(`/`))){
                return(as.character(-x[[2]][[2]] / x[[3]]))
              } else if (identical(x[[1]], quote(`+`))) {
                return(as.character(-x[[2]][[2]] + x[[3]]))
              } else if (identical(x[[1]], quote(`-`))) {
                return(as.character(-x[[2]][[2]] - x[[3]]))
              } else if (identical(x[[1]], quote(`*`))) {
                return(as.character(-x[[2]][[2]] * x[[3]]))
              }
            }
          }
          if (x[[2]][[2]] == 1 &&
                identical(x[[1]], quote(`*`))) {
            return(paste0("-", ..rxOpt(x[[3]])))
          }
        }
        if (is.atomic(x[[2]]) && is.atomic(x[[3]])){
          if (identical(x[[1]], quote(`/`))){
            return(as.character(x[[2]] / x[[3]]))
          } else if (identical(x[[1]], quote(`+`))) {
            return(as.character(x[[2]] + x[[3]]))
          } else if (identical(x[[1]], quote(`-`))) {
            return(as.character(x[[2]] - x[[3]]))
          } else if (identical(x[[1]], quote(`*`))) {
            return(as.character(x[[2]] * x[[3]]))
          }
        }
        if (is.atomic(x[[2]])){
          if (x[[2]] == 1) {
            if (identical(x[[1]], quote(`*`))) {
              return(..rxOpt(x[[3]]))
            }
          }
          if (x[[2]] == 0){
            if (identical(x[[1]], quote(`*`))){
              return("0")
            } else if (identical(x[[1]], quote(`+`))) {
              return(..rxOpt(x[[3]]))
            } else if (identical(x[[1]], quote(`-`))) {
              return(paste0("-", ..rxOpt(x[[3]])))
            } else if (identical(x[[1]], quote(`/`))) {
              return("0")
            }
          }
        }
        if (is.atomic(x[[3]])){
          if (x[[3]] == 1) {
            if (identical(x[[1]], quote(`*`))) {
              return(..rxOpt(x[[2]]))
            }
          }
          if (x[[3]] == 0){
            if (identical(x[[1]], quote(`*`))){
              return("0")
            } else if (identical(x[[1]], quote(`+`))) {
              return(as.character(x[[2]]))
            } else if (identical(x[[1]], quote(`-`))) {
              return(as.character(x[[2]]))
            } else if (identical(x[[2]], quote(`/`))) {
              stop("cannot divide by zero")
            }
          }
        }
        return(paste0(..rxOpt(x[[2]]), as.character(x[[1]]),
                      ..rxOpt(x[[3]])));
      } else {
        ## Unary Operators
        return(paste0(as.character(x[[1]]),
                     ..rxOpt(x[[2]])))
      }
    } else if (identical(x[[1]], quote(`~`)) ||
                 identical(x[[1]], quote(`=`)) ||
                 identical(x[[1]], quote(`<-`))) {
      .rxOptEnv$.new <- c();
      .ret <- paste0(..rxOpt(x[[2]]),
                     ifelse(identical(x[[1]], quote(`<-`)),
                            "=", as.character(x[[1]])),
                     .rxOptExpr(x[[3]]))
      .extra <- c()
      if (length(.rxOptEnv$.new) > 0) {
        for (.i in seq_along(.rxOptEnv$.rep)) {
          if (any(.rxOptEnv$.rep[[.i]] == .rxOptEnv$.new) &&
                !any(.rxOptEnv$.rep[[.i]] == .rxOptEnv$.added)) {
            .cur <- c(.rxOptEnv$.rep[[.i]]);
            if (.i != 1) {
              for (.j in seq(1, .i - 1)) {
                while (!any(.rxOptEnv$.rep[[.j]] == .rxOptEnv$.added) &&
                      regexpr(rex::rex(or(.cur)),
                              names(.rxOptEnv$.rep)[.j]) != -1) {
                  .extra <- c(.extra, paste0(.rxOptEnv$.rep[[.j]],
                                             "~", names(.rxOptEnv$.rep)[.j]));
                  .rxOptEnv$.added <- c(.rxOptEnv$.added,
                                        .rxOptEnv$.rep[.j]);
                }
              }
            }
            .extra <- c(.extra,
                        paste0(.rxOptEnv$.rep[[.i]], "~",
                               ..rxOpt(eval(parse(text=paste0("quote(", names(.rxOptEnv$.rep)[.i], ")"))))))
            .rxOptEnv$.added <- c(.rxOptEnv$.added,
                                  .rxOptEnv$.rep[.i]);
          }
        }
      }
      return(paste(c(.extra, .ret), collapse = "\n"));
    } else if (identical(x[[1]], quote(`[`))) {
      return(paste0(..rxOpt(x[[2]]), "[", ..rxOpt(x[[3]]), "]"));
    } else {
      .ret0 <- lapply(x, ..rxOpt);
      .ret <- paste0(.ret0[[1]], "(")
      if (.ret == "((") .ret <- "("
      .ret0 <- .ret0[-1];
      .ret <- paste0(.ret, paste(unlist(.ret0), collapse = ", "), ")");
      return(.ret)
    }
  }
}

##' Optimize RxODE for computer evaluation
##'
##' This optimizes RxODE code for computer evaluation by only
##' calculating redundant expressions once.
##'
##' @param x RxODE model that can be access by rxNorm
##'
##' @param msg This is the name of type of object that RxODE is
##'     optimizing that will in the message when optimizing.  For
##'     example "model" will produce the following message while
##'     optimizing the model:
##'
##'  finding duplicate expressions in model...
##'
##' @return Optimized RxODE model text.  The order and type lhs and
##'     state variables is maintained while the evaluation is sped up.
##'     While parameters names are maintained, their order may be
##'     modified.
##'
##' @author Matthew L. Fidler
##' @export
rxOptExpr <- function(x, msg="model") {
  .oldOpts <- options()
  options(digits=22)
  on.exit(options(.oldOpts))
  .rxOptEnv$.list <- list();
  .rxOptEnv$.rep <- list();
  .rxOptEnv$.added <- c();
  .rxOptEnv$.exclude <- "";
  message(sprintf("finding duplicate expressions in %s...", msg))
  .p <- eval(parse(text = paste0("quote({", x, "})")))
  .lines <- ..rxOpt(.p, progress = TRUE)
  .rxOptEnv$.list <- .rxOptEnv$.list[which(unlist(.rxOptEnv$.list) > 1L)];
  .exprs <- names(.rxOptEnv$.list)[order(nchar(names(.rxOptEnv$.list)))];
  .exprs <- .exprs[regexpr(rex::rex(start, regNum, end), .exprs,
                           perl = TRUE) == -1]
  .thetaEtaR <- rex::rex(start, or("THETA[", "ETA["), any_numbers, "]", end)
  .exprs <- .exprs[regexpr(.thetaEtaR, .exprs, perl = TRUE) == -1]
  if (length(.exprs) > 0) {
    ## Take out unary [-] that way
    ## expr1=-ka       #nolint
    ## expr2 = expr-ka #nolint
    ## will not become expr = exprka where exprka isn't defined.
    .exprs <- .exprs[regexpr("^[-]", .exprs) == -1];
    if (length(.exprs) == 0){
      return(x)
    }
    .rp <- rxOptRep_(.exprs)
    .rxOptEnv$.rep <- as.list(.rp[[1]]);
    .rxOptEnv$.exclude <- ""
    message(sprintf("optimizing duplicate expressions in %s...", msg))
    .opt <- ..rxOpt(.p, progress = TRUE);
    return(paste(.opt, collapse = "\n"))
  } else {
    return(paste(.lines, collapse ="\n"))
  }
}
