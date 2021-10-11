#' Replace comments with label()
#'
#' @param src source function
#' @return character vector with labels("") replaced in function
#' @author Matthew Fidler
#' @noRd
.rxReplaceCommentWithLabel <- function(src) {
  .env <- new.env(parent=emptyenv())
  .env$inIni <- FALSE
  .env$convertLabel <- FALSE
  .regIni <- rex::rex(boundary, or("ini", "lotri"), "(", any_spaces, "{")
  .regOther <- rex::rex(boundary, or(group("}", any_spaces, ")")), group(or("model"), "(", any_spaces, "{"))
  .regCommentOnBlankLine <- "^ *#+ *(.*) *$"
  .regLabel <- "^( *[^\n\"]+) *#+ *(.*) *$"
  .ret <- vapply(src,
                 function(line) {
                   if (regexpr(.regIni, line, perl=TRUE) != -1) {
                     .env$inIni <- TRUE
                   } else if (regexpr(.regOther, line, perl=TRUE) != -1) {
                     .env$inIni <- FALSE
                   } else if (.env$inIni) {
                     if (regexpr(.regCommentOnBlankLine, line) != -1) {
                     } else if (regexpr(.regLabel, line) != -1) {
                       .env$convertLabel <- TRUE
                       .label <- deparse1(sub(.regLabel, "\\2", line))
                       return(sub(.regLabel, paste0("\\1; label(", .label, ")"), line))
                     }
                   }
                   line
                 }, character(1), USE.NAMES = FALSE)
  if (.env$convertLabel) {
    cli::cli_alert_info("parameter labels from comments will be replaced by 'label()'")
  }
  .ret <- deparse(eval(parse(text=paste(.ret, collapse="\n"), keep.source=FALSE)))
  .ret
}

#' Convert RxODE/nlmixr model function to a string
#'
#' @param fun function name for parsing
#' @return Modified function
#' @author Matthew Fidler
#' @noRd
.rxFunction2string <- function(fun) {
  .srcRef <- attr(fun, "srcref")
  if (is.null(.srcRef)) {
    cli::cli_alert_info("parameter labels from comments are typically ignored in non-interactive mode")
    cli::cli_alert_info("Need to run with the source intact to parse comments")
    .ret <- deparse(fun)
  } else {
    .ret <- .rxReplaceCommentWithLabel(as.character(.srcRef, useSource = TRUE))
  }
  .ret
}

.rxFunction2ui <- function(fun) {
  .fun <- eval(parse(text=paste(.rxFunction2string(fun), collapse="\n")))
  .fun()
}

.lastIni <- NULL
.lastIniQ <- NULL
#' Ini block for RxODE/nlmixr models
#'
#' @param x expression
#' @param ... Other expressions for `ini()` function
#' @return Ini block
#' @author Matthew Fidler
#' @export
ini <- function(x, ..., envir = parent.frame()) {
  if (is(substitute(x), "{")) {
    .ini <- eval(bquote(lotri(.(substitute(x)))), envir=envir)
    assignInMyNamespace(".lastIni", .ini)
    assignInMyNamespace(".lastIniQ", bquote(.(substitute(x))))
    return(invisible(.ini))
  }
  UseMethod("ini")
}

#' @export
#' @rdname ini
ini.default <- function(x, ...) {

}

#' Model block for RxODE/nlmixr models
#'
#' @param x model expression
#' @param ... Other arguments
#' @inheritParams eval
#' @return Model block with ini information included.  `ini` must be called before `model` block
#' @author Matthew Fidler
#' @export
model <- function(x, ..., envir=parent.frame()) {
  if (is(substitute(x), "{")) {
    .ini <- .lastIni
    .iniQ <- .lastIniQ
    if (is.null(.ini)) {
      stop("ini({}) block must be called before the model block",
           call.=FALSE)
    }
    assignInMyNamespace(".lastIni", NULL)
    assignInMyNamespace(".lastIniQ", NULL)
    .mod <- .rxMuRef(eval(bquote(.errProcessExpression(quote(.(substitute(x))), .ini))))
    .mod$.ini <- .iniQ
    .meta <- new.env(parent=emptyenv())
    if (!identical(envir, globalenv())) {
      for (.i in ls(envir, all=TRUE)) {
        assign(.i, get(.i, envir), .meta)
      }
    }
    .mod$meta <- .meta
    class(.mod) <- "rxUi"
    return(.mod)
   }
  UseMethod("model")
}

#' @export
#' @rdname model
model.default <- function(x, ...) {

}

#' S3 for getting information from UI model
#'
#' @param x list of (UIenvironment, exact).  UI environment is the parsed function for RxODE.  `exact` is a boolean that says
#' @param ... Other arguments
#' @return value that was requested from the UI object
#' @author Matthew Fidler
#' @export
rxUiGet <- function(x, ...) {
  if (!inherits(x, "rxUiGet")) {
    stop("object is wrong type for `rxUiGet`")
  }
  UseMethod("rxUiGet")
}

#' @rdname rxUiGet
#' @export
rxUiGet.funPrint <- function(x, ...) {
  .x <- x[[1]]
  .ls <- ls(.x$meta, all=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(.x$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- .x$iniFun
  .ret[[.len + 3]] <- .x$modelFun
  .ret
}

#' @export
#' @rdname rxUiGet
rxUiGet.fun <- function(x, ...) {
  .ret <- rxUiGet.funPrint(x, ...)
  .ret2 <- function(){}
  body(.ret2) <- as.call(.ret)
  .ret2
}

#' @export
#' @rdname rxUiGet
rxUiGet.ini <- function(x, ...) {
  get("iniDf", x[[1]])
}

#'@export
#' @rdname rxUiGet
rxUiGet.iniFun <- function(x, ...) {
  .x <- x[[1]]
  .arg <- class(x)[1]
  bquote(ini(.(.x$.ini)))
}

#' @export
#' @rdname rxUiGet
rxUiGet.modelFun <- function(x, ...) {
  .x <- x[[1]]
  bquote(model(.(as.call(c(quote(`{`),.x$lstExpr)))))
}

#' @export
#' @rdname rxUiGet
rxUiGet.default <- function(x, ...) {
  .arg <- class(x)[1]
  if (!exists(.arg, envir=x[[1]])) return(NULL)
  get(.arg, x[[1]])
}

#' @export
`$.rxUi` <- function(obj, arg, exact = TRUE) {
  .lst <- list(obj, exact)
  class(.lst) <- c(arg, "rxUiGet")
  rxUiGet(.lst)
}


#' @export
.DollarNames.rxUi <- function(x, pattern) {
  .cmp <- vapply(as.character(methods("rxUiGet")), function(x){substr(x,9, nchar(x))}, character(1), USE.NAMES = FALSE)
  .cmp <- c(.cmp[.cmp != "default"], ls(x, all=TRUE))
  grep(pattern, .cmp, value = TRUE)
}

