##' Replace comments with label()
##'
##' @param src source function
##' @return character vector with labels("") replaced in function
##' @author Matthew Fidler
##' @noRd
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

##' Convert RxODE/nlmixr model function to a string
##'
##' @param fun function name for parsing
##' @return Modified function
##' @author Matthew Fidler
##' @noRd
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

.lastIni <- NULL
##' Ini block for RxODE/nlmixr models
##'
##' @param x expression
##' @param ... Other expressions for `ini()` function
##' @return Ini block
##' @author Matthew Fidler
##' @export
ini <- function(x, ..., envir = parent.frame()) {
  if (is(substitute(x), "{")) {
    .ini <- eval(bquote(lotri(.(substitute(x)))), envir=envir)
    assignInMyNamespace(".lastIni", .ini)
    return(invisible(.ini))
  }
  UseMethod("ini")
}

##' @export
##' @rdname ini
ini.default <- function(x, ...) {

}

##' Model block for RxODE/nlmixr models
##'
##' @param x model expression
##' @param ... Other arguments
##' @return Model block with ini information included.  `ini` must be called before `model` block
##' @author Matthew Fidler
##' @export
model <- function(x, ...) {
  if (is(substitute(x), "{")) {
    .ini <- .lastIni
    if (is.null(.ini)) {
      stop("ini({}) block must be called before the model block",
           call.=FALSE)
    }
    assignInMyNamespace(".lastIni", NULL)
    .mod <- eval(bquote(.errProcessExpression(quote(.(substitute(x))), .ini)))
    return(.mod)
   }
  UseMethod("model")
}

##' @export
##' @rdname model
model.default <- function(x, ...) {

}


