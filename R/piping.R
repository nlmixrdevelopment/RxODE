#'  This copies the RxODE UI object so it can be modifed
#'
#' @param ui Original UI object
#' @return Copied UI object
#' @author Matthew L. Fidler
#' @noRd
.copyUi <- function(ui) {
  .ret <- new.env(parent=emptyenv())
  lapply(ls(ui, envir=ui, all.names=TRUE), function(item){
    assign(item, get(item, envir=ui), envir=.ret)
  })
  class(.ret) <- class(ui)
  .ret
}

#' This expands a list of expressions
#'
#' @param lines These are the expressions as a list
#' @param brackets This is the indicator of the bracket lines ie `{}` that are expanded
#' @return Single list of expressions; `a=b` becomes `a<-b` in this expression
#' @author Matthew L. Fidler
#' @noRd
.quoteExpandBrackets <- function(lines, brackets) {
  if (length(brackets) == 0) return(lines)
  .expandedForm <- NULL
  .currentLine <- 1
  for (.b in brackets) {
    .bracketExpression <- lines[[.b]]
    .unlistedBrackets <- lapply(seq_along(.bracketExpression)[-1],
                                function(i) {
                                  .c <- .bracketExpression[[i]]
                                  if (identical(.c[[1]], quote(`=`))) {
                                    .c[[1]] <- quote(`<-`)
                                  }
                                  .c
                                })
    if (.currentLine == .b) {
      .expandedForm <- c(.expandedForm, .unlistedBrackets)
    } else {
      .expandedForm <- c(.expandedForm, lines[seq(.currentLine, .b - 1)],
                         .unlistedBrackets)
    }
    .currentLine <- .b + 1
  }
  if (.currentLine <= length(lines)) {
    .expandedForm <- c(.expandedForm, lines[seq(.currentLine, length(lines))])
  }
  .expandedForm
}


#'  Returns quoted call information
#'
#' @param callInfo Call information
#'
#' @return Quote call information.  for `name=expression`, change to
#'   `name<-expression` in quoted call list. For expressions that are
#'   within brackets ie `{}`, unlist the brackets as if they were
#'   called in one single sequence.
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.quoteCallInfoLines <- function(callInfo) {
  .bracket <- rep(FALSE, length=length(callInfo))
  .env <- environment()
  .ret <- lapply(seq_along(callInfo), function(i) {
    .name <- names(callInfo)[i]
    if (!is.null(.name)) {
      if (.name != "") {
        # Changed named items to
        return(as.call(list(quote(`<-`), .enQuote(.name),
                            eval(call("quote", callInfo[[i]])))))
      }
    }
    .quoted <- eval(call("quote", callInfo[[i]]))
    if (identical(.quoted[[1]], quote(`{`))) {
      .bracket[i] <- TRUE
      assign(".bracket", .bracket, envir=.env)
    }
    .quoted
  })
  .w <- which(.bracket)
  .quoteExpandBrackets(.ret, .w)
}
