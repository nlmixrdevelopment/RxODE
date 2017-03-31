## Generated code from build/refresh.R


##' @export
spread_.solveRxDll <- function (data, key_col, value_col, fill = NA, convert = FALSE,
    drop = TRUE, sep = NULL) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("spread_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
unite_.solveRxDll <- function (data, col, from, sep = "_", remove = TRUE) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("unite_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
separate_.solveRxDll <- function (data, col, into, sep = "[^[:alnum:]]+", remove = TRUE,
    convert = FALSE, extra = "warn", fill = "warn", ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("separate_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
gather_.solveRxDll <- function (data, key_col, value_col, gather_cols, na.rm = FALSE,
    convert = FALSE, factor_key = FALSE) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("gather_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
sample_frac.solveRxDll <- function (tbl, size = 1, replace = FALSE, weight = NULL, .env = parent.frame()) {
  call <- as.list(match.call())[-1];
  call$tbl <- dplyr::as.tbl(tbl);
  return(do.call(getFromNamespace("sample_frac","dplyr"), call, envir = parent.frame(1)));
}


##' @export
sample_n.solveRxDll <- function (tbl, size, replace = FALSE, weight = NULL, .env = parent.frame()) {
  call <- as.list(match.call())[-1];
  call$tbl <- dplyr::as.tbl(tbl);
  return(do.call(getFromNamespace("sample_n","dplyr"), call, envir = parent.frame(1)));
}


##' @export
group_by_.solveRxDll <- function (.data, ..., .dots, add = FALSE) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("group_by_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
rename_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("rename_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
arrange_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("arrange_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
summarise_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("summarise_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
transmute_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("transmute_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
mutate_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("mutate_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
distinct_.solveRxDll <- function (.data, ..., .dots, .keep_all = FALSE) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("distinct_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
rename_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("rename_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
select_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("select_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
arrange_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("arrange_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
slice_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("slice_","dplyr"), call, envir = parent.frame(1)));
}


##' @export
filter_.solveRxDll <- function (.data, ..., .dots) {
  call <- as.list(match.call())[-1];
  call$.data <- asTbl(.data);
  return(do.call(getFromNamespace("filter_","dplyr"), call, envir = parent.frame(1)));
}

##' @export
row.names.solveRxDll <- function (x) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("row.names.data.frame", call, envir = parent.frame(1)));
}

##' @export
by.solveRxDll <- function (data, INDICES, FUN, ..., simplify = TRUE) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- as.data.frame(data)
  return(do.call("by.data.frame", call, envir = parent.frame(1)));
}


##' @export
aggregate.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("aggregate.data.frame", call, envir = parent.frame(1)));
}

##' @export
anyDuplicated.solveRxDll <- function (x, incomparables = FALSE, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("anyDuplicated.data.frame", call, envir = parent.frame(1)));
}

##' @export
droplevels.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("droplevels.data.frame", call, envir = parent.frame(1)));
}

##' @export
duplicated.solveRxDll <- function (x, incomparables = FALSE, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("duplicated.data.frame", call, envir = parent.frame(1)));
}


##' @export
edit.solveRxDll <- function (name, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$name <- as.data.frame(name)
  return(do.call("edit.data.frame", call, envir = parent.frame(1)));
}


##' @export
Math.solveRxDll <- function(x, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("Math.data.frame", call, envir = parent.frame(1)));
}

##' @export
rowsum.solveRxDll <- function (x, group, reorder = TRUE, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("rowsum.data.frame", call, envir = parent.frame(1)));
}

##' @export
split.solveRxDll <- function (x, f, drop = FALSE, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("split.data.frame", call, envir = parent.frame(1)));
}

##' @export
subset.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("subset.data.frame", call, envir = parent.frame(1)));
}


##' @export
stack.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("stack.data.frame", call, envir = parent.frame(1)));
}


##' @export
unstack.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("unstack.data.frame", call, envir = parent.frame(1)));
}

##' @export
unique.solveRxDll <- function (x, incomparables = FALSE, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("unique.data.frame", call, envir = parent.frame(1)));
}

##' @export
within.solveRxDll <- function (data, expr, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- as.data.frame(data)
  return(do.call("within.data.frame", call, envir = parent.frame(1)));
}

##' @export
with.solveRxDll <- function (data, expr, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- as.data.frame(data)
  return(do.call("with.data.frame", call, envir = parent.frame(1)));
}

##' @export
dim.solveRxDll <- function (x) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.matrix(x)
  return(do.call("dim", call, envir = parent.frame(1)));
}

##' @export
dimnames.solveRxDll <- function (x) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.matrix(x)
  return(do.call("dimnames", call, envir = parent.frame(1)));
}

##' @export
t.solveRxDll <- function (x) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.matrix(x)
  return(do.call("t", call, envir = parent.frame(1)));
}

