## Generated code from build/refresh.R

##' @name spread_
##' @export spread_.solveRxDll
##'
##' @method spread_ solveRxDll
##'
##' @title spread_ for \code{solveRxDll} object
##' @description tidyr compatability layer for solveRxDll
##' @param data Solved ODE, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
spread_.solveRxDll <- function(data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("spread_","tidyr"), call, envir = parent.frame(1)));
}

##' @name unite_
##' @export unite_.solveRxDll
##'
##' @method unite_ solveRxDll
##'
##' @title unite_ for \code{solveRxDll} object
##' @description tidyr compatability layer for solveRxDll
##' @param data Solved ODE, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
unite_.solveRxDll <- function(data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("unite_","tidyr"), call, envir = parent.frame(1)));
}

##' @name separate_
##' @export separate_.solveRxDll
##'
##' @method separate_ solveRxDll
##'
##' @title separate_ for \code{solveRxDll} object
##' @description tidyr compatability layer for solveRxDll
##' @param data Solved ODE, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
separate_.solveRxDll <- function(data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("separate_","tidyr"), call, envir = parent.frame(1)));
}

##' @name gather_
##' @export gather_.solveRxDll
##'
##' @method gather_ solveRxDll
##'
##' @title gather_ for \code{solveRxDll} object
##' @description tidyr compatability layer for solveRxDll
##' @param data Solved ODE, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
gather_.solveRxDll <- function(data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$data <- dplyr::as.tbl(data)
  return(do.call(getFromNamespace("gather_","tidyr"), call, envir = parent.frame(1)));
}

##' @name sample_frac
##' @export sample_frac.solveRxDll
##'
##' @method sample_frac solveRxDll
##'
##' @title sample_frac for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param tbl Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
sample_frac.solveRxDll <- function(tbl, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$tbl <- dplyr::as.tbl(tbl)
  return(do.call(getFromNamespace("sample_frac","dplyr"), call, envir = parent.frame(1)));
}

##' @name sample_n
##' @export sample_n.solveRxDll
##'
##' @method sample_n solveRxDll
##'
##' @title sample_n for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param tbl Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
sample_n.solveRxDll <- function(tbl, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$tbl <- dplyr::as.tbl(tbl)
  return(do.call(getFromNamespace("sample_n","dplyr"), call, envir = parent.frame(1)));
}

##' @name group_by_
##' @export group_by_.solveRxDll
##'
##' @method group_by_ solveRxDll
##'
##' @title group_by_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
group_by_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("group_by_","dplyr"), call, envir = parent.frame(1)));
}

##' @name rename_
##' @export rename_.solveRxDll
##'
##' @method rename_ solveRxDll
##'
##' @title rename_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
rename_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("rename_","dplyr"), call, envir = parent.frame(1)));
}

##' @name arrange_
##' @export arrange_.solveRxDll
##'
##' @method arrange_ solveRxDll
##'
##' @title arrange_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
arrange_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("arrange_","dplyr"), call, envir = parent.frame(1)));
}

##' @name summarise_
##' @export summarise_.solveRxDll
##'
##' @method summarise_ solveRxDll
##'
##' @title summarise_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
summarise_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("summarise_","dplyr"), call, envir = parent.frame(1)));
}

##' @name transmute_
##' @export transmute_.solveRxDll
##'
##' @method transmute_ solveRxDll
##'
##' @title transmute_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
transmute_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("transmute_","dplyr"), call, envir = parent.frame(1)));
}

##' @name mutate_
##' @export mutate_.solveRxDll
##'
##' @method mutate_ solveRxDll
##'
##' @title mutate_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
mutate_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("mutate_","dplyr"), call, envir = parent.frame(1)));
}

##' @name distinct_
##' @export distinct_.solveRxDll
##'
##' @method distinct_ solveRxDll
##'
##' @title distinct_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
distinct_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("distinct_","dplyr"), call, envir = parent.frame(1)));
}

##' @name rename_
##' @export rename_.solveRxDll
##'
##' @method rename_ solveRxDll
##'
##' @title rename_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
rename_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("rename_","dplyr"), call, envir = parent.frame(1)));
}

##' @name select_
##' @export select_.solveRxDll
##'
##' @method select_ solveRxDll
##'
##' @title select_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
select_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("select_","dplyr"), call, envir = parent.frame(1)));
}

##' @name arrange_
##' @export arrange_.solveRxDll
##'
##' @method arrange_ solveRxDll
##'
##' @title arrange_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
arrange_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("arrange_","dplyr"), call, envir = parent.frame(1)));
}

##' @name slice_
##' @export slice_.solveRxDll
##'
##' @method slice_ solveRxDll
##'
##' @title slice_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
slice_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
  return(do.call(getFromNamespace("slice_","dplyr"), call, envir = parent.frame(1)));
}

##' @name filter_
##' @export filter_.solveRxDll
##'
##' @method filter_ solveRxDll
##'
##' @title filter_ for \code{solveRxDll} object
##' @description dplyr compatability layer for solveRxDll
##' @param .data Solved equation, an \code{solveRxDll} object.
##' @param ... Additional arguments
##'
filter_.solveRxDll <- function(.data, ...){
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$.data <- asTbl(.data)
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

##' @importFrom stats aggregate
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

##' @importFrom utils edit
##' @export
edit.solveRxDll <- function (name, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$name <- as.data.frame(name)
  return(do.call("edit.data.frame", call, envir = parent.frame(1)));
}

##' @importFrom methods Math
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

##' @importFrom utils stack
##' @export
stack.solveRxDll <- function (x, ...) {
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$x <- as.data.frame(x)
  return(do.call("stack.data.frame", call, envir = parent.frame(1)));
}

##' @importFrom utils unstack
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

