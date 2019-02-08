#' type_sum function for units
#' @name tibble
#' @param x see \link[pillar]{type_sum}
#' @param ... see \link[pillar]{type_sum}
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(pillar::type_sum, units)
#'   S3method(pillar::type_sum, mixed_units)
#'   S3method(pillar::pillar_shaft, units)
#'   S3method(pillar::pillar_shaft, mixed_units)
#'   S3method(pillar::format_type_sum, type_sum_units)
#' } else {
#'   export(type_sum.units)
#'   export(type_sum.mixed_units)
#'   export(pillar_shaft.units)
#'   export(pillar_shaft.mixed_units)
#'   export(format_type_sum.type_sum_units)
#' }
type_sum.units <- loadNamespace("units")$type_sum.units
#'@name tibble
format_type_sum.type_sum_units  <- loadNamespace("units")$format_type_sum.type_sum_units
#'@name tibble
pillar_shaft.units <- loadNamespace("units")$pillar_shaft.units

#'@name tibble
type_sum.mixed_units <- loadNamespace("units")$type_sum.mixed_units

#' @name tibble
pillar_shaft.mixed_units <- loadNamespace("units")$pillar_shaft.mixed_units
function(x, ...) {
  if (! requireNamespace("pillar", quietly = TRUE))
    stop("package pillar not available: install first?")
  out <- format(x, ...)
  pillar::new_pillar_shaft_simple(out, align = "right", min_width = 6)
}
