#' type_sum function for units
#' @name tibble
#' @param x see \link[pillar]{type_sum}
#' @param ... see \link[pillar]{type_sum}
#' @param width see \link[pillar]{type_sum}
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


## Now ggforce

##'@export
scale_x_unit <- ggforce::scale_x_unit

##'@export
scale_y_unit <- ggforce::scale_y_unit

##'@export
ScaleContinuousPositionUnit <- ggforce::ScaleContinuousPositionUnit

##'@export
scale_type <- ggplot2::scale_type

##'@export
ggplot <- ggplot2::ggplot

##'@export
aes <- ggplot2::aes

##'@export
geom_line <- ggplot2::geom_line

##'@export
facet_wrap <- ggplot2::facet_wrap

##'@export
scale_type.units <- loadNamespace("ggforce")$scale_type.units


#' @importFrom ggforce facet_wrap_paginate
#' @export
ggforce::facet_wrap_paginate


#' @importFrom ggforce facet_grid_paginate
#' @export
ggforce::facet_grid_paginate

#' @importFrom ggforce ScaleContinuousPositionUnit
#' @export
ggforce::ScaleContinuousPositionUnit

#' @importFrom ggplot2 aes
#' @export
ggplot2::aes

#' @importFrom ggplot2 facet_wrap
#' @export
ggplot2::facet_wrap


#' @importFrom ggplot2 geom_line
#' @export
ggplot2::geom_line


#' @importFrom ggplot2 ggplot
#' @export
ggplot2::ggplot

#' @importFrom ggplot2 scale_type
#' @export
ggplot2::scale_type

#' @importFrom ggforce scale_x_unit
#' @export
ggforce::scale_x_unit


#' @importFrom ggforce scale_y_unit
#' @export
ggforce::scale_y_unit


#' @importFrom units set_units
#' @export
units::set_units

#' @importFrom units drop_units
#' @export
units::drop_units

##'@importFrom ggplot2 ylab
##'@export
ggplot2::ylab

##'@importFrom ggplot2 xlab
##'@export
ggplot2::xlab

##'@importFrom lotri lotri
##'@export
lotri::lotri
