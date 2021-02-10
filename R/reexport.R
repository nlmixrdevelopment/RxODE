#' type_sum function for units
#' @name tibble
#' @param x see [type_sum][pillar::type_sum]
#' @param ... see [type_sum][pillar::type_sum]
#' @param width see [type_sum][pillar::type_sum]
#' @return type_sum units
type_sum.units <- loadNamespace("units")$type_sum.units

#' @name tibble
format_type_sum.type_sum_units <- loadNamespace("units")$format_type_sum.type_sum_units

#' @name tibble
pillar_shaft.units <- loadNamespace("units")$pillar_shaft.units

#' @name tibble
type_sum.mixed_units <- loadNamespace("units")$type_sum.mixed_units

#' @name tibble
pillar_shaft.mixed_units <- loadNamespace("units")$pillar_shaft.mixed_units
function(x, ...) {
  if (!requireNamespace("pillar", quietly = TRUE)) {
    stop("package 'xpillar' not available", .call = FALSE)
  }
  out <- format(x, ...)
  pillar::new_pillar_shaft_simple(out, align = "right", min_width = 6)
}

#' @export
scale_type <- ggplot2::scale_type

#' @export
ggplot <- ggplot2::ggplot

#' @export
aes <- ggplot2::aes

#' @export
geom_line <- ggplot2::geom_line

#' @export
facet_wrap <- ggplot2::facet_wrap


#' @importFrom ggplot2 %+replace%
`%+replace%`

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

#' @importFrom units set_units
#' @export
units::set_units

#' @importFrom units drop_units
#' @export
units::drop_units

#' @importFrom ggplot2 ylab
#' @export
ggplot2::ylab

#' @importFrom ggplot2 xlab
#' @export
ggplot2::xlab

#' @importFrom ggplot2 waiver
#' @export
ggplot2::waiver

#' Empty Guide
#'
#' This empty guide draws nothing; It is included in RxODE for
#' compatibility with ggplot 3.2
#'
#' @inheritParams ggplot2::guide_none
#' @return nothing, simply included to be compatible with ggplot 3.2
#' @export
guide_none <- function(title = waiver(), position = waiver()) {
  stop("needs \"ggplot2\" 3.3.0", call. = FALSE)
}

#' @importFrom lotri lotri
#' @export
lotri::lotri

#' @importFrom ggplot2  label_value
#' @export
ggplot2::label_value

#' @importFrom ggplot2 label_both
#' @export
ggplot2::label_both

#' @importFrom ggplot2 label_context
#' @export
ggplot2::label_context

#' @importFrom ggplot2 label_wrap_gen
#' @export
ggplot2::label_wrap_gen

#' @importFrom ggplot2 label_context
#' @export
ggplot2::label_context

#' @importFrom ggplot2 scale_x_discrete
#' @export
ggplot2::scale_x_discrete

#' @importFrom ggplot2 scale_y_discrete
#' @export
ggplot2::scale_y_discrete

#' @importFrom ggplot2 scale_x_continuous
#' @export
ggplot2::scale_x_continuous

#' @importFrom ggplot2 scale_y_continuous
#' @export
ggplot2::scale_y_continuous

#' @importFrom ggplot2 scale_x_date
#' @export
ggplot2::scale_x_date

#' @importFrom ggplot2 scale_y_date
#' @export
ggplot2::scale_y_date

#' @importFrom ggplot2 expand_limits
#' @export
ggplot2::expand_limits


.SD <- NULL
`:=` <- function(...) { ## nocov start
  stop("this is only used in 'data.table'", .call = FALSE)
} ## nocov end
