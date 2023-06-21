#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' Some functions from other packages
#' @importFrom pcutils lib_ps
#' @importFrom pcutils get_cols
#' @importFrom pcutils mmscale
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr left_join
#' @importFrom dplyr arrange
#' @importFrom dplyr count
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr rename
#' @importFrom dplyr distinct
#' @name init_ReporterScore
