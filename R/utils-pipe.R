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

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


deprecated <- function(old, new) {
    assign(old, new, envir = asNamespace(packageName()))
}

#' @export GRSA
deprecated("GRSA", reporter_score)
#' @export RSA
deprecated("RSA", reporter_score)
#' @export GRSA_by_cm
deprecated("GRSA_by_cm", RSA_by_cm)
#' @export plot_KOs_distribution
deprecated("plot_KOs_distribution", plot_features_distribution)
#' @export plot_KOs_in_pathway
deprecated("plot_KOs_in_pathway", plot_features_in_pathway)
#' @export plot_KOs_heatmap
deprecated("plot_KOs_heatmap", plot_features_heatmap)
#' @export plot_KOs_box
deprecated("plot_KOs_box", plot_features_box)
#' @export plot_KOs_network
deprecated("plot_KOs_network", plot_features_network)
#' @export get_KOs
deprecated("get_KOs", get_features)
#' @export plot_report_bar
deprecated("plot_report_bar", plot_report)
#' @export download_org_pathway
deprecated("download_org_pathway", get_org_pathway)

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
