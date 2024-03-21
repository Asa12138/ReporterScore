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
#' @examples
#' seq_len(5) %>% sum()
#'
NULL


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
#' @export get_KOs
deprecated("get_KOs", get_features)
#' @export plot_report_bar
deprecated("plot_report_bar", plot_report)
#' @export download_org_pathway
deprecated("download_org_pathway", update_org_pathway)
#' @export get_org_pathway
deprecated("get_org_pathway", update_org_pathway)
#' @export plot_KOs_network
deprecated("plot_KOs_network", plot_features_network)
