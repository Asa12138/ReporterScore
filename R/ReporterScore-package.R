#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom utils data head tail download.file untar packageName
#' @importFrom stats p.adjust reorder cor cor.test rnorm sd var na.omit setNames median quantile
#' @importFrom pcutils lib_ps get_cols strsplit2
#' @import dplyr
#' @import ggplot2
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
utils::globalVariables(c(
    "CARDinfo", "CPDlist", "Compound_htable",
    "GOlist", "KO_htable", "KOlist", "Module_htable",
    "Pathway_htable", "GOinfo"
))
