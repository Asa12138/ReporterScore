#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom utils data
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom stats p.adjust reorder
#' @importFrom foreach %dopar%
#' @importFrom pcutils lib_ps get_cols
#' @import dplyr
#' @import ggplot2
## usethis namespace: end
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
