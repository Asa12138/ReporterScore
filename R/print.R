#' Print reporter_score
#'
#' @param x reporter_score
#' @param ... add
#'
#' @return No value
#' @exportS3Method
#' @method print reporter_score
#'
print.reporter_score <- function(x, ...) {
  reporter_score_res <- x
  pcutils::dabiao("KO abundance table", print = TRUE)
  cat("With ", nrow(reporter_score_res$kodf), " KOs and ", ncol(reporter_score_res$kodf), " samples.\n")
  pcutils::dabiao("group", print = TRUE)
  vs_group <- attributes(reporter_score_res$reporter_s)$vs_group

  if (attributes(reporter_score_res$reporter_s)$mode == "directed") {
    title <- paste0(vs_group, collapse = "/ ")
  } else {
    title <- paste0(vs_group, collapse = "/ ")
  }
  cat("vs group: ", title, "\n", sep = "")
  pcutils::dabiao("parameter", print = TRUE)
  cat("use mode: ", attributes(reporter_score_res$reporter_s)$mode, "; use method: ", attributes(reporter_score_res$reporter_s)$method, "; the feature: ", attributes(reporter_score_res$reporter_s)$feature,
    sep = ""
  )
}


#' Print rs_by_cm
#'
#' @param x rs_by_cm
#' @param ... add
#'
#' @return No value
#' @exportS3Method
#' @method print rs_by_cm
#'
print.rs_by_cm <- function(x, ...) {
  rsa_cm_res <- x
  ncluster <- sum(grepl("Cluster", names(rsa_cm_res)))
  pcutils::dabiao("RSA result with ", ncluster, " clusters", print = TRUE)

  pcutils::dabiao("KO abundance table", print = TRUE)
  cat("With ", nrow(rsa_cm_res$kodf), " KOs and ", ncol(rsa_cm_res$kodf), " samples.\n")
  pcutils::dabiao("group", print = TRUE)
  vs_group <- attributes(rsa_cm_res$Cluster1$reporter_s)$vs_group

  title <- paste0(vs_group, collapse = "/ ")
  cat("vs group: ", title, "\n", sep = "")

  cat("Patterns: \n")
  print(rsa_cm_res$cm_res$centers)

  pcutils::dabiao("parameter", print = TRUE)
  cat("use method: ", attributes(rsa_cm_res$Cluster1$reporter_s)$method, "; the feature: ", attributes(rsa_cm_res$Cluster1$reporter_s)$feature, sep = "")
}
