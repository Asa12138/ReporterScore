% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculation.R
\name{plot_KOs_box}
\alias{plot_KOs_box}
\title{Plot KOs boxplot}
\usage{
plot_KOs_box(
  kodf,
  group = NULL,
  metadata = NULL,
  map_id = "map00780",
  select_ko = NULL,
  box_param = NULL,
  KOlist_file = NULL,
  modulelist = NULL
)
}
\arguments{
\item{kodf}{KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples}

\item{group}{The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.}

\item{metadata}{metadata}

\item{map_id}{the pathway or module id}

\item{select_ko}{select which ko}

\item{box_param}{parameters pass to \code{\link[pcutils]{group_box}}}

\item{KOlist_file}{default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}}

\item{modulelist}{NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.}
}
\description{
Plot KOs boxplot
}
\examples{
data("reporter_score_res")
plot_KOs_box(reporter_score_res,"Group",metadata,
     select_ko=c("K00059","K00208","K00647","K00652","K00833","K01012"))

}
