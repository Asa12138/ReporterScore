% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{reporter_score_res}
\alias{reporter_score_res}
\alias{reporter_score_res2}
\title{`reporter_score()` result from KO_abundance_test}
\format{
a list contain 7 elements.
\describe{
\item{kodf}{your input KO_abundance table}
\item{ko_stat}{ko statistics result contains p.value and z_score}
\item{reporter_s}{the reporter score in each pathway}
\item{modulelist}{default KOlist or customized modulelist dataframe}
\item{group}{The compare group (two category) in your data}
\item{metadata}{sample information dataframe contains group}
}
}
\description{
`reporter_score()` result from KO_abundance_test

`reporter_score()` result from KO_abundance_test
}
\seealso{
Other test_data: 
\code{\link{KO_abundance}},
\code{\link{genedf}}
}
\concept{test_data}
