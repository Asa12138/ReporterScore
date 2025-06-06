% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kegg_net.R
\name{c_net_from_pathway_xml}
\alias{c_net_from_pathway_xml}
\title{Create a network from KEGG pathway XML files}
\usage{
c_net_from_pathway_xml(pathway_xml)
}
\arguments{
\item{pathway_xml}{A `tbl_graph` or `igraph` object, or a file path to a KEGG XML file.}
}
\value{
A `metanet` object representing the pathway network.
}
\description{
Create a network from KEGG pathway XML files
}
\seealso{
plot_pathway_net
}
