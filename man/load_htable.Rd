% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_files.R
\name{load_htable}
\alias{load_htable}
\alias{load_KOlist}
\alias{load_CPDlist}
\alias{load_KO_desc}
\alias{load_KO_htable}
\alias{load_Pathway_htable}
\alias{load_Enzyme_htable}
\alias{load_Module_htable}
\alias{load_Compound_htable}
\alias{load_org_pathway}
\title{Load the specific table (from 'KEGG')}
\usage{
load_htable(type, verbose = TRUE)

load_KOlist(verbose = TRUE)

load_CPDlist(verbose = TRUE)

load_KO_desc(verbose = TRUE)

load_KO_htable(verbose = TRUE)

load_Pathway_htable(verbose = TRUE)

load_Enzyme_htable(verbose = TRUE)

load_Module_htable(verbose = TRUE)

load_Compound_htable(verbose = TRUE)

load_org_pathway(org = "hsa", verbose = TRUE)
}
\arguments{
\item{type}{"ko", "module", "pathway", "compound" ...}

\item{verbose}{logical}

\item{org}{kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html, default, "hsa"}
}
\value{
KO_htable

KOlist

CPDlist

KO description

KO_htable

Pathway_htable

Enzyme_htable

Module_htable

Compound_htable

KOlist
}
\description{
Load the specific table (from 'KEGG')

Load the KOlist (from 'KEGG')

Load the CPDlist (from 'KEGG')

Load the KO description (from 'KEGG')

Load the KO_htable (from 'KEGG')

Load the Pathway_htable (from 'KEGG')

Load the Enzyme_htable (from 'KEGG')

Load the Module_htable (from 'KEGG')

Load the Compound_htable (from 'KEGG')

Load the pathway information for an organism (from 'KEGG')
}
\examples{
Pathway_htable <- load_htable("pathway")
head(Pathway_htable)
}
