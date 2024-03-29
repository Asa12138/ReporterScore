#' @title The KOlist used for enrichment.
#' @description an list contains two data.frame named pathway and module.
#'
#' @family data
#' @docType data
#' @name KOlist
#' @format four columns in each data.frame.
#' \describe{
#' \item{id}{"map0010" or "M00001"}
#' \item{K_num}{contians how many KOs in this pathway or module}
#' \item{KOs}{KOs name}
#' \item{Description}{the description of this pathway or module}
#' }
#'
NULL

#' @title The CPDlist used for enrichment.
#' @description an list contains two data.frame named pathway and module.
#'
#' @family data
#' @docType data
#' @name CPDlist
#' @format four columns in each data.frame.
#' \describe{
#' \item{id}{"map0010" or "M00001"}
#' \item{K_num}{contians how many Compounds in this pathway or module}
#' \item{KOs}{Compounds name}
#' \item{Description}{the description of this pathway or module}
#' }
#'
NULL

#' @title The GOlist used for enrichment.
#' @description an list contains three data.frame named BP, CC, MF.
#'
#' @family data
#' @docType data
#' @name GOlist
#' @format four columns in each data.frame.
#' \describe{
#' \item{id}{"map0010" or "M00001"}
#' \item{K_num}{contians how many Genes in this GO term}
#' \item{KOs}{Genes name}
#' \item{Description}{the description of this GO term}
#' }
#'
NULL

#' @title KO htable from 'KEGG'
#'
#' @family data
#' @docType data
#' @name KO_htable
#'
NULL

#' @title Module htable from 'KEGG'
#'
#' @family data
#' @docType data
#' @name Module_htable
#'
NULL

#' @title Pathway htable from 'KEGG'
#'
#' @family data
#' @docType data
#' @name Pathway_htable
#'
NULL

#' @title Compound htable from 'KEGG'
#'
#' @family data
#' @docType data
#' @name Compound_htable
#'
NULL

#' @title The KOs abundance table and group table.
#'
#' @family test_data
#' @docType data
#' @name KO_abundance
#'
NULL

#' @title The KOs abundance table and group table.
#'
#' @docType data
#' @name metadata
#' @rdname KO_abundance
NULL

#' @title `reporter_score()` result from KO_abundance_test
#'
#' @family test_data
#' @docType data
#' @name reporter_score_res
#' @format a list contain 7 elements.
#' \describe{
#' \item{kodf}{your input KO_abundance table}
#' \item{ko_stat}{ko statistics result contains p.value and z_score}
#' \item{reporter_s}{the reporter score in each pathway}
#' \item{modulelist}{default KOlist or customized modulelist dataframe}
#' \item{group}{The compare group (two category) in your data}
#' \item{metadata}{sample information dataframe contains group}
#' }
NULL

#' @title `reporter_score()` result from KO_abundance_test
#'
#' @docType data
#' @name reporter_score_res2
#' @rdname reporter_score_res
NULL

#' @title pathway information for "hsa"
#'
#' @family data
#' @docType data
#' @name hsa_kegg_pathway
NULL

#' @title pathway information for "mmu"
#'
#' @family data
#' @docType data
#' @name mmu_kegg_pathway
NULL

#' @title human gene table
#'
#' @family test_data
#' @docType data
#' @name genedf
NULL
