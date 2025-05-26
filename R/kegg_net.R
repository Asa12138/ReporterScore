#' Download KEGG pathway XML files and create networks
#'
#' @param download_dir Directory to save the downloaded XML files.
#'
#' @returns No value
#' @export
#'
update_pathway_xml_ls <- function(download_dir) {
  # 检查依赖包
  if (!requireNamespace("ggkegg", quietly = TRUE)) {
    stop("Package 'ggkegg' is required but not installed. Run: install.packages('ggkegg')")
  }

  # 创建目录
  if (!dir.exists(download_dir)) {
    dir.create(download_dir, recursive = TRUE)
  }
  pack_dir <- tools::R_user_dir("ReporterScore")
  if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

  # 加载通路ID表
  Pathway_htable <- load_Pathway_htable() # 假设该函数已定义
  pathway_ids <- gsub("map", "ko", Pathway_htable$Pathway_id)

  # 初始化结果列表和失败记录
  pathway_xml_ls <- list()
  failed_paths <- character()

  # 遍历所有通路ID
  for (path in pathway_ids) {
    message("Processing pathway: ", path)

    # 直接尝试，报错则跳过
    result <- tryCatch(
      {
        pathway_xml_ls[[path]] <- ggkegg::pathway(path, directory = download_dir)
        message("Success: ", path)
      },
      error = function(e) {
        failed_paths <<- c(failed_paths, path)
        message("Failed to download ", path, ": ", e$message)
        return(NULL)
      }
    )
  }

  # 保存结果（含失败记录）
  attributes(pathway_xml_ls)$download_time <- Sys.time()
  attributes(pathway_xml_ls)$failed_paths <- failed_paths

  save_path <- paste0(pack_dir, "/pathway_xml_ls.rda")
  save(pathway_xml_ls, file = save_path)
  message("\nResults saved to: ", save_path)

  # 打印失败汇总
  if (length(failed_paths) > 0) {
    warning(
      "Failed to download ", length(failed_paths), " pathways:\n",
      paste(failed_paths, collapse = ", ")
    )
  }
}

#' Create a network from KEGG pathway XML files
#'
#' @param pathway_xml A `tbl_graph` or `igraph` object, or a file path to a KEGG XML file.
#'
#' @returns A `metanet` object representing the pathway network.
#' @export
#' @seealso plot_pathway_net
c_net_from_pathway_xml <- function(pathway_xml) {
  lib_ps("MetaNet", "igraph", "ggkegg", library = FALSE)
  name <- NULL
  if (is.null(pathway_xml)) {
    return(NULL)
  }
  if (identical(class(pathway_xml), c("tbl_graph", "igraph"))) {
    g <- pathway_xml
  } else if (is.character(pathway_xml) && file.exists(pathway_xml)) {
    pid <- basename(pathway_xml)
    pid <- gsub("\\.xml$", "", pid) # Remove .xml extension if present
    g <- ggkegg::pathway(pid, directory = dirname(pathway_xml))
  } else {
    stop("pathway_xml must be a tbl_graph, igraph object or a file path to a KEGG XML file.")
  }

  node1 <- as.data.frame(g)
  pathway_id <- unique(node1$pathway_id)[1]
  edge1 <- igraph::as_data_frame(g, what = "edges")


  path_gene_compound <- dplyr::filter(node1, type %in% c("ortholog", "gene", "compound")) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  if (nrow(edge1) > 0) {
    path_net <- dplyr::filter(
      edge1, from %in% path_gene_compound$name,
      to %in% path_gene_compound$name
    )
  } else {
    path_net <- edge1
    path_net$subtype_name <- character()
  }
  MetaNet::c_net_from_edgelist(path_net, vertex_df = path_gene_compound, direct = TRUE) -> path_net_c
  # igraph::V(path_net_c)$label <- gsub(" .*", "", igraph::V(path_net_c)$name)
  igraph::V(path_net_c)$label <- igraph::V(path_net_c)$graphics_name
  MetaNet::c_net_set(path_net_c,
    vertex_group = "type", vertex_class = "type",
    edge_type = "subtype_name"
  ) -> path_net_c

  attributes(path_net_c)$pathway_id <- pathway_id
  path_net_c
}

#' Plot a KEGG pathway network
#'
#' @param path_net_c A `metanet` object representing the pathway network.
#' @param ... Additional arguments passed to `MetaNet::c_net_plot`.
#' @param simplify Logical, whether to simplify the network by removing multiple edges and loops. Default is `FALSE`.
#' @param plot_depth Logical, whether to plot the network as a tree layout. Default is `FALSE`.
#'
#' @returns A plot of the pathway network.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("MetaNet") && requireNamespace("ggkegg")) {
#'   tmp_dir <- tempdir()
#'   pcutils::download2("https://rest.kegg.jp/get/ko01521/kgml", file.path(tmp_dir, "ko01521.xml"))
#'   path_net_c <- c_net_from_pathway_xml(file.path(tmp_dir, "ko01521.xml"))
#'   plot_pathway_net(path_net_c)
#'   pathway_net_index(path_net_c)
#' }
#' }
plot_pathway_net <- function(path_net_c, simplify = FALSE, plot_depth = FALSE, ...) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  name <- NULL

  if (simplify) igraph::simplify(path_net_c, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "first") -> path_net_c

  default_arg <- list(
    labels_num = "all", main = attributes(path_net_c)$pathway_id,
    edge.arrow.size = 0.3, edge.arrow.width = 0.8
  )

  if (plot_depth) {
    default_arg <- append(default_arg, list(
      coors = igraph::as_tree(), rescale = TRUE
    ))
    do.call(MetaNet::c_net_plot, append(
      list(go = path_net_c),
      pcutils::update_param(default_arg, list(...))
    ))
  }
  MetaNet::get_v(path_net_c) -> tmp_v
  if (any(is.na(tmp_v$x)) || any(is.na(tmp_v$y))) {
    warning("Some node coordinates are NA. This may affect the plot layout.")
  }
  dplyr::filter(tmp_v, !is.na(x) & !is.na(y)) -> tmp_v
  MetaNet::c_net_filter(path_net_c, name %in% tmp_v$name) -> path_net_c
  do.call(MetaNet::c_net_plot, append(
    list(go = path_net_c),
    pcutils::update_param(default_arg, list(...))
  ))
}

#' Create a pathway network index data.frame
#'
#' @param path_net_c A `tbl_graph` or `igraph` object representing a KEGG pathway network.
#'
#' @returns A data frame containing vertex indices and their attributes, including in-degree and out-degree.
#' @export
#' @seealso plot_pathway_net
pathway_net_index <- function(path_net_c) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  message("Removing multiple edges and loops from the network.")
  igraph::simplify(path_net_c, remove.multiple = TRUE, remove.loops = TRUE) -> path_net_c2

  igraph::vertex.attributes(path_net_c2) <- igraph::vertex.attributes(path_net_c2)[c("name", "type")]
  MetaNet::net_par(path_net_c2, mode = "v")$v_index -> path_net_index
  path_net_index$In_degree <- igraph::degree(path_net_c2, mode = "in")
  path_net_index$Out_degree <- igraph::degree(path_net_c2, mode = "out")
  attributes(path_net_index)$pathway_id <- attributes(path_net_c)$pathway_id

  # 计算深度
  igraph::layout_as_tree(path_net_c) -> layout_tree

  data.frame(
    row.names = NULL,
    Pathway_id = attributes(path_net_c)$pathway_id,
    path_net_index,
    Depth = layout_tree[, 2]
  )
}


#' Create an index for all KEGG pathway networks
#'
#' @param pathway_xml_ls A list of KEGG pathway XML files, where each element is a `tbl_graph` or `igraph` object.
#'
#' @returns A data frame containing indices for all pathways, including pathway IDs and their attributes.
#' @export
#'
get_all_pathway_net_index <- function(pathway_xml_ls = NULL) {
  lib_ps("MetaNet", "igraph", library = FALSE)
  if (is.null(pathway_xml_ls)) pathway_xml_ls <- load_something("pathway_xml_ls", with_new = FALSE)

  load_Pathway_htable(verbose = FALSE) -> Pathway_htable
  Pathway_htable$Pathway_id <- gsub("map", "ko", Pathway_htable$Pathway_id)
  Pathway_htable2 <- filter(Pathway_htable, !level2_name %in% c("Global and overview maps"))

  Pathway_htable2$Pathway_id -> pathway_ids

  pathway_net <- list()
  pathway_index_list <- list()
  for (path in pathway_ids) {
    message("Processing pathway: ", path)
    pathway_net[[path]] <- c_net_from_pathway_xml(pathway_xml_ls[[path]])
    if (is.null(pathway_net[[path]])) {
      message("Failed to create network for pathway: ", path)
      next
    }
    pathway_index_list[[path]] <- pathway_net_index(pathway_net[[path]])
  }
  do.call(rbind, pathway_index_list) -> all_pathway_index
  all_pathway_index
}
