#' Perform enrichment analysis
#'
#' This function performs KO enrichment analysis using the `clusterProfiler` package.
#'
#' @param ko_stat ko_stat dataframe from \code{\link[ReporterScore]{ko.test}}.
#' @param padj_threshold p.adjust threshold to determine whether a feature significant or not. p.adjust < padj_threshold, default: 0.05
#' @param logFC_threshold logFC threshold to determine whether a feature significant or not. abs(logFC)>logFC_threshold, default: NULL
#' @param add_mini add_mini when calculate the logFC. e.g (10+0.1)/(0+0.1), default 0.05*min(avg_abundance)
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @param feature one of "ko", "gene", "compound"
#' @param type "pathway" or "module" for default KOlist_file.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#' @family common_enrich
#' @return A data frame containing the enrichment results.
#' @export
#' @examples
#' ## use `enricher` from the `clusterProfiler` package.
#' if (requireNamespace("clusterProfiler")) {
#'   data("reporter_score_res")
#'   enrich_res <- KO_enrich(reporter_score_res)
#'   plot(enrich_res)
#' }
KO_enrich <- function(ko_stat, padj_threshold = 0.05,
                      logFC_threshold = NULL, add_mini = NULL, p.adjust.method = "BH",
                      type = c("pathway", "module")[1], feature = "ko",
                      modulelist = NULL, verbose = TRUE) {
  res.dt <- path2name <- path2ko <- sig_KO <- KOs <- NULL

  pcutils::lib_ps("clusterProfiler", library = FALSE)
  KO_enrich_internal(ko_stat, padj_threshold,
    logFC_threshold, add_mini, p.adjust.method,
    type, feature,
    modulelist, verbose,
    mode = 1
  )

  path2ko <- dplyr::filter(path2ko, KOs %in% res.dt$KO_id)
  e <- clusterProfiler::enricher(
    gene = sig_KO, TERM2GENE = path2ko, TERM2NAME = path2name,
    pAdjustMethod = p.adjust.method, pvalueCutoff = 1, qvalueCutoff = 1
  )

  if (verbose) pcutils::dabiao("`clusterProfiler::enricher` done")

  GO_res <- as.data.frame(e)
  GO_res <- rename(GO_res, "p.value" = "pvalue", "Significant_K_num" = "Count")
  GO_res <- GO_res[, c(seq_len(6), 9)]
  GO_res$Exist_K_num <- pcutils::strsplit2(GO_res$BgRatio, "/")[, 1] %>% as.numeric()
  GO_res$Exist_K_num <- as.integer(GO_res$Exist_K_num)
  GO_res$Significant_K_num <- as.integer(GO_res$Significant_K_num)
  class(GO_res) <- c("enrich_res", class(GO_res))
  attributes(GO_res)$method <- "enricher"
  attributes(GO_res)$type <- type
  return(GO_res)
}

KO_enrich_internal <- function(ko_stat, padj_threshold = 0.05,
                               logFC_threshold = NULL, add_mini = NULL, p.adjust.method = "BH",
                               type = c("pathway", "module")[1], feature = "ko",
                               modulelist = NULL, verbose = TRUE, mode = 1, weight = "logFC") {
  KO_id <- KOs <- p.adjust <- logFC <- NULL

  if (inherits(ko_stat, "reporter_score")) {
    reporter_res <- ko_stat
    ko_stat <- reporter_res$ko_stat
    modulelist <- reporter_res$modulelist
    if (is.character(modulelist)) {
      GOlist <- load_GOlist()
      modulelist <- eval(parse(text = modulelist))
    }
    type <- attributes(reporter_res$reporter_s)$type
  }
  res.dt <- ko_stat
  if (is.null(modulelist)) {
    modulelist <- get_modulelist(type, feature, verbose = verbose)
  }
  if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) stop("check your modulelist format!")

  if (!all(c("KO_id", "p.adjust") %in% colnames(res.dt))) {
    stop("check if p.adjust in your ko_stat dataframe!")
  }
  if ("origin_p.adjust" %in% colnames(res.dt)) {
    message("detect the origin_p.adjust, use the origin_p.adjust.")
    res.dt$p.adjust <- res.dt$origin_p.adjust
  }

  # modulelist内的KO才考虑
  path2ko <- pcutils::explode(modulelist[, c("id", "KOs")], 2, split = ",")
  res.dt <- dplyr::filter(res.dt, KO_id %in% path2ko$KOs)

  path2name <- modulelist[, c("id", "Description")]
  # set background
  # 这个跟指定universe的结果一致
  {
    path2ko <- dplyr::filter(path2ko, KOs %in% res.dt$KO_id)
  }

  if (!is.null(logFC_threshold)) {
    if (!"logFC" %in% colnames(res.dt)) {
      message("No logFC in the data.frame, calculate.")
      vs_group <- grep("average", colnames(res.dt), value = TRUE)
      if (length(vs_group) != 2) stop("logFC only available for two groups")
      tmp <- c(res.dt[, vs_group[1]], res.dt[, vs_group[2]])

      if (is.null(add_mini)) {
        add_mini <- min(tmp[tmp > 0]) * 0.05
      }
      res.dt$logFC <- log2((res.dt[, vs_group[2]] + add_mini) / (res.dt[, vs_group[1]] + add_mini))
    }
    sig_KO <- dplyr::filter(res.dt, p.adjust < padj_threshold, abs(logFC) > logFC_threshold) %>% dplyr::pull(KO_id)
  } else {
    sig_KO <- dplyr::filter(res.dt, p.adjust < padj_threshold) %>%
      dplyr::pull(KO_id) %>%
      unique()
  }

  if (length(sig_KO) < 1) {
    return(NULL)
  }


  assign("modulelist", modulelist, parent.frame())
  assign("sig_KO", sig_KO, parent.frame())
  assign("path2ko", path2ko, parent.frame())
  assign("path2name", path2name, parent.frame())
  assign("res.dt", res.dt, parent.frame())
}

#' Perform fisher's exact enrichment analysis
#' @inheritParams KO_enrich
#' @family common_enrich
#'
#' @return data.frame
#' @export
#'
#' @examples
#' ## use `fisher.test` from the `stats` package.
#' data("reporter_score_res")
#' fisher_res <- KO_fisher(reporter_score_res)
KO_fisher <- function(ko_stat, padj_threshold = 0.05,
                      logFC_threshold = NULL, add_mini = NULL, p.adjust.method = "BH",
                      type = c("pathway", "module")[1], feature = "ko",
                      modulelist = NULL, verbose = TRUE) {
  res.dt <- path2name <- path2ko <- sig_KO <- KOs <- Exist_K_num <- id <- p.value <- KO_id <- NULL

  KO_enrich_internal(ko_stat, padj_threshold,
    logFC_threshold, add_mini, p.adjust.method,
    type, feature,
    modulelist, verbose,
    mode = 2
  )
  sig_path_id <- dplyr::filter(path2ko, KOs %in% sig_KO) %>%
    dplyr::pull(id) %>%
    unique()

  modulelist <- dplyr::filter(modulelist, id %in% sig_path_id)

  sig_K_num <- length(sig_KO)
  nosig_K_num <- nrow(res.dt) - sig_K_num

  reps <- nrow(modulelist)

  lapply(seq_len(reps), \(i){
    tmp_kos <- strsplit(modulelist$KOs[i], ",")[[1]]

    z <- res.dt[res.dt$KO_id %in% tmp_kos, ] %>% dplyr::distinct(KO_id, .keep_all = TRUE)
    exist_KO <- nrow(z)
    significant_KO <- sum(z$KO_id %in% sig_KO)

    p_value <- stats::fisher.test(
      matrix(c(
        significant_KO, exist_KO - significant_KO,
        sig_K_num - significant_KO, nosig_K_num - (exist_KO - significant_KO)
      ), 2, 2, byrow = TRUE),
      alternative = "greater"
    )

    c(exist_KO, significant_KO, p_value$p.value)
  }) %>% do.call(rbind, .) -> res
  colnames(res) <- c("Exist_K_num", "Significant_K_num", "p.value")

  if (verbose) pcutils::dabiao("`fisher.test` done")

  fisher_res <- data.frame(
    ID = modulelist$id,
    Description = modulelist$Description,
    K_num = modulelist$K_num, res
  )

  fisher_res <- dplyr::filter(fisher_res, Exist_K_num > 0) %>% dplyr::arrange(p.value)
  fisher_res$p.adjust <- stats::p.adjust(fisher_res$p.value, method = p.adjust.method)

  fisher_res$Exist_K_num <- as.integer(fisher_res$Exist_K_num)
  fisher_res$Significant_K_num <- as.integer(fisher_res$Significant_K_num)
  class(fisher_res) <- c("enrich_res", class(fisher_res))
  attributes(fisher_res)$method <- "fisher.test"
  attributes(fisher_res)$type <- type
  return(fisher_res)
}


#' as enrich_res object
#'
#' @param gsea_res gsea_res from KO_gsea
#'
#' @return enrich_res object
#' @export
#' @rdname KO_enrich
#'
as.enrich_res <- function(gsea_res) {
  res <- data.frame(gsea_res, check.names = FALSE)
  class(res) <- c("enrich_res", class(res))
  res$Exist_K_num <- res$setSize
  attributes(res)$type <- attributes(gsea_res)$type
  res
}

#' Plot enrich_res
#'
#' @param enrich_res enrich_res object
#' @param mode plot style: 1~2
#' @param str_width default: 50
#' @param show_ID show pathway id
#' @param Pathway_description show KO description rather than KO id.
#' @param ... add
#' @param padj_threshold p.adjust threshold
#' @param facet_level facet plot if the type is "pathway" or "module"
#' @param facet_str_width str width for facet label
#' @param facet_anno annotation table for facet, two columns, first is level summary, second is pathway id.
#'
#' @family common_enrich
#' @return ggplot
#' @export
plot_enrich_res <- function(enrich_res, mode = 1, padj_threshold = 0.05,
                            show_ID = FALSE, Pathway_description = TRUE,
                            facet_level = FALSE, facet_anno = NULL, str_width = 50,
                            facet_str_width = 15, ...) {
  Description <- Significant_K_num <- order_value <- x <- fill <- size <- size_lab <- x_lab <- NULL
  GO <- pre_enrich_res(enrich_res, padj_threshold, show_ID, Pathway_description, facet_level, facet_anno)

  GO$Description <- factor(GO$Description, levels = dplyr::arrange(GO, -order_value) %>% dplyr::pull(Description))
  # 经典图
  if (mode == 1) {
    p <- ggplot(data = GO, aes(y = Description, x = x, fill = fill)) +
      geom_bar(stat = "identity", width = 0.7, position = "dodge")
  }
  if (mode == 2) {
    p <- ggplot(data = GO, aes(y = Description, x = x, fill = fill, size = size)) +
      geom_point(shape = 21, position = position_dodge(width = 0.7)) +
      scale_size(name = size_lab, range = c(3, 8))
  }

  if (is.numeric(GO$fill)) {
    p <- p + scale_fill_gradient(name = "p.adjust", low = "red", high = "blue", limits = c(0, padj_threshold))
  } else if (any(c("NES", "GSA.scores") %in% colnames(GO))) {
    p <- p + scale_fill_manual(name = "Group", values = c("red", "blue"))
  } else {
    p <- p + scale_fill_discrete(name = "Group")
  }

  p <- p + scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width)) +
    theme_bw() + labs(y = NULL, x = x_lab)

  if ("facet_level" %in% colnames(GO)) {
    p <- p + facet_grid(facet_level ~ ., scales = "free_y", space = "free", labeller = label_wrap_gen(facet_str_width)) +
      theme(
        strip.text.y = element_text(angle = 0),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", color = "black")
      )
  }
  return(p)
}

pre_enrich_res <- function(enrich_res, padj_threshold = 0.05,
                           show_ID = FALSE, Pathway_description = TRUE,
                           facet_level = FALSE, facet_anno = NULL) {
  flag <- FALSE
  if (inherits(enrich_res, "enrich_res")) {
    GO <- enrich_res
  } else if (is.list(enrich_res) & all(vapply(enrich_res, \(i)inherits(i, "enrich_res"), logical(1)))) {
    multi_enrich_res <- enrich_res
    if (is.null(names(multi_enrich_res))) names(multi_enrich_res) <- paste0("Res", seq_along(multi_enrich_res))
    GO <- lapply(
      names(multi_enrich_res),
      \(i){
        data.frame(multi_enrich_res[[i]][, c("ID", "Description", "Exist_K_num", "Significant_K_num", "p.adjust")], Group = i, row.names = NULL)
      }
    ) %>%
      do.call(rbind, .)
    attributes(GO)$type <- attributes(multi_enrich_res[[1]])$type
    flag <- TRUE
  } else {
    GO <- enrich_res
  }

  GO <- dplyr::filter(GO, p.adjust <= padj_threshold)
  if (nrow(GO) < 1) stop("No pathway p.adjst less than ", padj_threshold)

  if (facet_level) {
    tmpdf <- get_facet_anno(GO, facet_anno)
    if (!is.null(tmpdf)) {
      GO <- dplyr::left_join(GO, tmpdf, by = c("ID"))
    }
  }

  if (show_ID) GO$Description <- paste0(GO$ID, ": ", GO$Description)
  if (!Pathway_description) GO$Description <- GO$ID

  if (flag) {
    GO$order_value <- -GO$p.adjust
    GO$x <- -log(GO$p.adjust)
    GO$fill <- GO$Group
    x_lab <- "-log(p.adjust)"
  } else if ("NES" %in% colnames(GO)) {
    GO$order_value <- -GO$NES
    GO$x <- GO$NES
    GO$Group <- ifelse(GO$NES > 0, "Up", "Down")
    GO$fill <- GO$Group
    x_lab <- "NES"
  } else if ("GSA.scores" %in% colnames(GO)) {
    GO$order_value <- -GO$GSA.scores
    GO$x <- GO$GSA.scores
    GO$Group <- ifelse(GO$`GSA.scores` > 0, "Up", "Down")
    GO$fill <- GO$Group
    x_lab <- "GSA.scores"
  } else {
    GO$order_value <- -GO$p.adjust
    if ("Significant_K_num" %in% colnames(GO)) {
      GO$x <- GO$Significant_K_num / GO$Exist_K_num
      x_lab <- "Significant_K_Ratio"
    } else {
      GO$x <- GO$Exist_K_num
      x_lab <- "Exist_K_num"
    }
    GO$fill <- GO$p.adjust
  }

  if ("Significant_K_num" %in% colnames(GO)) {
    size_lab <- "Significant_K_num"
  } else {
    size_lab <- "Exist_K_num"
  }

  GO$size <- GO[, size_lab, drop = TRUE]
  assign("size_lab", size_lab, parent.frame())
  assign("x_lab", x_lab, parent.frame())
  return(GO)
}

#' Plot enrich_res
#' @param x enrich_res object
#' @rdname plot_enrich_res
#'
#' @return ggplot
#' @exportS3Method
#' @method plot enrich_res
plot.enrich_res <- function(x, mode = 1, padj_threshold = 0.05,
                            show_ID = FALSE, Pathway_description = TRUE,
                            facet_level = FALSE, facet_anno = NULL, str_width = 50, facet_str_width = 15, ...) {
  plot_enrich_res(
    x, mode, padj_threshold, show_ID,
    Pathway_description,
    facet_level, facet_anno, str_width, facet_str_width, ...
  )
}

#' Perform gene set enrichment analysis
#'
#' @param weight the metric used for ranking, default: logFC
#' @inheritParams KO_enrich
#' @family common_enrich
#' @export
#'
#' @return DOSE object
#' @examples
#' message("The following example require some time to run:")
#' \donttest{
#' ## use `GSEA` from the `clusterProfiler` package.
#' if (requireNamespace("clusterProfiler")) {
#'   data("reporter_score_res")
#'   gsea_res <- KO_gsea(reporter_score_res, p.adjust.method = "none")
#'   enrichplot::gseaplot(gsea_res, geneSetID = data.frame(gsea_res)$ID[1])
#'   gsea_res_df <- as.enrich_res(gsea_res)
#'   plot(gsea_res_df)
#' }
#' }
KO_gsea <- function(ko_stat, weight = "logFC", add_mini = NULL,
                    p.adjust.method = "BH",
                    type = c("pathway", "module")[1], feature = "ko",
                    modulelist = NULL, verbose = TRUE) {
  res.dt <- path2name <- path2ko <- sig_KO <- KOs <- Exist_K_num <- id <- p.value <- logFC <- NULL

  pcutils::lib_ps("clusterProfiler", library = FALSE)
  KO_enrich_internal(ko_stat,
    padj_threshold = 1,
    logFC_threshold = NULL, add_mini, p.adjust.method,
    type, feature,
    modulelist, verbose, mode = 3, weight
  )
  if (!weight %in% colnames(res.dt)) {
    message("Use logFC as the weight.")
    if (!"logFC" %in% colnames(res.dt)) {
      message("No logFC in the data.frame, calculate.")
      vs_group <- grep("average", colnames(res.dt), value = TRUE)
      if (length(vs_group) != 2) stop("logFC only available for two groups")
      tmp <- c(res.dt[, vs_group[1]], res.dt[, vs_group[2]])

      if (is.null(add_mini)) {
        add_mini <- min(tmp[tmp > 0]) * 0.05
      }
      res.dt$logFC <- log2((res.dt[, vs_group[2]] + add_mini) / (res.dt[, vs_group[1]] + add_mini))
    }
    res.dt_sort <- res.dt %>% dplyr::arrange(-logFC)
    kos <- res.dt_sort$logFC
    names(kos) <- res.dt_sort$KO_id
  } else {
    kos <- res.dt[, weight, drop = TRUE]
    names(kos) <- res.dt$KO_id
    kos <- sort(kos, decreasing = TRUE)
  }
  kos <- na.omit(kos)
  e <- clusterProfiler::GSEA(kos,
    TERM2GENE = path2ko, TERM2NAME = path2name, verbose = FALSE,
    pvalueCutoff = 1, pAdjustMethod = p.adjust.method
  )
  if (verbose) pcutils::dabiao("`clusterProfiler::GSEA` done")
  attributes(e)$type <- type
  return(e)
}

pre_rs <- function(reporter_res, mode = 1, verbose = TRUE) {
  Exist_K_num <- NULL
  stopifnot(inherits(reporter_res, "reporter_score"))

  modulelist <- reporter_res$modulelist
  if (is.null(modulelist)) stop("no modulelist")
  if (is.character(modulelist)) {
    GOlist <- load_GOlist()
    modulelist <- eval(parse(text = modulelist))
  }
  if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) stop("check your modulelist format!")

  group <- reporter_res$group
  metadata <- reporter_res$metadata
  kodf <- reporter_res$kodf
  if (verbose) pcutils::dabiao("Checking group")
  if (!is.null(metadata)) {
    if (length(group) != 1) stop("'group' should be one column name of metadata when metadata exsit!")
    idx <- rownames(metadata) %in% colnames(kodf)
    metadata <- metadata[idx, , drop = FALSE]
    kodf <- kodf[, rownames(metadata), drop = FALSE]
    if (verbose) message(nrow(metadata), " samples are matched for next step.")
    if (length(idx) < 2) stop("too less common samples")
    sampFile <- data.frame(group = metadata[, group], row.names = row.names(metadata))
  } else {
    if (length(group) != ncol(kodf)) stop("'group' length should equal to columns number of kodf when metadata is NULL!")
    sampFile <- data.frame(row.names = colnames(kodf), group = group)
  }
  modulelist$Exist_K_num <- vapply(transform_modulelist(modulelist), \(i){
    sum(rownames(kodf) %in% i)
  }, numeric(1))

  modulelist <- dplyr::filter(modulelist, Exist_K_num > 0)

  envir <- parent.frame()

  assign("sampFile", sampFile, envir)
  assign("modulelist", modulelist, envir)
  if (mode == 1) {
    if (verbose) pcutils::dabiao("Removing all-zero rows: ", sum(rowSums(abs(kodf)) == 0))
    kodf <- kodf[rowSums(abs(kodf)) > 0, ]
    assign("kodf", kodf, envir)
  } else if (mode == 2) {
    ko_stat <- reporter_res$ko_stat
    if (!all(c("KO_id", "p.adjust") %in% colnames(ko_stat))) {
      stop("check if p.adjust in your ko_stat dataframe!")
    }
    if ("origin_p.adjust" %in% colnames(ko_stat)) {
      message("detect the origin_p.adjust, use the origin_p.adjust.")
      ko_stat$p.adjust <- ko_stat$origin_p.adjust
    }
    assign("ko_stat", ko_stat, envir)
  }
}

#' Perform gene set analysis
#'
#' @param reporter_res reporter_res
#' @param method Problem type: "quantitative" for a continuous parameter; "Two class unpaired" ; "Survival" for censored survival outcome; "Multiclass" : more than 2 groups, coded 1,2,3...; "Two class paired" for paired outcomes, coded -1,1 (first pair), -2,2 (second pair), etc
#' @param p.adjust.method "BH"
#' @param verbose TRUE
#' @param perm 1000
#' @param ... additional parameters to \code{\link[GSA]{GSA}}
#'
#' @return enrich_res object
#' @export
#'
#' @family common_enrich
#' @examples
#' \donttest{
#' ## use `GSA` from the `GSA` package.
#' if (requireNamespace("GSA")) {
#'   data("reporter_score_res")
#'   gsa_res <- KO_gsa(reporter_score_res, p.adjust.method = "none", perm = 200)
#'   plot(gsa_res)
#' }
#' }
KO_gsa <- function(reporter_res, method = "Two class unpaired", p.adjust.method = "BH", verbose = TRUE, perm = 1000, ...) {
  kodf <- modulelist <- kodf <- modulelist <- p.value <- NULL
  pre_rs(reporter_res, mode = 1, verbose = verbose)

  tkodf <- pcutils::t2(kodf)
  sampFile$group <- as.numeric(as.factor(sampFile$group))

  genesets <- transform_modulelist(modulelist)
  lib_ps("GSA", library = FALSE)
  GSA.obj <- GSA::GSA(as.matrix(kodf), sampFile$group,
    genenames = rownames(kodf), genesets = genesets,
    resp.type = method, nperms = perm, ...
  )

  gsa_res <- as.data.frame(GSA.obj[c("GSA.scores", "pvalues.lo", "pvalues.hi")])
  gsa_res$p.value <- apply(gsa_res, 1, \(i)ifelse(i[1] > 0, i[3], i[2]))
  gsa_res$p.adjust <- p.adjust(gsa_res$p.value, method = p.adjust.method)

  path_res <- data.frame(
    row.names = modulelist$id, ID = modulelist$id,
    Description = modulelist$Description,
    K_num = modulelist$K_num,
    Exist_K_num = modulelist$Exist_K_num,
    gsa_res
  )
  path_res <- dplyr::arrange(path_res, p.value)
  class(path_res) <- c("enrich_res", class(path_res))
  path_res
}


#' Perform Gene Set Variation Analysis
#'
#' @param reporter_res reporter_res
#' @param verbose verbose
#' @param method see \code{\link{ko.test}}
#' @param p.adjust.method p.adjust.method
#' @param ... additional parameters to \code{\link[GSVA]{gsva}}
#'
#' @return enrich_res
#' @export
#' @family common_enrich
#' @examples
#' \donttest{
#' ## use `gsva` from the `GSVA` package.
#' if (requireNamespace("GSVA")) {
#'   data("reporter_score_res")
#'   gsva_res <- KO_gsva(reporter_score_res, p.adjust.method = "none")
#' }
#' }
KO_gsva <- function(reporter_res, verbose = TRUE, method = "wilcox.test", p.adjust.method = "BH", ...) {
  kodf <- sampFile <- modulelist <- p.value <- NULL
  pre_rs(reporter_res, mode = 1, verbose = verbose)

  geneSets <- transform_modulelist(modulelist)
  lib_ps("GSVA", library = FALSE)
  gsva_es <- GSVA::gsva(as.matrix(kodf), geneSets, mx.diff = 1, verbose = verbose, ...)

  gsva_res <- ko.test(gsva_es,
    group = sampFile$group, method = method,
    p.adjust.method1 = p.adjust.method, verbose = verbose
  )

  gsva_res <- data.frame(
    modulelist[
      match(rownames(gsva_res), modulelist$id),
      c("id", "Description", "K_num", "Exist_K_num")
    ],
    gsva_res[, c("p.value", "p.adjust")]
  )
  colnames(gsva_res)[1] <- "ID"
  gsva_res <- dplyr::arrange(gsva_res, p.value)
  gsva_res
}


#' Perform Simultaneous Enrichment Analysis
#'
#'
#' @param reporter_res The input reporter result.
#' @param verbose If TRUE, print verbose messages. Default is TRUE.
#' @param ... Additional parameters to be passed to \code{\link[rSEA]{SEA}} function.
#'
#' @return enrich_res
#'
#' @export
#' @family common_enrich
#' @examples
#' \donttest{
#' ## use `SEA` from the `rSEA` package.
#' if (requireNamespace("rSEA")) {
#'   data("reporter_score_res")
#'   sea_res <- KO_sea(reporter_score_res, verbose = TRUE)
#' }
#' }
KO_sea <- function(reporter_res, verbose = TRUE, ...) {
  ko_stat <- kodf <- sampFile <- modulelist <- p.value <- NULL
  pre_rs(reporter_res, mode = 2, verbose = verbose)

  pathlist <- transform_modulelist(modulelist)
  lib_ps("rSEA", library = FALSE)
  sea_res <- rSEA::SEA(ko_stat$p.adjust, ko_stat$KO_id, pathlist = pathlist, ...)
  sea_res <- data.frame(modulelist[, c("id", "Description", "K_num", "Exist_K_num")], sea_res[, -seq_len(3)])
  colnames(sea_res)[1] <- "ID"
  sea_res$p.adjust <- apply(sea_res, 1, \(i)min(i[8:9])) %>% as.numeric()
  sea_res <- dplyr::arrange(sea_res, p.adjust)
  sea_res
}

#' Perform Significance Analysis of Function and Expression
#'
#'
#' @param reporter_res The input reporter result.
#' @param verbose If TRUE, print verbose messages. Default is TRUE.
#' @param perm The number of permutations. Default is 1000.
#' @param C.matrix The contrast matrix. Default is NULL, and it will be generated from the module list.
#' @param p.adjust.method Method for p-value adjustment. Default is "BH".
#' @param ... Additional parameters to be passed to \code{\link[safe]{safe}} function.
#'
#' @return A data frame containing SAFE results for KO enrichment.
#'
#' @export
#' @family common_enrich
#' @examples
#' \donttest{
#' ## use `safe` from the `safe` package.
#' if (requireNamespace("safe")) {
#'   data("reporter_score_res")
#'   safe_res <- KO_safe(reporter_score_res,
#'     verbose = TRUE,
#'     perm = 200, p.adjust.method = "none"
#'   )
#' }
#' }
KO_safe <- function(reporter_res, verbose = TRUE, perm = 1000,
                    C.matrix = NULL, p.adjust.method = "BH", ...) {
  kodf <- sampFile <- modulelist <- p.value <- NULL
  pre_rs(reporter_res, mode = 1, verbose = verbose)

  if (is.null(C.matrix)) C.matrix <- transform_modulelist(modulelist, mode = 3) %>% as.matrix()
  C.matrix <- C.matrix[intersect(rownames(kodf), rownames(C.matrix)), ]
  kodf <- kodf[intersect(rownames(kodf), rownames(C.matrix)), ]
  lib_ps("safe", library = FALSE)
  safe_res <- safe::safe(X.mat = kodf, y.vec = sampFile$group, C.mat = C.matrix, Pi.mat = perm, ...)
  safe_res <- data.frame(modulelist[, c("id", "Description", "K_num", "Exist_K_num")],
    stat = safe_res@global.stat,
    p.value = safe_res@global.pval,
    error = safe_res@global.error
  )
  colnames(safe_res)[1] <- "ID"
  safe_res$p.adjust <- p.adjust(safe_res$p.value, method = p.adjust.method)
  safe_res <- dplyr::arrange(safe_res, p.value)
  safe_res
}

#' Perform Pathway Analysis with Down-weighting of Overlapping Genes (PADOG)
#'
#' @param reporter_res The input reporter result.
#' @param verbose If TRUE, print verbose messages. Default is TRUE.
#' @param p.adjust.method Method for p-value adjustment. Default is "BH".
#' @param ... Additional parameters to be passed to \code{\link[PADOG]{padog}} function.
#' @param perm The number of permutations. Default is 1000.
#'
#' @return A data frame containing PADOG results for KO enrichment.
#' @return A data frame with columns "ID," "Description," "K_num," "Exist_K_num," "p.value," and "p.adjust."
#'
#' @export
#' @family common_enrich
#' @examples
#' \donttest{
#' ## use `PADOG` from the `PADOG` package.
#' if (requireNamespace("PADOG")) {
#'   data("reporter_score_res")
#'   padog_res <- KO_padog(reporter_score_res,
#'     verbose = TRUE,
#'     perm = 200, p.adjust.method = "none"
#'   )
#' }
#' }
KO_padog <- function(reporter_res, verbose = TRUE, perm = 1000, p.adjust.method = "BH", ...) {
  kodf <- sampFile <- modulelist <- p.value <- NULL
  pre_rs(reporter_res, mode = 1, verbose = verbose)

  lib_ps("PADOG", library = FALSE)
  geneSets <- transform_modulelist(modulelist)
  padog_res <- PADOG::padog(as.matrix(kodf),
    group = factor(sampFile$group, labels = c("c", "d")),
    gslist = geneSets, NI = perm, ...
  )

  padog_res <- data.frame(
    modulelist[
      match(rownames(padog_res), modulelist$id),
      c("id", "Description", "K_num", "Exist_K_num")
    ],
    padog_res[, -seq_len(2)],
    row.names = NULL
  )
  colnames(padog_res)[c(1, 9)] <- c("ID", "p.value")
  padog_res$p.adjust <- p.adjust(padog_res$p.value, method = p.adjust.method)
  padog_res <- dplyr::arrange(padog_res, p.value)
  padog_res
}
