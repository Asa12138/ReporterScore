reporter_color <- c("#e31a1c", "#1f78b4", "#b15928", "#810f7c", "#FFF021", "#b2df8a")
reporter_theme <- {
    ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
            axis.text = element_text(color = "black"),
            plot.margin = grid::unit(rep(0.5, 4), "lines"),
            strip.background = ggplot2::element_rect(fill = NA)
        )
}

#' Plot the reporter_res
#'
#' @param reporter_res result of `get_reporter_score` or `reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param y_text_size y_text_size
#' @param str_width str_width to wrap
#' @param mode 1~2 plot style.
#' @param show_ID show pathway id
#' @param Pathway_description show KO description rather than KO id.
#' @param facet_level facet plot if the type is "pathway" or "module"
#' @param facet_str_width str width for facet label
#' @param facet_anno annotation table for facet, two columns, first is level summary, second is pathway id.
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @aliases plot_report_bar
#' @examples
#' data("reporter_score_res")
#' plot_report(reporter_score_res, rs_threshold = c(2.5, -2.5), y_text_size = 10, str_width = 40)
plot_report <- function(reporter_res, rs_threshold = 1.64, mode = 1, y_text_size = 13, str_width = 100, show_ID = FALSE,
                        Pathway_description = TRUE, facet_level = FALSE, facet_anno = NULL, facet_str_width = 15) {
    Group <- Description <- ReporterScore <- Exist_K_num <- NULL
    reporter_res2=cols1=title=breaks=NULL

    reporter_res <- pre_reporter_res(reporter_res, rs_threshold)
    flag <- attributes(reporter_res)$flag
    filter_report(reporter_res, rs_threshold)

    if (flag) {
        reporter_res2$Group <- reporter_res2$Cluster
        cols1 <- setNames(
            pcutils::get_cols(length(unique(reporter_res2$Group))),
            unique(reporter_res2$Group)
        )
    }

    if (facet_level) {
        tmpdf <- get_facet_anno(reporter_res, facet_anno)
        if (!is.null(tmpdf)) {
            reporter_res2 <- dplyr::left_join(reporter_res2, tmpdf, by = c("ID"))
        }
    }

    reporter_res2 <- reporter_res2[stats::complete.cases(reporter_res2), ]
    # rownames(reporter_res2)=reporter_res2$ID

    if (show_ID) reporter_res2$Description <- paste0(reporter_res2$ID, ": ", reporter_res2$Description)
    if (!Pathway_description) reporter_res2$Description <- reporter_res2$ID

    if (mode == 1) {
        if (flag) {
            reporter_res2$Description <- factor(reporter_res2$Description,
                levels = dplyr::arrange(reporter_res2, Group, ReporterScore) %>%
                    dplyr::pull(Description) %>% unique()
            )
            p <- ggplot(reporter_res2, aes(ReporterScore, Description, fill = Group))
        } else {
            p <- ggplot(reporter_res2, aes(ReporterScore, stats::reorder(Description, ReporterScore), fill = Group))
        }
        p <- p + geom_bar(stat = "identity", position = "dodge") +
            scale_fill_manual(values = cols1) +
            theme_light()
    }
    if (mode == 2) {
        if (flag) {
            reporter_res2$Description <- factor(reporter_res2$Description,
                levels = dplyr::arrange(reporter_res2, Group, ReporterScore) %>%
                    dplyr::pull(Description) %>% unique()
            )
            p <- ggplot(reporter_res2, aes(ReporterScore, Description, size = Exist_K_num, fill = Exist_K_num))
        } else {
            p <- ggplot(reporter_res2, aes(ReporterScore, stats::reorder(Description, ReporterScore), size = Exist_K_num, fill = Exist_K_num))
        }
        p <- p +
            geom_point(shape = 21, position = "dodge") +
            scale_fill_gradient(low = "#FF000033", high = "red", guide = "legend") + theme_light()
    }

    p <- p + labs(y = "") +
        scale_y_discrete(labels = label_wrap_gen(width = str_width)) +
        scale_x_continuous(breaks = breaks) +
        theme(
            axis.text.x = element_text(colour = "black", size = 13),
            axis.text.y = element_text(colour = "black", size = y_text_size)
        )
    if (facet_level) {
        p <- p + facet_grid(facet_level ~ ., scales = "free_y", space = "free", labeller = label_wrap_gen(facet_str_width)) +
            theme(
                strip.text.y = element_text(angle = 0),
                strip.background = element_rect(fill = "grey90"),
                strip.text = element_text(face = "bold", color = "black")
            )
    }
    if (attributes(reporter_res)$mode == "directed" & is.null(attributes(reporter_res)$pattern)) {
        if (any(reporter_res2$ReporterScore > rs_threshold[2])) p <- p + geom_vline(xintercept = rs_threshold[2], linetype = 2)
        if (any(reporter_res2$ReporterScore < rs_threshold[1])) p <- p + geom_vline(xintercept = rs_threshold[1], linetype = 2)
    } else {
        p <- p + geom_vline(xintercept = rs_threshold[2], linetype = 2)
    }
    # if(length(attributes(reporter_res)$vs_group)==2)p=p+labs(title = paste(attributes(reporter_res)$vs_group,collapse = " vs "))
    p <- p + labs(title = title)
    return(p)
}

pre_reporter_res <- function(reporter_res, rs_threshold) {
    reporter_res2=cols1=title=breaks=NULL
    if (inherits(reporter_res, "reporter_score")) {
        reporter_res <- reporter_res$reporter_s
    }
    if (is.data.frame(reporter_res)) {}
    attributes(reporter_res)$flag <- FALSE
    if (inherits(reporter_res, "rs_by_cm")) {
        rsa_cm_res <- reporter_res
        ncluster <- sum(grepl("Cluster", names(rsa_cm_res)))
        clusters_name <- grep("Cluster", names(rsa_cm_res), value = TRUE)
        reporter_res <- lapply(
            clusters_name,
            \(i){
                data.frame(rsa_cm_res[[i]]$reporter_s, Cluster = i, row.names = NULL)
            }
        ) %>%
            do.call(rbind, .)
        attributes(reporter_res) <- pcutils::update_param(attributes(rsa_cm_res[["Cluster1"]]$reporter_s), attributes(reporter_res))
        attributes(reporter_res)$flag <- TRUE
    } else if (all(vapply(reporter_res, \(i){
        inherits(i, "reporter_score")
    }, logical(1)))) {
        multi_reporter_res <- reporter_res
        if (is.null(names(multi_reporter_res))) names(multi_reporter_res) <- paste0("Res", seq_along(multi_reporter_res))
        ncluster <- length(multi_reporter_res)
        clusters_name <- names(multi_reporter_res)

        reporter_res <- lapply(
            clusters_name,
            \(i){
                filter_report(multi_reporter_res[[i]]$reporter_s, rs_threshold)
                data.frame(reporter_res2, Cluster = i, row.names = NULL)
            }
        ) %>%
            do.call(rbind, .)
        attributes(reporter_res) <- pcutils::update_param(attributes(multi_reporter_res[[clusters_name[1]]]$reporter_s), attributes(reporter_res))

        attributes(reporter_res)$flag <- TRUE
    }
    return(reporter_res)
}

filter_report <- function(reporter_res, rs_threshold) {
    reporter_res <- na.omit(reporter_res)
    if (length(rs_threshold) == 1) rs_threshold <- c(rs_threshold, -rs_threshold)

    vs_group <- attributes(reporter_res)$vs_group

    rs_threshold <- sort(rs_threshold)

    if ((attributes(reporter_res)$mode == "directed") & is.null(attributes(reporter_res)$pattern)) {
        reporter_res2 <- reporter_res[(reporter_res$ReporterScore >= rs_threshold[2]) | (reporter_res$ReporterScore <= rs_threshold[1]), ]
        if (nrow(reporter_res2) < 1) stop("No pathway left.")
        if (length(vs_group) == 2) {
            reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0,
                paste0("Enriched in ", vs_group[2]),
                paste0("Enriched in ", vs_group[1])
            )
            cols1 <- setNames(c("P" = "orange", "N" = "seagreen"), paste0("Enriched in ", vs_group[2:1]))
        } else {
            reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0, "Increase", "Decrease")
            cols1 <- setNames(c("P" = "orange", "N" = "seagreen"), c("Increase", "Decrease"))
        }
        title <- paste0(vs_group, collapse = "/ ")
        breaks <- c(scales::breaks_extended(3)(range(c(0, reporter_res2$ReporterScore))), rs_threshold)
    } else {
        reporter_res2 <- reporter_res[reporter_res$ReporterScore >= rs_threshold[2], ]
        if (nrow(reporter_res2) < 1) stop("No pathway left.")
        reporter_res2$Group <- "Significant"
        cols1 <- c("Significant" = "red2")
        title <- paste0(vs_group, collapse = "/")
        breaks <- c(scales::breaks_extended(3)(range(c(0, reporter_res2$ReporterScore))), rs_threshold[2])
    }

    envir <- parent.frame()
    assign("reporter_res2", reporter_res2, envir)
    assign("cols1", cols1, envir)
    assign("title", title, envir)
    assign("breaks", breaks, envir)
    assign("rs_threshold", rs_threshold, envir)
    # return(list(reporter_res2=reporter_res2,cols1=cols1,title=title,breaks=breaks))
}

get_facet_anno <- function(reporter_res, facet_anno, mode = c("bar", "circle")[1]) {
    if (mode == "bar") {
        if (!is.null(facet_anno)) {
            tmpdf <- facet_anno
        } else {
            if (is.null(attributes(reporter_res)$type)) {
                warning("No attributes(reporter_res)$type found.")
                envir <- parent.frame()
                assign("facet_level", FALSE, envir)
                return(NULL)
            } else if (attributes(reporter_res)$type == "pathway") {
                load_Pathway_htable(envir = environment())
                tmpdf <- Pathway_htable[, c("level1_name", "Pathway_id")]
            } else if (attributes(reporter_res)$type == "module") {
                load_Module_htable(envir = environment())
                tmpdf <- Module_htable[c("module2_name", "Module_id")]
            } else if (attributes(reporter_res)$type == "ALL") {
                tmpdf <- reporter_res[c("ONT", "ID")]
            } else {
                # other organisms
                load_Pathway_htable(envir = environment(), verbose = FALSE)
                tmpdf <- Pathway_htable[, c("level1_name", "Pathway_id")]
                tmpdf$Pathway_id <- gsub("map", attributes(reporter_res)$type, tmpdf$Pathway_id)
            }
        }
        colnames(tmpdf) <- c("facet_level", "ID")
        return(tmpdf)
    } else if (mode == "circle") {
        if (!is.null(facet_anno)) {
            tmpdf <- facet_anno
            tmpdf$add_col <- "add"
            node_name <- colnames(facet_anno)[ncol(facet_anno)]
        } else {
            if (is.null(attributes(reporter_res)$type)) {
                stop("No attributes(reporter_res)$type found.")
            } else if (attributes(reporter_res)$type == "pathway") {
                load_Pathway_htable(envir = environment())
                tmpdf <- Pathway_htable
                node_name <- "Pathway"
            } else if (attributes(reporter_res)$type == "module") {
                load_Module_htable(envir = environment())
                tmpdf <- Module_htable
                node_name <- "Module"
            } else if (attributes(reporter_res)$type == "ALL") {
                tmpdf <- reporter_res[c("ONT", "ID", "Description")]
                node_name <- "GO term"
            } else {
                load_Pathway_htable(envir = environment())
                tmpdf <- Pathway_htable
                tmpdf$Pathway_id <- gsub("map", attributes(reporter_res)$type, tmpdf$Pathway_id)
                node_name <- "Pathway"
            }
        }
        attributes(tmpdf)$node_name <- node_name
        return(tmpdf)
    }
}

#' Plot the reporter_res as circle_packing
#'
#' @param reporter_res result of `get_reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param str_width str_width to wrap
#' @param mode 1~2 plot style.
#' @param Pathway_description show KO description rather than KO id.
#' @param facet_anno annotation table for facet, more two columns, last is pathway name, last second is pathway id.
#' @param show_ID show pathway id
#' @param show_level_name show the level name?
#' @param show_tip_label show the tip label?
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_report_circle_packing(reporter_score_res, rs_threshold = c(2, -2), str_width = 40)
plot_report_circle_packing <- function(reporter_res, rs_threshold = 1.64, mode = 2, facet_anno = NULL,
                                       show_ID = FALSE, Pathway_description = TRUE,
                                       str_width = 10, show_level_name = "all", show_tip_label = TRUE) {
    reporter_res2=cols1=title=breaks=ID=NULL
    reporter_res <- pre_reporter_res(reporter_res, rs_threshold)
    flag <- attributes(reporter_res)$flag
    filter_report(reporter_res, rs_threshold)

    if (flag) {
        if (mode == 1) message("mode can just set to 2 when reporter_res is a `rs_by_cm` object")
        mode <- 2
        reporter_res2$Group <- reporter_res2$Cluster
        cols1 <- setNames(
            pcutils::get_cols(length(unique(reporter_res2$Group))),
            unique(reporter_res2$Group)
        )
    }

    tmpdf <- get_facet_anno(reporter_res, facet_anno, mode = "circle")
    node_name <- attributes(tmpdf)$node_name
    reporter_res2 <- reporter_res2[stats::complete.cases(reporter_res2), ]

    reporter_res2 <- dplyr::distinct(reporter_res2, ID, .keep_all = TRUE)
    rownames(reporter_res2) <- reporter_res2$ID
    # prepare the hierarchy table
    f_tax <- data.frame(tmpdf[-ncol(tmpdf)], row.names = tmpdf[, ncol(tmpdf) - 1, drop = TRUE], check.names = FALSE)
    colnames(f_tax)[ncol(f_tax)] <- node_name

    tmp_weight <- reporter_res2[, c("ID", "ReporterScore")]
    tmp_weight$weight <- abs(tmp_weight$ReporterScore)

    f_tax2 <- dplyr::right_join(f_tax, tmp_weight[, c("ID", "weight")], by = setNames("ID", node_name))

    # prepare the anno table
    anno <- reporter_res2[, c("ID", "Description", "ReporterScore", "Group")]
    if (show_ID) anno$Description <- paste0(anno$ID, ": ", anno$Description)
    if (Pathway_description) {
        anno$Description <- anno$Description
    }

    p <- pcutils::my_circle_packing(f_tax2, anno,
        mode = mode, Group = "Group", Score = "ReporterScore", label = "Description",
        show_level_name = show_level_name, show_tip_label = show_tip_label, str_width = str_width
    )

    if (mode == 1) {
        if (length(cols1) == 2) {
            p <- p + scale_fill_gradient2(low = "seagreen", high = "orange", na.value = NA, name = "ReporterScore")
        } else {
            p <- p + scale_fill_gradient(low = "white", high = "red2", na.value = NA, name = "ReporterScore", limits = c(0, max(tmp_weight$weight)))
        }
    }
    if (mode == 2) {
        p <- p + scale_fill_manual(values = pcutils::add_alpha(cols1, 0.8), na.translate = FALSE)
    }
    p + scale_color_manual(values = get_cols(ncol(f_tax2), reporter_color)) + coord_fixed()
}


#' Plot the significance of pathway
#' @param reporter_res result of `get_reporter_score` or `reporter_score`
#' @param map_id the pathway or module id
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_significance(reporter_score_res, map_id = c("map05230", "map03010"))
plot_significance <- function(reporter_res, map_id) {
    value=ID=NULL
    if (inherits(reporter_res, "reporter_score")) reporter_res <- reporter_res$reporter_s

    lib_ps("ggrepel", library = FALSE)
    if (!is.null(attributes(reporter_res)$perm)) {
        perm <- attributes(reporter_res)$perm
    } else {
        perm <- 1000
    }
    reporter_res2 <- na.omit(dplyr::filter(reporter_res, ID %in% map_id))
    if (nrow(reporter_res2) < 1) {
        return(NULL)
    }
    reporter_res2 <- reporter_res2[, c("ID", "BG_Mean", "BG_Sd", "Z_score", "p.value", "ReporterScore")]
    # 生成正态分布数据
    # set.seed(123)  # 设置随机种子以确保结果可重复
    df <- apply(reporter_res2, 1, \(i)rnorm(perm, mean = as.numeric(i[2]), sd = as.numeric(i[3])))
    df2 <- reshape2::melt(as.data.frame(df), measure.vars = colnames(df), variable.name = "ID")

    # 绘制直方图
    p <- ggplot(df2, aes(x = value)) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        labs(title = paste0("Background distribution of Pathway Z-score; Permutation: ", perm), x = "Pathway Z-score", y = "") +
        facet_grid(~ID)
    p

    # 找到最高柱子的高度
    # max_height <- max(ggplot_build(p)$data[[1]]$count)
    reddf <- data.frame(reporter_res2, y = perm / 30)

    # 添加红色小棒子
    p <- p +
        geom_segment(data = reddf, aes(x = Z_score, xend = Z_score, y = 0, yend = y), color = "red", size = 1) +
        geom_point(data = reddf, aes(x = Z_score, y = y), color = "red", size = 4) +
        ggrepel::geom_text_repel(
            data = reddf,
            aes(x = Z_score, y = y * 1.4, label = paste0("ReporterScore: ", round(ReporterScore, 2), "\np-value: ", p.value)),
            color = "red", size = 3
        )
    return(p)
}

#' plot the Z-score of features distribution
#'
#' @param reporter_res result of `reporter_score`
#' @param map_id the pathway or module id
#' @param text_size text_size=4
#' @param rug_length rug_length=0.04
#' @param text_position text_position, e.g. c(x=3,y=0.4)
#'
#' @return ggplot
#' @export
#' @aliases plot_KOs_distribution
#' @examples
#' data("reporter_score_res")
#' plot_features_distribution(reporter_score_res, map_id = c("map05230", "map03010"))
plot_features_distribution <- function(reporter_res, map_id, text_size = 4, text_position = NULL, rug_length = 0.04) {
    ID=value=Z_score=Significantly=NULL
    stopifnot(inherits(reporter_res, "reporter_score"))

    ko_stat <- reporter_res$ko_stat

    A <- lapply(map_id, \(i)data.frame(ID = i, get_KOs(map_id = i, ko_stat = ko_stat, modulelist = reporter_res$modulelist))) %>% do.call(rbind, .)
    A$ID <- factor(A$ID, levels = map_id)
    if (nrow(A) < 1) {
        return(NULL)
    }

    df <- data.frame(value = ko_stat$Z_score)
    reddf <- reporter_res$reporter_s %>% dplyr::filter(ID %in% map_id)
    reddf$ID <- factor(reddf$ID, levels = map_id)

    if (is.null(text_position)) text_position <- c(x = max(df$value), y = 0.4)

    p <- ggplot() +
        geom_density(data = df, aes(x = value), fill = "lightblue", color = "black") +
        geom_vline(xintercept = median(df$value), linetype = 2) +
        geom_segment(data = A, aes(x = Z_score, xend = Z_score, y = 0, yend = -rug_length, color = Significantly)) +
        facet_grid(~ID) +
        labs(title = paste0("Z-score distribution of all features"), x = "Feature Z-score", y = "Density") +
        ggrepel::geom_text_repel(
            data = reddf,
            aes(x = text_position[1], y = text_position[2], label = paste0(
                "ReporterScore: ",
                round(ReporterScore, 2),
                "\np-value: ", signif(p.value, 4),
                "\np.adjust: ", signif(p.adjust, 4)
            )),
            color = "black", size = text_size
        ) +
        theme_classic() +
        scale_color_manual(name = "", values = c("Depleted" = "seagreen", "Enriched" = "orange", "None" = "grey", "Significant" = "red2"))
    p
}

#' Plot features trend in one pathway or module
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param box_color box and point color, default: c("#e31a1c","#1f78b4")
#' @param line_color line color, default: c("Depleted"="seagreen","Enriched"="orange","None"="grey")
#' @param show_number show the numbers.
#' @param scale scale the data by row.
#' @param feature_type show in the title ,default: KOs
#'
#' @import ggplot2
#' @return ggplot
#' @export
#' @aliases plot_KOs_in_pathway
#' @examples
#' data("reporter_score_res")
#' plot_features_in_pathway(ko_stat = reporter_score_res, map_id = "map00860")
plot_features_in_pathway <- function(ko_stat, map_id = "map00780",
                                     modulelist = NULL, select_ko = NULL,
                                     box_color = reporter_color, show_number = TRUE, scale = FALSE, feature_type = "KOs",
                                     line_color = c("Depleted" = "seagreen", "Enriched" = "orange", "None" = "grey", "Significant" = "red2")) {
    Group <- value <- KO_id <- Group2 <- value2 <- Group1 <- value1 <- type <- Significantly <- KOlist <- p.adjust <- n <- NULL
    if (is.null(names(line_color))) names(line_color) <- c("Depleted", "Enriched", "None", "Significant")[seq_along(line_color)]
    pcutils::lib_ps("ggnewscale", "reshape2", library = FALSE)
    flag <- FALSE
    if (inherits(ko_stat, "reporter_score")) {
        reporter_res <- ko_stat
        kodf <- reporter_res$kodf
        ko_stat <- reporter_res$ko_stat
        group <- reporter_res$group
        metadata <- reporter_res$metadata
        modulelist <- reporter_res$modulelist
        if (is.character(modulelist)) {
            load_GOlist(envir = environment())
            modulelist <- eval(parse(text = modulelist))
        }
        flag <- TRUE
        RS <- reporter_res$reporter_s[reporter_res$reporter_s$ID == map_id, "ReporterScore"]
    }

    if (is.null(select_ko)) {
        A <- get_KOs(map_id = map_id, ko_stat = ko_stat, modulelist = modulelist)
    } else {
        A <- ko_stat[ko_stat$KO_id %in% select_ko, ]
    }

    if (nrow(A) < 1) {
        return(NULL)
    }

    if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) stop("check your modulelist format!")
    Description <- modulelist[modulelist$id == map_id, "Description"]

    vs_group <- attributes(ko_stat)$vs_group
    if (identical(vs_group, "Numeric variable")) {
        if (!flag) stop("please input a reporter_score object when vs_group is a numeric variable")
        kodf_A <- kodf[rownames(A), ]
        if (scale) kodf_A <- pcutils::trans(kodf_A, method = "standardize", margin = 1)
        kodf_A$KO_id <- rownames(kodf_A)
        line_df <- reshape2::melt(kodf_A, id.vars = "KO_id", variable.name = "Sample_id")
        line_df <- dplyr::left_join(line_df,
                                    data.frame(Sample_id = rownames(metadata),
                                               metadata[group], check.names = FALSE),
                                    by = "Sample_id")
        line_df <- dplyr::left_join(line_df, A[, c("KO_id", "Significantly")], by = "KO_id")
        if (show_number) {
            num <- dplyr::count(A, Significantly) %>% dplyr::mutate(label = paste0(Significantly, ": ", n))
            line_df$Significantly <- setNames(num$label, num$Significantly)[line_df$Significantly]
            # line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
            names(line_color) <- setNames(num$label, num$Significantly)[names(line_color)]
        }
        p <- ggplot(data = line_df, aes(x = get(group), y = value, color = Significantly)) +
            geom_point(show.legend = FALSE, alpha = 0.5) +
            # geom_line(aes(group=KO_id))+
            geom_smooth(aes(group = KO_id), se = FALSE) +
            labs(
                title = ifelse(is.null(select_ko),
                    ifelse(map_id == Description,
                        paste0(feature_type, " in ", map_id),
                        paste0(feature_type, " in ", map_id, " (", Description, ")")
                    ), "Selected KOs"
                ),
                x = group, y = "Abundance"
            ) +
            scale_color_manual(values = line_color) +
            theme_classic(base_size = 13) +
            theme(axis.text = element_text(color = "black"))
        if (flag) p <- p + labs(subtitle = paste0("ReporterScore: ", round(RS, 3)))
        return(p)
    }

    colnames(A)[colnames(A) %in% paste0("average_", vs_group)] <- vs_group
    if (scale) A[, c(vs_group)] <- pcutils::trans(A[, c(vs_group)], method = "standardize", margin = 1)

    box_df <- reshape2::melt(A[, c("KO_id", vs_group)], id.vars = "KO_id", variable.name = "Group")
    box_df$Group <- factor(box_df$Group, levels = (vs_group))
    line_df <- data.frame()

    for (i in seq_len(length(vs_group) - 1)) {
        tmp <- A[, c("KO_id", vs_group[i:(i + 1)], "Significantly")]
        colnames(tmp) <- c("KO_id", "value1", "value2", "Significantly")
        tmp$Group1 <- vs_group[i]
        tmp$Group2 <- vs_group[i + 1]
        line_df <- rbind(line_df, tmp)
    }
    if (show_number) {
        num <- dplyr::count(A, Significantly) %>% dplyr::mutate(label = paste0(Significantly, ": ", n))
        line_df$Significantly <- setNames(num$label, num$Significantly)[line_df$Significantly]
        # line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
        names(line_color) <- setNames(num$label, num$Significantly)[names(line_color)]
    }

    p <- ggplot() +
        geom_boxplot(data = box_df, aes(x = Group, y = value, color = Group), show.legend = FALSE) +
        geom_point(data = box_df, aes(x = Group, y = value, color = Group), show.legend = FALSE) +
        scale_color_manual(values = pcutils::get_cols(nlevels(box_df$Group), box_color)) +
        labs(
            title = ifelse(is.null(select_ko),
                ifelse(map_id == Description,
                    paste0(feature_type, " in ", map_id),
                    paste0(feature_type, " in ", map_id, " (", Description, ")")
                ), "Selected KOs"
            ),
            x = "", y = "Abundance"
        ) +
        ggnewscale::new_scale_color() +
        geom_segment(
            data = line_df,
            aes(x = Group2, y = value2, xend = Group1, yend = value1, color = Significantly)
        ) +
        scale_color_manual(values = line_color) +
        theme_classic(base_size = 13) +
        theme(axis.text = element_text(color = "black"))
    if (flag) p <- p + labs(subtitle = paste0("ReporterScore: ", round(RS, 3)))
    p
}

#' Plot features boxplot
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples. or result of `get_reporter_score`
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param box_param parameters pass to \code{\link[pcutils]{group_box}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param KO_description show KO description rather than KO id.
#' @param str_width str_width to wrap
#' @param only_sig only show the significant features
#'
#' @aliases plot_KOs_box
#' @export
#' @return ggplot
#' @examples
#' data("reporter_score_res")
#' plot_features_box(reporter_score_res,
#'     select_ko = c("K00059", "K00208", "K00647", "K00652", "K00833", "K01012")
#' )
#' plot_features_box(reporter_score_res, select_ko = "K00059", KO_description = TRUE)
plot_features_box <- function(kodf, group = NULL, metadata = NULL,
                              map_id = "map00780", select_ko = NULL, only_sig = FALSE,
                              box_param = NULL,
                              modulelist = NULL,
                              KO_description = FALSE, str_width = 50) {
    Significantly=NULL
    flag <- FALSE
    if (inherits(kodf, "reporter_score")) {
        reporter_res <- kodf
        kodf <- reporter_res$kodf
        group <- reporter_res$group
        metadata <- reporter_res$metadata
        modulelist <- reporter_res$modulelist
        if (is.character(modulelist)) {
            load_GOlist(envir = environment())
            modulelist <- eval(parse(text = modulelist))
        }
        flag <- TRUE
    }

    if (is.null(select_ko)) select_ko <- get_KOs(map_id = map_id, modulelist = modulelist)
    if (only_sig & flag) {
        sig_names <- reporter_res$ko_stat %>%
            dplyr::filter(Significantly != "None") %>%
            rownames()
        select_ko <- intersect(select_ko, sig_names)
    }
    tkodf <- kodf %>%
        t() %>%
        as.data.frame()
    cols <- which(colnames(tkodf) %in% select_ko)

    if (length(cols) == 0) stop("No select KOs! check map_id or select_ko")
    if (length(cols) > 36) {
        message(("Too many KOs, do you still want to plot?"))
        flag <- readline("yes/no(y/n)?")
        if (!tolower(flag) %in% c("yes", "y")) {
            return(NULL)
        }
    }

    plotdat <- tkodf[, cols, drop = FALSE]

    if (KO_description) {
        load_KO_desc(envir = environment())
        if (grepl("C\\d{5}", colnames(plotdat)[1])) {
            load_Compound_htable(envir = environment())
            ko_desc <- Compound_htable[, c("Compound_id", "Compound_name")]
            colnames(ko_desc) <- c("KO_id", "KO_name")
        }
        newname <- ko_desc[match(colnames(plotdat), ko_desc$KO_id), "KO_name", drop = TRUE]
        if (all(is.na(newname))) warning("No description for KO found, are you sure rownames of kodf are KOs?")
        colnames(plotdat) <- ifelse(is.na(newname), colnames(plotdat), newname) %>% stringr::str_wrap(., width = str_width)
    }

    if (is.numeric(metadata[, group])) {
        p <- do.call(
            pcutils::my_lm,
            append(list(tab = plotdat, var = group, metadata = metadata), box_param)
        ) +
            theme_classic(base_size = 13) + theme(axis.text = element_text(color = "black")) +
            scale_fill_manual(values = pcutils::get_cols(nlevels(metadata[, group]), reporter_color)) +
            scale_color_manual(values = pcutils::get_cols(nlevels(metadata[, group]), reporter_color))
        return(p)
    }

    metadata[, group] <- factor(metadata[, group], levels = levels(factor(metadata[, group])))
    do.call(
        pcutils::group_box,
        append(
            list(tab = plotdat, group = group, metadata = metadata),
            pcutils::update_param(list(p_value1 = TRUE, trend_line = TRUE), box_param)
        )
    ) +
        theme_classic(base_size = 13) + theme(axis.text = element_text(color = "black")) +
        scale_fill_manual(values = pcutils::get_cols(nlevels(metadata[, group]), reporter_color)) +
        scale_color_manual(values = pcutils::get_cols(nlevels(metadata[, group]), reporter_color))
}

#' Plot features heatmap
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples. or result of `get_reporter_score`
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param heatmap_param parameters pass to \code{\link[pheatmap]{pheatmap}}
#' @param KO_description show KO description rather than KO id.
#' @param str_width str_width to wrap
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param only_sig only show the significant KOs
#' @param columns change columns
#'
#' @aliases plot_KOs_heatmap
#' @return ggplot
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_features_heatmap(reporter_score_res, map_id = "map00780")
plot_features_heatmap <- function(kodf, group = NULL, metadata = NULL,
                                  map_id = "map00780", select_ko = NULL,
                                  only_sig = FALSE, columns = NULL,
                                  modulelist = NULL,
                                  KO_description = FALSE, str_width = 50,
                                  heatmap_param = list()) {
    Significantly=NULL
    pcutils::lib_ps("pheatmap", library = FALSE)
    flag <- FALSE
    if (inherits(kodf, "reporter_score")) {
        reporter_res <- kodf
        kodf <- reporter_res$kodf
        group <- reporter_res$group
        metadata <- reporter_res$metadata
        modulelist <- reporter_res$modulelist
        if (is.character(modulelist)) {
            load_GOlist(envir = environment())
            modulelist <- eval(parse(text = modulelist))
        }
        flag <- TRUE
    }

    if (is.null(select_ko)) select_ko <- get_KOs(map_id = map_id, modulelist = modulelist)

    if (only_sig & flag) {
        sig_names <- reporter_res$ko_stat %>%
            dplyr::filter(Significantly != "None") %>%
            rownames()
        select_ko <- intersect(select_ko, sig_names)
    }

    cols <- which(rownames(kodf) %in% select_ko)

    if (length(cols) == 0) stop("No select KOs! check map_id or select_ko")
    plotdat <- kodf[cols, , drop = FALSE]
    if (!is.null(columns)) plotdat <- plotdat[, columns]

    if (KO_description) {
        load_KO_desc(envir = environment())
        if (grepl("C\\d{5}", rownames(plotdat)[1])) {
            load_Compound_htable(envir = environment())
            ko_desc <- Compound_htable[, c("Compound_id", "Compound_name")]
            colnames(ko_desc) <- c("KO_id", "KO_name")
        }
        newname <- ko_desc[match(rownames(plotdat), ko_desc$KO_id), "KO_name", drop = TRUE]
        if (all(is.na(newname))) warning("No description for KO found, are you sure rownames of kodf are KOs?")
        rownames(plotdat) <- ifelse(is.na(newname), rownames(plotdat), newname) %>% stringr::str_wrap(., width = str_width)
    }

    if (is.numeric(metadata[, group])) {
        annotation_colors <- NULL
    } else {
        metadata[, group] <- factor(metadata[, group], levels = levels(factor(metadata[, group])))
        annotation_colors <- list(pcutils::get_cols(nlevels(factor(metadata[, group])), reporter_color))
        names(annotation_colors) <- group
        names(annotation_colors[[group]]) <- levels(factor(metadata[, group]))
    }

    do.call(
        pheatmap::pheatmap,
        pcutils::update_param(
            list(
                mat = plotdat, cluster_cols = FALSE,
                color = pcutils::get_cols(n = 100, pal = "bluered"),
                annotation_colors = annotation_colors,
                annotation_col = metadata[group], scale = "row"
            ),
            heatmap_param
        )
    )
}


#' Plot c_means result
#'
#' @param rsa_cm_res a cm_res object
#' @param ... additional
#' @param filter_membership filter membership 0~1.
#' @param mode 1~2
#' @param show.clust.cent show cluster center?
#' @param show_num show number of each cluster?
#'
#' @rdname RSA_by_cm
#' @export
#' @return ggplot
plot_c_means <- function(rsa_cm_res, filter_membership, mode = 1, show.clust.cent = TRUE, show_num = TRUE, ...) {
    stopifnot(inherits(rsa_cm_res, "rs_by_cm"))
    plot(rsa_cm_res$cm_res, filter_membership = filter_membership, show.clust.cent = show.clust.cent, mode = mode, show_num = show_num, ...)
}

#' Plot c_means result
#'
#' @param x a cm_res object
#' @param ... additional
#' @param filter_membership filter membership
#' @param mode 1~2
#' @param show.clust.cent show cluster center?
#' @param show_num show number of each cluster?
#'
#' @return ggplot
#' @exportS3Method
#' @method plot cm_res
plot.cm_res <- function(x, filter_membership, mode = 1, show.clust.cent = TRUE, show_num = TRUE, ...) {
    Group= value= Name= Cluster = Membership=NULL
    lib_ps("factoextra", library = FALSE)

    cm_data <- x$cm_data
    pcutils::dabiao("filter clusters, Membership >= ", filter_membership)
    cm_group <- cm_data[cm_data$Membership >= filter_membership, ] # 筛选部分显著被聚类的项
    if (show_num) {
        tmp <- cm_group %>%
            count(Cluster) %>%
            mutate(new_cluster = paste0(Cluster, ": ", n))
        cm_group$Cluster <- setNames(tmp$new_cluster, tmp$Cluster)[cm_group$Cluster]
    }
    data_scaled <- x$data_scaled

    if (mode == 2) {
        # show the cluster
        p <- factoextra::fviz_cluster(list(data = data_scaled[rownames(cm_group), ], cluster = cm_group$Cluster),
                                      geom = c("point"),
                                      ellipse = TRUE,
                                      ellipse.alpha = 0.3, # used to be 0.6 if only points are plotted.
                                      ellipse.type = "norm",
                                      ellipse.level = 0.68,
                                      repel = TRUE, show.clust.cent = show.clust.cent
        ) + reporter_theme
    }
    if (mode == 1) {
        cm_group.melt <- reshape2::melt(cm_group, id.vars = c("Cluster", "Membership", "Name", "Weight"), variable.name = "Group")
        cm_group.melt$Cluster <- factor(cm_group.melt$Cluster)

        p <- ggplot() +
            geom_line(
                data = cm_group.melt,
                aes(x = Group, y = value, group = Name, color = Cluster, alpha = Membership),
                linewidth = 0.8
            ) +
            reporter_theme +
            scale_x_discrete(expand = c(0, 0)) +
            theme(plot.margin = unit(c(1, 2, 1, 1), "lines"))
        if (show.clust.cent) {
            centers <- x$centers %>% as.data.frame()

            if (show_num) {
                centers$Cluster <- tmp$new_cluster
            } else {
                centers$Cluster <- rownames(centers)
            }

            centers_dat <- reshape2::melt(centers, id.vars = "Cluster", variable.name = "Group")

            p <- p +
                scale_alpha_continuous(range = c(0.2, 0.5)) +
                geom_line(data = centers_dat, aes(x = Group, y = value, group = Cluster, color = Cluster), linewidth = 3)
        }
    }
    cols1 <- pcutils::get_cols(length(unique(cm_group$Cluster)))
    return(p + scale_color_manual(values = cols1) + scale_fill_manual(values = cols1))
}

#' Plot htable levels
#'
#' @param select select ids
#' @param type "ko", "module", "pathway", "compound"
#' @param htable custom a htable
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data("KO_abundance_test")
#' plot_htable(select = rownames(KO_abundance))
plot_htable <- function(type = "ko", select = NULL, htable = NULL) {
    n=level2_name=level1_name=NULL

    title <- ""
    if (is.null(htable)) {
        type <- match.arg(type, c("ko", "module", "pathway", "compound"))
        switch(type,
            "ko" = {
                prefix <- "KO_htable"
                title <- "KEGG KO annotation"
            },
            "module" = {
                prefix <- "Module_htable"
                title <- "KEGG module annotation"
            },
            "pathway" = {
                prefix <- "Pathway_htable"
                title <- "KEGG pathway annotation"
            },
            "compound" = {
                prefix <- "Compound_htable"
                title <- "Compounds with biological roles"
            }
        )
        load_htable(type, envir = environment())
        htable <- get(prefix, envir = environment())
        switch(type,
            "ko" = {
                htable <- htable %>% select(level1_name, level2_name, KO_id)
            },
            "module" = {
                htable <- htable %>% select(module2_name, module3_name, Module_id)
            },
            "pathway" = {
                htable <- htable %>% select(level1_name, level2_name, Pathway_id)
            },
            "compound" = {
                htable <- htable %>%
                    filter(Class == "Compounds with biological roles") %>%
                    select(compound1_name, compound2_name, Compound_id)
            }
        )
    }
    if (ncol(htable) != 3) stop("htable should be a three-columns table.")
    colnames(htable) <- c("level1_name", "level2_name", "KO_id")

    if (!is.null(select)) htable <- dplyr::filter(htable, KO_id %in% select)

    if (nrow(htable) < 1) {
        return(NULL)
    }

    dplyr::count(htable, level1_name, level2_name) -> a
    a$level2_name <- factor(a$level2_name, levels = a$level2_name)
    patt <- setNames(pcutils::get_cols(length(unique(a$level1_name))), unique(a$level1_name))

    ggplot(a) +
        geom_col(aes(x = n, y = level2_name, fill = level1_name)) +
        scale_fill_manual(values = patt) +
        geom_text(aes(x = n, y = level2_name, label = n),
            nudge_x = max(a$n) * 1 / 200, hjust = 0
        ) +
        labs(y = NULL, x = "Items number", title = title) +
        theme_classic() +
        theme(axis.text.y = element_text(color = patt[a$level1_name])) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, max(a$n) * 1.1))
}


#' plot_KEGG_map
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param map_id the pathway or module id
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param feature one of "ko", "gene", "compound"
#' @param type "pathway" or "module" for default KOlist for microbiome, "CC", "MF", "BP", "ALL" for default GOlist for homo sapiens. And org in listed in 'https://www.genome.jp/kegg/catalog/org_list.html' such as "hsa" (if your kodf is come from a specific organism, you should specify type here).
#' @param color_var use which variable to color
#' @param save_dir where to save the png files
#' @param color color
#'
#' @references https://zhuanlan.zhihu.com/p/357687076
#' @return png files
#' @export
#'
#' @examples
#' message("The following example will download some files:")
#' \dontrun{
#' data("reporter_score_res")
#' plot_KEGG_map(reporter_score_res$ko_stat, map_id = "map00780",
#'      type = "pathway", feature = "ko", color_var = "Z_score")
#' }
plot_KEGG_map <- function(ko_stat, map_id = "map00780", modulelist = NULL, type = "pathway", feature = "ko",
                          color_var = "Z_score", save_dir = "ReporterScore_temp_download/", color = c("blue", "grey", "red")) {
    if (inherits(ko_stat, "reporter_score")) {
        reporter_res <- ko_stat
        ko_stat <- reporter_res$ko_stat
        modulelist <- reporter_res$modulelist
        if (is.character(modulelist)) {
            load_GOlist(envir = environment())
            modulelist <- eval(parse(text = modulelist))
        }
        flag <- TRUE
        RS <- reporter_res$reporter_s[reporter_res$reporter_s$ID == map_id, "ReporterScore"]
    }
    if (is.null(modulelist)) modulelist <- get_modulelist(type = type, feature = feature, verbose = FALSE)

    A <- get_KOs(map_id = map_id, ko_stat = ko_stat, modulelist = modulelist)
    if (nrow(A) < 1) A <- ko_stat

    lib_ps("pathview")
    if (type == "pathway") type <- "ko"
    if (type == "module") stop("need pathway")
    pathway.id <- gsub("[^0-9]", "", map_id)

    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- paste0(type, pathway.id, ".", color_var, ".png")

    if (is.numeric(A[, color_var, drop = TRUE])) {
        discrete <- FALSE
        limit <- max(abs(A[, color_var, drop = TRUE]))
        both.dirs <- TRUE
    } else {
        stop("need numeric")
    }
    params <- list(
        pathway.id = pathway.id, species = type, kegg.dir = save_dir, out.suffix = color_var, res = 500,
        discrete = discrete, limit = limit, both.dirs = both.dirs, low = color[1], mid = color[2], high = color[3]
    )

    if (feature %in% c("ko")) {
        do.call(pathview::pathview, pcutils::update_param(list(gene.data = A[color_var], gene.idtype = "KEGG"), params))
    }
    if (feature %in% "gene") {
        message('please make sure the rownames of your input is entrezid\ndo a transfer like\n`a=clusterProfiler::bitr(rownames(ko_stat),"SYMBOL","ENTREZID",OrgDb = org.Hs.eg.db::org.Hs.eg.db)\nrownames(ko_stat)=a$ENTREZID`')

        do.call(pathview::pathview, pcutils::update_param(list(gene.data = A[color_var], gene.idtype = "entrez"), params))
    }
    if (feature == "compound") {
        do.call(pathview::pathview, pcutils::update_param(list(cpd.data = A[color_var], cpd.idtype = "kegg"), params))
    }
    file.copy(filename, file.path(save_dir, filename), overwrite = TRUE)
    file.remove(filename)
    message("result have been saved in ", file.path(save_dir, filename))
    # pcutils::read.file(file.path(save_dir, filename))
}
