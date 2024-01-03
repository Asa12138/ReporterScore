#' Update files from KEGG
#' @export
#' @return No value
#' @param download_dir where to save the .keg file?
#'
#' @examples
#' message("The following example will download some files:")
#' \donttest{
#' # update_KEGG()
#' }
update_KEGG <- function(download_dir = NULL) {
    pcutils::dabiao("1.`update_KO_file`")
    update_KO_file(download_dir = download_dir)
    pcutils::dabiao("2.`update_KO_htable`")
    update_htable(type = "ko", download_dir = download_dir)
    pcutils::dabiao("3.`update_Pathway_htable`")
    update_htable(type = "pathway", download_dir = download_dir)
    pcutils::dabiao("3.`update_Module_htable`")
    update_htable(type = "module", download_dir = download_dir)
    pcutils::dabiao("3.`update_Compound_htable`")
    update_htable(type = "compound", download_dir = download_dir)
}

#' Update the KO2Pathway and CPD2Pathway files from KEGG
#'
#' @param download_dir the ReporterScore user dir location, detect automatically by `tools::R_user_dir("ReporterScore")`.
#' @param RDSfile saved KO_files.RDS file
#'
#' @export
#' @description
#' Download links:
#'
#' \code{https://rest.kegg.jp/list/pathway}
#' \code{https://rest.kegg.jp/link/pathway/ko}
#' \code{https://rest.kegg.jp/link/pathway/compound}
#' \code{https://rest.kegg.jp/list/module}
#' \code{https://rest.kegg.jp/link/module/ko}
#' \code{https://rest.kegg.jp/link/module/compound}
#' @rdname update_KEGG
update_KO_file <- function(download_dir = NULL, RDSfile = NULL) {
    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    if (is.null(download_dir)) {
        dd <- "ReporterScore_temp_download"
    } else {
        dd <- download_dir
    }
    if (!dir.exists(dd)) dir.create(dd)

    if (is.null(RDSfile)) {
        pcutils::dabiao("Trying to download files from https://rest.kegg.jp/ .")

        pathway_desc <- KEGGREST::keggList("pathway")
        pathway_desc <- data.frame(id = names(pathway_desc), pathway = pathway_desc)
        module_desc <- KEGGREST::keggList("module")
        module_desc <- data.frame(id = names(module_desc), pathway = module_desc)

        path2ko <- KEGGREST::keggLink("pathway", "ko")
        path2ko <- data.frame(pathway = (path2ko), ko = names(path2ko), row.names = NULL)
        module2ko <- KEGGREST::keggLink("module", "ko")
        module2ko <- data.frame(pathway = (module2ko), ko = names(module2ko), row.names = NULL)
        path2compound <- KEGGREST::keggLink("pathway", "compound")
        path2compound <- data.frame(pathway = (path2compound), ko = names(path2compound), row.names = NULL)
        module2compound <- KEGGREST::keggLink("module", "compound")
        module2compound <- data.frame(pathway = (module2compound), ko = names(module2compound), row.names = NULL)

        KO_files <- list(
            pathway_desc = pathway_desc, module_desc = module_desc,
            path2ko = path2ko, module2ko = module2ko, path2compound = path2compound, module2compound = module2compound
        )
        saveRDS(KO_files, file = file.path(dd, paste0("KO_files", Sys.time(), ".RDS")))
    } else {
        KO_files <- readRDS(RDSfile)
    }

    KOlist <- list()
    KOlist$pathway <- KO_files$path2ko %>%
        dplyr::mutate(pathway = gsub("path:", "", pathway), ko = gsub("ko:", "", ko)) %>%
        dplyr::filter(grepl("map", pathway)) %>%
        custom_modulelist(., KO_files$pathway_desc)
    KOlist$module <- KO_files$module2ko %>%
        dplyr::mutate(pathway = gsub("md:", "", pathway), ko = gsub("ko:", "", ko)) %>%
        custom_modulelist(., KO_files$module_desc)
    CPDlist <- list()
    CPDlist$pathway <- KO_files$path2compound %>%
        dplyr::mutate(pathway = gsub("path:", "", pathway), ko = gsub("cpd:", "", ko)) %>%
        dplyr::filter(grepl("map", pathway)) %>%
        custom_modulelist(., KO_files$pathway_desc)
    CPDlist$module <- KO_files$module2ko %>%
        dplyr::mutate(pathway = gsub("md:", "", pathway), ko = gsub("cpd:", "", ko)) %>%
        custom_modulelist(., KO_files$module_desc)
    attributes(CPDlist)$download_time <- attributes(KOlist)$download_time <- Sys.time()

    save(KOlist, file = paste0(pack_dir, "/new_KOlist.rda"))
    save(CPDlist, file = paste0(pack_dir, "/new_CPDlist.rda"))
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
}


#' Update some htable from KEGG
#'
#' @param type "ko", "module", "pathway", "compound" ...
#' @param keg_file path of a .keg file, such as ko00001.keg from https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext.
#' @param download save the .keg file?
#' @param download_dir where to save the .keg file?
#'
#' @export
#' @rdname update_KEGG
#'
update_htable <- function(type, keg_file = NULL, download = FALSE, download_dir = NULL) {
    type <- match.arg(type, c("ko", "module", "pathway", "compound"))

    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    pcutils::lib_ps("KEGGREST", library = FALSE)
    if (download) {
        if (is.null(download_dir)) {
            dd <- "ReporterScore_temp_download"
        } else {
            dd <- download_dir
        }
    }

    if (type == "ko") {
        br_id <- "ko00001"
        if (!is.null(keg_file)) {
            if (!grepl(br_id, keg_file)) stop("keg_file should be path of ", br_id, ".keg.")
            lines <- readLines(keg_file, warn = FALSE)
        } else {
            KEGGREST::keggGet(paste0("br:", br_id)) -> a
            if (download) writeLines(a, con = file.path(dd, paste0(br_id, "_", Sys.time(), ".keg")), sep = "")
            lines <- strsplit(a, "\n")[[1]]
        }
        df_res <- brite2df(lines)

        Cdf <- sub(" ", "\t", df_res[["C"]]) %>% strsplit2(., "\t", colnames = c("path_id", "path_name"))
        KO_htable <- data.frame(
            sub("\\d* ", "", df_res$A),
            sub("\\d* ", "", df_res$B),
            Cdf,
            df_res[, c("id", "name")]
        )

        colnames(KO_htable) <- c(
            "level1_name", "level2_name",
            "level3_id", "level3_name", "KO_id", "KO_name"
        )

        KO_htable$level3_id <- ifelse(KO_htable$level1_name %in% c("Brite Hierarchies", "Not Included in Pathway or Brite"),
            KO_htable$level3_id, paste0("map", KO_htable$level3_id)
        )

        KO_htable$level3_name <- sub(" \\[.*:ko.*\\]", "", KO_htable$level3_name)
        attributes(KO_htable)$download_time <- attributes(KO_htable)$build_time <- Sys.time()
        save(KO_htable, file = paste0(pack_dir, "/new_KO_htable.rda"))
    }
    if (type == "module") {
        br_id <- "ko00002"
        if (!is.null(keg_file)) {
            if (!grepl(br_id, keg_file)) stop("keg_file should be path of ", br_id, ".keg.")
            lines <- readLines(keg_file, warn = FALSE)
        } else {
            KEGGREST::keggGet(paste0("br:", br_id)) -> a
            if (download) writeLines(a, con = file.path(dd, paste0(br_id, "_", Sys.time(), ".keg")), sep = "")
            lines <- strsplit(a, "\n")[[1]]
        }
        df_res <- brite2df(lines)

        df_res$A <- gsub("</?b>", "", df_res$A)
        df_res$B <- gsub("</?b>", "", df_res$B)
        Module_htable <- df_res
        colnames(Module_htable) <- c(
            "module1_name", "module2_name",
            "module3_name", "Module_id", "Module_name"
        )

        attributes(Module_htable)$download_time <- attributes(Module_htable)$build_time <- Sys.time()
        save(Module_htable, file = paste0(pack_dir, "/new_Module_htable.rda"))
    }
    if (type == "pathway") {
        br_id <- "br08901"
        if (!is.null(keg_file)) {
            if (!grepl(br_id, keg_file)) stop("keg_file should be path of ", br_id, ".keg.")
            lines <- readLines(keg_file, warn = FALSE)
        } else {
            KEGGREST::keggGet(paste0("br:", br_id)) -> a
            if (download) writeLines(a, con = file.path(dd, paste0(br_id, "_", Sys.time(), ".keg")), sep = "")
            lines <- strsplit(a, "\n")[[1]]
        }
        df_res <- brite2df(lines)
        Pathway_htable <- df_res
        colnames(Pathway_htable) <- c("level1_name", "level2_name", "Pathway_id", "Pathway_name")
        Pathway_htable$Pathway_id <- paste0("map", Pathway_htable$Pathway_id)
        attributes(Pathway_htable)$download_time <- attributes(Pathway_htable)$build_time <- Sys.time()
        save(Pathway_htable, file = paste0(pack_dir, "/new_Pathway_htable.rda"))
    }
    if (type == "compound") {
        compounds_brites <- brite_file2df(br_id = "br08902") %>% dplyr::filter(B == "Compounds")
        if (!is.null(keg_file)) {
            if (!dir.exists(keg_file)) {
                stop(
                    "when type='compound', keg_file should be path of a directory contains: ",
                    paste0(compounds_brites$id, ".keg", collapse = ", "),
                    " , as there are more than .keg files needed."
                )
            }
            keg_files <- list.files("ReporterScore_temp_download/",
                full.names = TRUE,
                pattern = paste0(compounds_brites$id, collapse = "|")
            ) %>%
                grep("*.keg", x = ., value = TRUE)
        }

        compounds_list <- lapply(compounds_brites$id, \(i){
            if (!is.null(keg_file) & any(grepl(i, keg_files))) {
                keg_file_tmp <- keg_files[grepl(i, keg_files)][1]
                message("use", keg_file_tmp, "\n")
                lines <- readLines(keg_file_tmp, warn = FALSE)
            } else {
                KEGGREST::keggGet(paste0("br:", i)) -> a
                if (download) writeLines(a, con = file.path(dd, paste0(i, "_", Sys.time(), ".keg")), sep = "")
                lines <- strsplit(a, "\n")[[1]]
            }
            brite2df2(lines)
        })
        df_res <- lapply(seq_along(compounds_list), \(i)compounds_list[[i]] <- cbind("Class" = compounds_brites$name[i], compounds_list[[i]]))
        df_res <- do.call(plyr::rbind.fill, df_res)

        Compound_htable <- df_res[, c("Class", "A", "B", "C", "D", "id", "name")]
        Compound_htable <- dplyr::mutate_all(Compound_htable, \(i)gsub(" \\[Fig\\]", "", i))

        colnames(Compound_htable) <- c(
            "Class", "compound1_name", "compound2_name",
            "compound3_name", "compound4_name", "Compound_id", "Compound_name"
        )

        attributes(Compound_htable)$download_time <- attributes(Compound_htable)$build_time <- Sys.time()
        save(Compound_htable, file = paste0(pack_dir, "/new_Compound_htable.rda"))
    }
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
}

#' Load the specific table (from KEGG)
#'
#' @param type "ko", "module", "pathway", "compound" ...
#' @param envir `.GlobalEnv` (default) or `environment()`
#' @param verbose logical
#'
#' @export
#' @return KO_htable in `.GlobalEnv`
#' @examples
#' load_htable("pathway")
#' head(KO_htable)
load_htable <- function(type, envir = .GlobalEnv, verbose = TRUE) {
    type <- match.arg(type, c("ko", "module", "pathway", "compound"))
    switch(type,
        "ko" = {
            prefix <- "KO_htable"
        },
        "module" = {
            prefix <- "Module_htable"
        },
        "pathway" = {
            prefix <- "Pathway_htable"
        },
        "compound" = {
            prefix <- "Compound_htable"
        }
    )

    KO_htable_file <- file.path(tools::R_user_dir("ReporterScore"), paste0("new_", prefix, ".rda"))

    if (file.exists(KO_htable_file)) {
        load(KO_htable_file, envir = envir)
    } else {
        data(list = prefix, package = "ReporterScore", envir = envir)
    }
    if (verbose) {
        pcutils::dabiao("load ", prefix)
        if (!is.null(attributes(get(prefix, envir = envir))$"download_time")) {
            pcutils::dabiao(paste0(prefix, " download time: ", attributes(get(prefix, envir = envir))$"download_time"))
            message("If you want to update ", prefix, ", use `update_htable(type='", type, "')`")
        }
    }
}

#' Load the KOlist (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return KOlist in `.GlobalEnv`
load_KOlist <- function(envir = .GlobalEnv, verbose = TRUE) {
    if (TRUE) {
        # KOlist_file=system.file("data","new_KOlist.rda",package = "ReporterScore")
        KOlist_file <- file.path(tools::R_user_dir("ReporterScore"), "new_KOlist.rda")
    }
    if (file.exists(KOlist_file)) {
        load(KOlist_file, envir = envir)
    } else {
        data("KOlist", package = "ReporterScore", envir = envir)
    }
    if (verbose) {
        pcutils::dabiao("load KOlist")
        if (!is.null(attributes(get("KOlist", envir = envir))$"download_time")) {
            pcutils::dabiao(paste0("KOlist download time: ", attributes(get("KOlist", envir = envir))$"download_time"))
            message("If you want to update KOlist, use `update_KO_file()`")
        }
    }
}

#' Load the CPDlist (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return CPDlist in `.GlobalEnv`
load_CPDlist <- function(envir = .GlobalEnv, verbose = TRUE) {
    if (TRUE) {
        # CPDlist_file=system.file("data","new_CPDlist.rda",package = "ReporterScore")
        CPDlist_file <- file.path(tools::R_user_dir("ReporterScore"), "new_CPDlist.rda")
    }
    if (file.exists(CPDlist_file)) {
        load(CPDlist_file, envir = envir)
    } else {
        data("CPDlist", package = "ReporterScore", envir = envir)
    }
    if (verbose) {
        pcutils::dabiao("load CPDlist")
        if (!is.null(attributes(get("CPDlist", envir = envir))$"download_time")) {
            pcutils::dabiao(paste0("CPDlist download time: ", attributes(get("CPDlist", envir = envir))$"download_time"))
            message("If you want to update CPDlist, use `update_KO_file()`")
        }
    }
}


#' Load the KO description (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return KO description in `.GlobalEnv`
load_KO_desc <- function(envir = .GlobalEnv, verbose = TRUE) {
    load_KO_htable(envir = environment(), verbose = verbose)
    ko_desc <- dplyr::distinct(KO_htable[, c("KO_id", "KO_name")]) %>% dplyr::arrange(KO_id)
    assign("ko_desc", ko_desc, envir = envir)
}

#' Load the KO_htable (from KEGG)
#'
#' @rdname load_htable
#' @export
#' @return KO_htable in `.GlobalEnv`
load_KO_htable <- function(envir = .GlobalEnv, verbose = TRUE) {
    load_htable(type = "ko", envir = envir, verbose = verbose)
}

#' Load the Pathway_htable (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return Pathway_htable in `.GlobalEnv`
load_Pathway_htable <- function(envir = .GlobalEnv, verbose = TRUE) {
    load_htable(type = "pathway", envir = envir, verbose = verbose)
}

#' Load the Module_htable (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return Module_htable in `.GlobalEnv`
load_Module_htable <- function(envir = .GlobalEnv, verbose = TRUE) {
    load_htable(type = "module", envir = envir, verbose = verbose)
}

#' Load the Compound_htable (from KEGG)
#'
#' @rdname load_htable
#'
#' @export
#' @return Compound_htable in `.GlobalEnv`
load_Compound_htable <- function(envir = .GlobalEnv, verbose = TRUE) {
    load_htable(type = "compound", envir = envir, verbose = verbose)
}

brite2df <- function(lines) {
    # 当确定查询的tip都在同一等级时（比如D）用这个
    lines <- lines[!grepl("^#", lines)]

    if (any(lines == "!")) {
        contents <- lines[(which(lines == "!")[1] + 1):(which(lines == "!")[2] - 1)]
    } else {
        contents <- lines
    }

    all_levels <- substr(contents, 1, 1)
    levels <- unique(all_levels)
    tip <- levels[length(levels)]

    edges <- as.list(levels)
    names(edges) <- levels
    node <- names(edges)[-length(levels)]

    df <- vector("list", sum(all_levels == tip))
    tip_n <- 1
    # Process the input file
    for (line in contents) {
        for (i in node) {
            if (grepl(paste0("^", i), line)) edges[[i]] <- line
            next
        }
        if (grepl(paste0("^", tip), line)) {
            edges[[tip]] <- line
            df[[tip_n]] <- edges
            tip_n <- tip_n + 1
        }
    }
    df_res <- do.call(rbind, lapply(df, as.data.frame))
    # 遇到一个很坑的情况，一个文件中不止一个tip，有的tip同时分布在B，C，D...
    # 应该先按照A分割处理后合并

    for (i in levels) {
        df_res[[i]] <- gsub(paste0("^", i, "\\s*"), "", df_res[[i]])
    }
    tipdf <- sub("  ", "\t", df_res[[tip]]) %>% strsplit2(., "\t", colnames = c("id", "name"))
    cbind(df_res[, node], tipdf)
}
brite2df2 <- function(lines, id_pattern = "\\s+[CG]\\d{5}  ") {
    # 当确定查询的tip的pattern时用这个
    lines <- lines[!grepl("^#", lines)]
    if (any(lines == "!")) {
        contents <- lines[(which(lines == "!")[1] + 1):(which(lines == "!")[2] - 1)]
    } else {
        contents <- lines
    }

    all_levels <- substr(contents, 1, 1)
    levels <- unique(all_levels)

    edges <- as.list(levels)
    names(edges) <- levels

    df <- vector("list", sum(grepl(id_pattern, contents)))
    tip_n <- 1
    # Process the input file
    for (line in contents) {
        perfix <- substr(line, 1, 1)
        edges[[perfix]] <- line
        if (grepl(id_pattern, line)) {
            df[[tip_n]] <- edges[seq_len(which(levels == perfix))]
            tip_n <- tip_n + 1
        }
    }
    df <- lapply(df, as.data.frame)

    ncols <- vapply(df, ncol, numeric(1))
    df_res <- list()
    for (col in unique(ncols)) {
        tmp <- df[ncols == col]
        tmp <- do.call(rbind, tmp)
        for (i in levels[seq_len(col)]) {
            tmp[[i]] <- gsub(paste0("^", i, "\\s*"), "", tmp[[i]])
        }
        tipdf <- sub("  ", "\t", tmp[[levels[col]]]) %>% strsplit2(., "\t", colnames = c("id", "name"))
        df_res[[i]] <- cbind(tmp[, seq_len(col - 1)], tipdf)
    }

    # 遇到一个很坑的情况，一个文件中不止一个tip，有的tip同时分布在B，C，D...
    # 应该先按照A分割处理后合并
    df_res1 <- do.call(plyr::rbind.fill, df_res)
    df_res1[, c(levels[-length(levels)], "id", "name")]
}
brite_file2df <- function(keg_file = NULL, br_id = NULL, id_pattern = NULL) {
    some_pattern <- c(
        "br08001" = "\\s+[CG]\\d{5}"
    )
    if (!is.null(keg_file)) {
        lines <- readLines(keg_file, warn = FALSE)
    } else if (!is.null(br_id)) {
        pcutils::lib_ps("KEGGREST", library = FALSE)
        KEGGREST::keggGet(paste0("br:", br_id)) -> a
        lines <- strsplit(a, "\n")[[1]]
    }
    if (is.null(id_pattern)) {
        return(brite2df(lines))
    } else {
        return(brite2df2(lines, id_pattern = id_pattern))
    }
}
# brite_file2df(br_id = "br08902")


pre_gene_in_pathway <- function(pathway) {
    genes <- pathway[["GENE"]]
    if (is.null(genes)) {
        return(NULL)
    }
    kegg_gene_id <- genes[seq(1, length(genes), 2)]
    input_strings <- genes[seq(2, length(genes), 2)]

    data <- data.frame(
        gene_symbol = sub(";.*", "", input_strings),
        gene_desc = sub(".*;\\s*(.*?)\\s+\\[.*", "\\1", input_strings),
        KO_id = sub(".*\\[KO:(.*?)\\].*", "\\1", input_strings)
        # EC = sub(".*\\[EC:(.*?)\\].*", "\\1", input_strings)
    )
    data <- data.frame(pathway_id = pathway$ENTRY, kegg_gene_id = kegg_gene_id, data, row.names = NULL)
    data
}
pre_compound_in_pathway <- function(pathway) {
    compound <- pathway[["COMPOUND"]]
    if (is.null(compound)) {
        return(NULL)
    }

    data <- data.frame(pathway_id = pathway$ENTRY, kegg_compound_id = names(compound), kegg_compound_name = compound, row.names = NULL)
    data
}

#' Get pathway information for an organism from kegg
#'
#' @param org kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html, default, "hsa"
#' @param RDS_file path of a org.RDS file if you saved before.
#' @param download save the .keg file?
#' @param download_dir where to save the .keg file?
#'
#' @aliases download_org_pathway
#' @export
#' @rdname update_KEGG
update_org_pathway <- function(org = "hsa", RDS_file = NULL, download = TRUE, download_dir = NULL) {
    lib_ps("KEGGREST", library = FALSE)
    pcutils::dabiao("get pathway of '", org, "' from KEGG.")

    if (download) {
        if (is.null(download_dir)) {
            dd <- "ReporterScore_temp_download"
        } else {
            dd <- download_dir
        }
    }
    if (!is.null(RDS_file)) {
        if (!grepl(org, RDS_file)) stop("keg_file should be path of ", br_id, ".RDS.")
        hsa_genels <- readRDS(RDS_file)
        org_pathway <- hsa_genels$org_pathway
        hsa_gene <- hsa_genels$hsa_gene
    } else {
        org_pathway <- KEGGREST::keggList("pathway", org) # 获取KEGG数据库中所有人类通路
        # 暂时获取不了module，只有pathway
        hsa_path <- data.frame(org_pathway) # 转成数据框,方便后续分析
        hsa_path$pathID <- rownames(hsa_path) # 提取pathway ID
        org_pathway <- hsa_path
        message("total ", nrow(hsa_path), " pathways")
        pcutils::dabiao("get each pathway information from KEGG.")
        hsa_gene <- list()
        hsa_gene <- lapply(seq_len(nrow(hsa_path)), \(i) {
            if (i %% 50 == 0) message(i, " pathways done.")
            hsa_info <- KEGGREST::keggGet(hsa_path[i, "pathID"])
            hsa_info[[1]]
        })
        names(hsa_gene) <- hsa_path$pathID
        if (download) saveRDS(list(org_pathway = org_pathway, hsa_gene = hsa_gene), file = file.path(dd, paste0(org, "_", Sys.time(), ".RDS")))
    }

    all_hsa_gene <- lapply(hsa_gene, pre_gene_in_pathway)
    all_hsa_compound <- lapply(hsa_gene, pre_compound_in_pathway)

    all_org_gene <- do.call(rbind, all_hsa_gene)
    rownames(all_org_gene) <- NULL
    all_org_compound <- do.call(rbind, all_hsa_compound)
    rownames(all_org_compound) <- NULL

    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    a <- list(org_pathway = org_pathway, all_org_gene = all_org_gene, all_org_compound = all_org_compound)
    attr(a, "download_time") <- Sys.time()
    assign(paste0(org, "_kegg_pathway"), a, envir = environment())
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
    save(list = paste0(org, "_kegg_pathway"), file = file.path(pack_dir, paste0(org, "_kegg_pathway.rda")))
}

#' Load the pathway information for an organism (from KEGG)
#'
#' @param org kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html, default, "hsa"
#'
#' @export
#' @return KOlist in `.GlobalEnv`
#' @rdname load_htable
load_org_pathway <- function(org = "hsa", envir = .GlobalEnv, verbose = TRUE) {
    pack_dir <- tools::R_user_dir("ReporterScore")
    path_file <- file.path(pack_dir, paste0(org, "_kegg_pathway.rda"))
    if (file.exists(path_file)) {
        load(path_file, envir = envir)
    } else if (org %in% c("hsa", "mmu")) {
        data(paste0(org, "_kegg_pathway"), package = "ReporterScore", envir = envir)
    } else {
        message("No pathway information for organism '", org, "', download?\n")
        flag <- readline("yes/no(y/n)?")
        if (tolower(flag) %in% c("yes", "y")) {
            update_org_pathway(org = org)
            load(path_file, envir = envir)
        } else {
            stop("No pathway information for organism '", org, "', please use `update_org_pathway('", org, "')` to download.")
        }
    }

    if (verbose) {
        pcutils::dabiao("load ", org, " pathway")
        a <- get(paste0(org, "_kegg_pathway"), envir = envir)
        if (!is.null(attributes(a)$"download_time")) {
            pcutils::dabiao(paste0(org, " pathway download time: ", attributes(a)$"download_time"))
            message("If you want to update ", org, " pathway,", " use `update_org_pathway('", org, "')`")
        }
    }
}


#' Update the GO2gene files (from GO database)
#'
#'
#' @param download_dir download_dir
#' @param GO_file GO_file
#'
#' @export
#' @return No value
#' @description
#' Download links:
#' \code{http://geneontology.org/docs/download-ontology/}
#' \code{https://asa12138.github.io/FileList/GOlist.rda}
#' @examples
#' message("The following example will download some files:")
#' \donttest{
#' # update_GO_file()
#' }
update_GO_file <- function(download_dir = NULL, GO_file = NULL) {
    if (is.null(download_dir)) download_dir <- "ReporterScore_temp_download"
    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    if (is.null(GO_file)) {
        pcutils::dabiao("Trying to download files from https://asa12138.github.io/FileList/GOlist.rda")
        ori_time <- getOption("timeout")
        on.exit(options(timeout = ori_time))

        dir.create(download_dir, recursive = TRUE)
        options(timeout = 300)
        tryCatch(expr = {
            download.file("https://asa12138.github.io/FileList/GOlist.rda", destfile = file.path(download_dir, "GOlist.rda"))
        }, error = function(e) {
            stop("Try download yourself from https://asa12138.github.io/FileList/GOlist.rda")
        })
        GO_file <- file.path(download_dir, "GOlist.rda")
    } else {
        if (!(file.exists(GO_file) & grepl("GOlist.rda", GO_file))) stop("Wrong file: ", GO_file)
    }
    if (file.exists(paste0(pack_dir, "/new_GOlist.rda"))) {
        load(paste0(pack_dir, "/new_GOlist.rda"), envir = environment())
        old_time <- attributes(GOlist)$download_time
        load(GO_file, envir = environment())
        new_time <- attributes(GOlist)$download_time
        if (!(new_time > old_time)) {
            message("There are no changes of GO_file, do not need to update.")
            return(invisible())
        }
    }
    file.copy(GO_file, paste0(pack_dir, "/new_GOlist.rda"), overwrite = TRUE)
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
}


#' Load the GOlist (from GO database)
#'
#' @param envir `.GlobalEnv` (default) or `environment()`
#' @param verbose logical
#'
#' @export
#' @return GOlist in `.GlobalEnv`
#' @examples
#' load_GOlist()
#'
load_GOlist <- function(envir = .GlobalEnv, verbose = TRUE) {
    prefix <- "GOlist"
    GOlist_file <- file.path(tools::R_user_dir("ReporterScore"), paste0("new_", prefix, ".rda"))

    if (file.exists(GOlist_file)) {
        load(GOlist_file, envir = envir)
    } else {
        message("Not find GOlist_file, please run `update_GO_file()` first!")
    }

    if (verbose) {
        pcutils::dabiao("load ", prefix)
        if (!is.null(attributes(get("GOlist", envir = envir))$"download_time")) {
            pcutils::dabiao(paste0(prefix, " download time: ", attributes(get("GOlist", envir = envir))$"download_time"))
            message("If you want to update ", prefix, ", use `update_GO_file()`")
        }
    }
}

#' update_GOinfo
#'
#' @param download_dir download_dir
#' @param obo_file obo_file from http://current.geneontology.org/ontology/go.obo
#'
#' @export
#' @rdname update_GO_file
update_GOinfo <- function(download_dir = NULL, obo_file = NULL) {
    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    if (is.null(download_dir)) download_dir <- "ReporterScore_temp_download"
    if (is.null(obo_file)) {
        pcutils::dabiao("Trying to download files from http://current.geneontology.org/ontology/go.obo")
        ori_time <- getOption("timeout")
        on.exit(options(timeout = ori_time))

        dir.create(download_dir, recursive = TRUE)
        options(timeout = 300)
        tryCatch(expr = {
            download.file("http://current.geneontology.org/ontology/go.obo", destfile = file.path(download_dir, "go.obo"))
        }, error = function(e) {
            stop("Try download yourself from http://current.geneontology.org/ontology/go.obo")
        })
        obo_file <- file.path(download_dir, "go.obo")
    } else {
        if (!(file.exists(obo_file) & grepl("go.obo", obo_file))) stop("Wrong file: ", obo_file)
    }

    # 读取文本文件
    data_text <- readLines(obo_file)

    # 初始化一个空数据框
    df <- data.frame(
        id = character(0), name = character(0), namespace = character(0),
        def = character(0), synonym = character(0), alt_id = character(0),
        stringsAsFactors = FALSE
    )

    # 初始化变量以存储当前条目的信息
    cnt <- 0

    # 遍历每一行文本
    for (line in data_text) {
        # 如果行为空白，表示当前条目结束，将其添加到数据框中
        if (grepl("\\[Term\\]", line)) {
            if (cnt > 0) df <- rbind(df, current_entry)
            cnt <- cnt + 1
            if (cnt %% 1000 == 0) pcutils::dabiao(cnt, "done, please wait")
            # 重置current_entry
            current_entry <- list(id = NA, name = NA, namespace = NA, def = NA, synonym = NA, alt_id = NA)
        } else {
            # 使用正则表达式提取信息
            if (grepl("^id:", line)) {
                current_entry$id <- gsub("^id: ", "", line)
            } else if (grepl("^name:", line)) {
                current_entry$name <- gsub("^name: ", "", line)
            } else if (grepl("^namespace:", line)) {
                current_entry$namespace <- gsub("^namespace: ", "", line)
            } else if (grepl("^def:", line)) {
                current_entry$def <- gsub("^def: \"(.+)\" .+$", "\\1", line)
            } else if (grepl("^synonym:", line)) {
                if (is.na(current_entry$synonym)) {
                    current_entry$synonym <- gsub("^synonym: \"(.+)\" .+$", "\\1", line)
                } else {
                    current_entry$synonym <- paste0(current_entry$synonym, ";", gsub("^synonym: \"(.+)\" .+$", "\\1", line))
                }
            } else if (grepl("^alt_id:", line)) {
                if (is.na(current_entry$alt_id)) {
                    current_entry$alt_id <- gsub("^alt_id: ", "", line)
                } else {
                    current_entry$alt_id <- paste0(current_entry$alt_id, ";", gsub("^alt_id: ", "", line))
                }
            }
        }
    }
    colnames(df)[3] <- "ONT"

    GOinfo <- df
    attributes(GOinfo)$download_time <- Sys.time()
    save(GOinfo, file = paste0(pack_dir, "/new_GOinfo.rda"))
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
}

#' Load the GOinfo (from GO)
#'
#' @rdname load_GOlist
#'
#' @export
#' @return GOinfo in `.GlobalEnv`
load_GOinfo <- function(envir = .GlobalEnv, verbose = TRUE) {
    prefix <- "GOinfo"
    GOinfo_file <- file.path(tools::R_user_dir("ReporterScore"), paste0("new_", prefix, ".rda"))

    if (file.exists(GOinfo_file)) {
        load(GOinfo_file, envir = envir)
    } else {
        stop("use `update_GOinfo()` first")
    }

    if (verbose) {
        pcutils::dabiao("load ", prefix)
        if (!is.null(attributes(get("GOinfo", envir = envir))$"download_time")) {
            pcutils::dabiao(paste0(prefix, " download time: ", attributes(get("GOinfo", envir = envir))$"download_time"))
            message("If you want to update ", prefix, ", use `update_GOinfo()`")
        }
    }
}

# metacyc数据库



# CARD数据库

#' update CARDinfo from (from CARD database)
#'
#' @param download_dir download_dir
#' @param card_data card_data from https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2
#'
#' @export
#' @return No value
#' @examples
#' message("The following example will download some files:")
#' \donttest{
#' # update_CARDinfo()
#' }
update_CARDinfo <- function(download_dir = NULL, card_data = NULL) {
    pack_dir <- tools::R_user_dir("ReporterScore")
    if (!dir.exists(pack_dir)) dir.create(pack_dir, recursive = TRUE)

    if (is.null(download_dir)) download_dir <- "ReporterScore_temp_download"

    if (is.null(card_data)) {
        pcutils::dabiao("Trying to download files from https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2 ")
        ori_time <- getOption("timeout")
        on.exit(options(timeout = ori_time))

        dir.create(download_dir, recursive = TRUE)
        options(timeout = 300)
        tryCatch(expr = {
            download.file("https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2", destfile = file.path(download_dir, "card-data.tar.bz2"))
        }, error = function(e) {
            stop("Try download yourself from https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2 ")
        })
        card_data <- file.path(download_dir, "card-data.tar.bz2")
    } else {
        if (!(file.exists(card_data) & grepl("card-data.tar.bz2", card_data))) stop("Wrong file: ", card_data)
    }

    R.utils::bunzip2(card_data, destname = file.path(download_dir, "card-data.tar"), remove = FALSE, overwrite = TRUE)
    untar(tarfile = file.path(download_dir, "card-data.tar"), exdir = file.path(download_dir, "card-data"))

    {
        ARO_index <- readr::read_delim(file.path(download_dir, "card-data", "aro_index.tsv"), delim = "\t", progress = FALSE) %>% as.data.frame()
    } %>% suppressMessages()
    {
        antibiotics <- readr::read_delim(file.path(download_dir, "card-data", "shortname_antibiotics.tsv"), delim = "\t") %>% as.data.frame()
    } %>% suppressMessages()

    ARO_index <- dplyr::distinct(ARO_index, `ARO Accession`, .keep_all = TRUE)
    rownames(ARO_index) <- gsub("ARO:", "", (ARO_index$`ARO Accession`))
    CARDinfo <- list(ARO_index = ARO_index, antibiotics = antibiotics)
    attributes(CARDinfo)$download_time <- Sys.time()
    save(CARDinfo, file = paste0(pack_dir, "/new_CARDinfo.rda"))
    pcutils::dabiao(paste0("Update done at ", Sys.time()))
}

#' Load the CARDinfo (from CARD database)
#'
#' @param envir `.GlobalEnv` (default) or `environment()`
#' @param verbose logical
#'
#' @export
#' @return CARDinfo in `.GlobalEnv`
#' @examples
#' load_CARDinfo()
load_CARDinfo <- function(envir = .GlobalEnv, verbose = TRUE) {
    prefix <- "CARDinfo"
    GOinfo_file <- file.path(tools::R_user_dir("ReporterScore"), paste0("new_", prefix, ".rda"))

    if (file.exists(GOinfo_file)) {
        load(GOinfo_file, envir = envir)
    } else {
        stop("use `update_CARDinfo()` first")
    }

    if (verbose) {
        pcutils::dabiao("load ", prefix)
        if (!is.null(attributes(CARDinfo)$"download_time")) {
            pcutils::dabiao(paste0(prefix, " download time: ", attributes(CARDinfo)$"download_time"))
            message("If you want to update ", prefix, ", use `update_CARDinfo()`")
        }
    }
}
