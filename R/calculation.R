#' One step to get the reporter score of your KO abundance table.
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples.
#' @param group The comparison groups (at least two categories) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf. And you can use factor levels to change order.
#' @param metadata sample information data.frame contains group
#' @param mode 'mixed' or 'directed' (default, only for two groups differential analysis or multi-groups correlation analysis.), see details in \code{\link{pvalue2zs}}.
#' @param verbose logical
#' @param method the type of test. Default is `wilcox.test`. Allowed values include:
#' \itemize{
#' \item \code{\link[stats]{t.test}} (parametric) and \code{\link[stats]{wilcox.test}} (non-parametric). Perform comparison between two groups of samples. If the grouping variable contains more than two levels, then a pairwise comparison is performed.
#' \item \code{\link[stats]{anova}} (parametric) and \code{\link[stats]{kruskal.test}} (non-parametric). Perform one-way ANOVA test comparing multiple groups.
#' \item 'pearson', 'kendall', or 'spearman' (correlation), see \code{\link[stats]{cor}}.}
#' @param pattern a named vector matching the group, e.g. c('G1'=1,'G2'=3,'G3'=2), use the correlation analysis with specific pattern to calculate p-value.
#' @param feature one of 'ko', 'gene', 'compound'
#' @param type 'pathway' or 'module' for default KOlist for microbiome, 'CC', 'MF', 'BP', 'ALL' for default GOlist for homo sapiens. And org in listed in 'https://www.genome.jp/kegg/catalog/org_list.html' such as 'hsa' (if your kodf is come from a specific organism, you should specify type here).
#' @param threads default 1
#' @param p.adjust.method1 p.adjust.method for `ko.test`, see \code{\link[stats]{p.adjust}}
#' @param p.adjust.method2 p.adjust.method for the correction of ReporterScore, see \code{\link[stats]{p.adjust}}
#' @param modulelist NULL or customized modulelist dataframe, must contain 'id','K_num','KOs','Description' columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param perm permutation number, default: 4999.
#' @param min_exist_KO min exist KO number in a pathway (default, 3, when a pathway contains KOs less than 3, there will be no RS)
#' @param max_exist_KO max exist KO number in a pathway (default, 600, when a pathway contains KOs more than 600, there will be no RS)
#'
#' @aliases GRSA
#' @aliases RSA
#'
#' @return reporter_score object:
#' \item{kodf}{your input KO_abundance table}
#' \item{ko_stat}{ko statistics result contains p.value and z_score}
#' \item{reporter_s}{the reporter score in each pathway}
#' \item{modulelist}{default KOlist or customized modulelist dataframe}
#' \item{group}{The comparison groups in your data}
#' \item{metadata}{sample information dataframe contains group}
#'
#' for the `reporter_s` in result, whose columns represent:
#' \item{ID}{pathway id}
#' \item{Description}{pathway description}
#' \item{K_num}{total number of KOs/genes in the pathway}
#' \item{Exist_K_num}{number of KOs/genes in your inputdata that exist in the pathway}
#' \item{Significant_K_num}{number of kos/genes in your inputdata that are significant in the pathway}
#' \item{Z_score}{\eqn{Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}}}
#' \item{BG_Mean}{Background mean, \eqn{\mu _k}}
#' \item{BG_Sd}{Background standard deviation, \eqn{\sigma _k}}
#' \item{ReporterScore}{ReporterScore of the pathway, \eqn{ReporterScore=(Z_{pathway}-\mu _k)/\sigma _k}}
#' \item{p.value}{p.value of the ReporterScore}
#' \item{p.adjust}{adjusted p.value by p.adjust.method2}
#'
#' @family GRSA
#' @export
#' @examples
#' data("KO_abundance_test")
#' reporter_score_res <- reporter_score(KO_abundance, "Group", metadata,
#'   mode = "directed", perm = 499
#' )
#' head(reporter_score_res$reporter_s)
#' message("The following example require some time to run:")
#' \donttest{
#' reporter_score_res2 <- reporter_score(KO_abundance, "Group2", metadata,
#'   mode = "mixed",
#'   method = "kruskal.test", p.adjust.method1 = "none", perm = 499
#' )
#' reporter_score_res3 <- reporter_score(KO_abundance, "Group2", metadata,
#'   mode = "directed",
#'   method = "pearson", pattern = c("G1" = 1, "G2" = 3, "G3" = 2), perm = 499
#' )
#' }
reporter_score <- function(
    kodf, group, metadata = NULL, method = "wilcox.test", pattern = NULL,
    p.adjust.method1 = "none", mode = c("directed", "mixed")[1], verbose = TRUE,
    feature = "ko", type = c("pathway", "module")[1], p.adjust.method2 = "BH",
    modulelist = NULL, threads = 1, perm = 4999, min_exist_KO = 3, max_exist_KO = 600) {
  check_kodf_modulelist(kodf, type, feature, modulelist, verbose, mode = 1)

  stopifnot(mode %in% c("mixed", "directed"))
  method <- match.arg(method, c("t.test", "wilcox.test", "kruskal.test", "anova", "pearson", "kendall", "spearman", "lm"))
  if (!is.null(pattern)) {
    if (!method %in% c("pearson", "kendall", "spearman", "lm")) {
      stop("method should be one of \"pearson\", \"kendall\", \"spearman\" when you specify a pattern")
    }
    if (mode != "directed") {
      stop("method should be \"directed\" when you specify a pattern")
    }
  }

  if (!is.null(metadata)) {
    if (length(group) != 1) {
      stop("'group' should be one column name of metadata when metadata exsit!")
    }
    idx <- rownames(metadata) %in% colnames(kodf)
    metadata <- metadata[idx, , drop = FALSE]
    kodf <- kodf[, rownames(metadata), drop = FALSE]
    if (verbose) {
      message(nrow(metadata), " samples are matched for next step.")
    }
    if (length(idx) < 2) {
      stop("too less common samples")
    }
    sampFile <- data.frame(group = metadata[, group], row.names = row.names(metadata))
  } else {
    if (length(group) != ncol(kodf)) {
      stop("'group' length should equal to columns number of kodf when metadata is NULL!")
    }
    sampFile <- data.frame(row.names = colnames(kodf), group = group)
  }

  if (verbose) {
    pcutils::dabiao("Removing all-zero rows: ", sum(rowSums(abs(kodf)) == 0))
  }
  kodf <- kodf[rowSums(abs(kodf)) > 0, ]

  if (verbose) {
    pcutils::dabiao("1.KO test")
  }
  ko_pvalue <- ko.test(kodf, group, metadata, method = method, pattern = pattern, threads = threads, p.adjust.method1 = p.adjust.method1, verbose = verbose)
  if (verbose) {
    pcutils::dabiao("2.Transfer p.value to z-score")
  }
  ko_stat <- pvalue2zs(ko_pvalue, mode = mode, p.adjust.method1 = p.adjust.method1)
  attributes(ko_stat)$check <- TRUE
  if (verbose) {
    pcutils::dabiao("3.Calculating reporter score")
  }
  reporter_s <- get_reporter_score(ko_stat,
    type = type, feature = feature, threads = threads, p.adjust.method2 = p.adjust.method2, modulelist = modulelist, perm = perm,
    verbose = verbose, min_exist_KO = min_exist_KO, max_exist_KO = max_exist_KO
  )
  if (verbose) {
    pcutils::dabiao("All done")
  }

  if (is.null(modulelist)) {
    modulelist <- get_modulelist(type, feature, verbose = FALSE, chr = TRUE)
  }

  res <- list(kodf = kodf, ko_stat = ko_stat, reporter_s = reporter_s, modulelist = modulelist, group = group, metadata = metadata)
  class(res) <- "reporter_score"
  res
}


#' Differential analysis or Correlation analysis for KO-abundance table
#'
#' @inheritParams reporter_score
#'
#' @return ko_pvalue data.frame
#' @export
#' @family GRSA
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' ko_pvalue <- ko.test(KO_abundance, "Group", metadata)
#' }
ko.test <- function(kodf, group, metadata = NULL, method = "wilcox.test", pattern = NULL,
                    p.adjust.method1 = "none", threads = 1, verbose = TRUE) {
  i <- NULL
  t1 <- Sys.time()

  method <- match.arg(method, c("t.test", "wilcox.test", "kruskal.test", "anova", "pearson", "kendall", "spearman", "lm"))
  if (!is.null(pattern)) {
    if (!method %in% c("pearson", "kendall", "spearman")) {
      stop("method should be one of \"pearson\", \"kendall\", \"spearman\" when you specify a pattern")
    }
  }

  if (verbose) {
    pcutils::dabiao("Checking group")
  }
  if (!is.null(metadata)) {
    if (length(group) != 1) {
      stop("'group' should be one column name of metadata when metadata exsit!")
    }
    idx <- rownames(metadata) %in% colnames(kodf)
    metadata <- metadata[idx, , drop = FALSE]
    kodf <- kodf[, rownames(metadata), drop = FALSE]
    if (verbose) {
      message(nrow(metadata), " samples are matched for next step.")
    }
    if (length(idx) < 2) {
      stop("too less common samples")
    }
    sampFile <- data.frame(group = metadata[, group], row.names = row.names(metadata))
  } else {
    if (length(group) != ncol(kodf)) {
      stop("'group' length should equal to columns number of kodf when metadata is NULL!")
    }
    sampFile <- data.frame(row.names = colnames(kodf), group = group)
  }

  if (verbose) {
    pcutils::dabiao("Removing all-zero rows: ", sum(rowSums(abs(kodf)) == 0))
  }
  kodf <- kodf[rowSums(abs(kodf)) > 0, ]
  # calculate each
  if (verbose) {
    pcutils::dabiao("Calculating each KO")
  }
  if (verbose) {
    pcutils::dabiao("Using method: ", method)
  }
  tkodf <- t(kodf) %>%
    as.data.frame()
  group <- sampFile$group
  group2 <- NULL
  if (method %in% c("pearson", "kendall", "spearman", "lm")) {
    if (is.null(pattern)) {
      if (verbose) {
        message("Using correlation analysis: ", method, ", the groups will be transform to numeric, note the factor feature of group.")
      }
      if (is.numeric(group)) {
        group2 <- group
      } else {
        group2 <- as.numeric(factor(group))
      }
    } else {
      if (is.numeric(group)) {
        # stop('Can not use 'pattern' when group is a numeric variable.')
        if ((!is.numeric(pattern)) | length(group) != length(pattern)) {
          stop("pattern should be a numeric vector whose length equal to the group, e.g. c(1,2,3,5,3)")
        }
      } else {
        if ((!is.numeric(pattern)) | (is.null(names(pattern))) | (length(pattern) != nlevels(factor(group)))) {
          stop("pattern should be a named numeric vector matching the group, e.g. c(\"G1\"=1.5,\"G2\"=0.5,\"G3\"=2)")
        }
        group2 <- pattern[group]
      }
      if (verbose) {
        message("Using pattern")
      }
    }
  }

  res.dt <- data.frame(KO_id = rownames(kodf), row.names = rownames(kodf))

  if (is.numeric(sampFile$group)) {
    # stop('group should be a category variable.')
    vs_group <- "Numeric variable"
    if (!method %in% c("pearson", "kendall", "spearman", "lm")) {
      stop("group is a numeric variable, try \"pearson\", \"kendall\", \"spearman\" method")
      res.dt$cor <- cor(tkodf, group2, method = method)[, 1]
    }
  } else {
    vs_group <- levels(factor(sampFile$group))
    if (length(vs_group) == 1) {
      stop("'group' should be at least two elements factor")
    }
    if (length(vs_group) > 2) {
      if (method %in% c("t.test", "wilcox.test")) {
        stop("group\" more than two elements, try \"kruskal.test\", \"anova\" or \"pearson\", \"kendall\", \"spearman\"")
      }
    }

    for (i in vs_group) {
      tmpdf <- data.frame(average = apply(kodf[, which(group == i)], 1, mean), sd = apply(kodf[, which(group == i)], 1, sd))
      colnames(tmpdf) <- paste0(colnames(tmpdf), "_", i)
      res.dt <- cbind(res.dt, tmpdf)
    }
    if (length(vs_group) == 2) {
      # update, make sure the control group is first one.
      res.dt$diff_mean <- res.dt[, paste0("average_", vs_group[2])] - res.dt[, paste0("average_", vs_group[1])]
    }
    high_group <- apply(res.dt[, paste0("average_", vs_group)], 1, function(a) which(a == max(a))[[1]])
    res.dt$Highest <- vs_group[high_group]
  }
  if (method %in% c("pearson", "kendall", "spearman")) {
    res.dt$cor <- cor(tkodf, group2, method = method)[, 1]
  }
  # parallel
  reps <- nrow(kodf)

  tkodf <- t(kodf)
  # main function
  loop <- function(i) {
    val <- tkodf[, i]
    if (method == "wilcox.test") {
      pval <- stats::wilcox.test(val ~ group)$p.value
    }
    if (method == "t.test") {
      if (sd(val) == 0) {
        pval <- 1
      } else {
        pval <- stats::t.test(val ~ group)$p.value
      }
    }
    if (method == "kruskal.test") {
      pval <- stats::kruskal.test(val ~ group)$p.value
    }
    if (method == "anova") {
      pval <- stats::lm(val ~ group) %>%
        stats::anova(.) %>%
        .$`Pr(>F)` %>%
        .[1]
    }
    if (method == "lm") {
      # pval <- stats::lm(val~group) %>% stats::anova(.) %>%.$`Pr(>F)` %>% .[1]
      lm_res <- stats::lm(val ~ group)
      lm_res_summ <- summary(lm_res)
      lm_res_anova <- stats::anova(lm_res)
      r2 <- lm_res_summ %>%
        .$adj.r.squared %>%
        .[1]
      pval <- c(lm_res$coefficients[2], r2, as.numeric(lm_res_anova[1, ]))
    }
    if (method %in% c("pearson", "kendall", "spearman")) {
      # 用p-value还是直接效应量呢？?  Pval <- 1-abs(stats::cor(val,group2,method = method))
      pval <- stats::cor.test(val, group2, method = method)$p.value
    }
    if (verbose & (i %% 1000 == 0)) {
      message(i, " features done.")
    }
    if (is.na(pval)) {
      pval <- 1
    }
    pval
  }

  {
    if (threads > 1) {
      pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
      if (verbose) {
        pb <- utils::txtProgressBar(max = reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      } else {
        opts <- NULL
      }
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::`%dopar%`(
        foreach::foreach(
          i = seq_len(reps), .options.snow = opts,
          .export = c("group", "group2")
        ),
        suppressWarnings(loop(i))
      )
      snow::stopCluster(cl)
      gc()
    } else {
      res <- suppressWarnings(lapply(seq_len(reps), loop))
    }
  }
  # simplify method
  if (method == "lm") {
    res <- do.call(rbind, res)
    colnames(res) <- c("b-coefficient", "adj-r2", "df", "sum_sq", "mean_sq", "F-value", "p.value")
    res.dt <- cbind(res.dt, res)
  } else {
    res <- do.call(c, res)
    res.dt$p.value <- res
  }

  t2 <- Sys.time()
  stime <- sprintf("%.3f", t2 - t1)

  res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method1)

  resinfo <- paste0(
    "\nCompared groups: ", paste(vs_group, collapse = ", "), "\n", "Total KO number: ", reps, "\n", "Compare method: ", method, "\n", "Time use: ",
    stime, attr(stime, "units"), "\n"
  )
  if (verbose) message(resinfo)

  attributes(res.dt)$vs_group <- vs_group
  attributes(res.dt)$method <- method
  attributes(res.dt)$p.adjust.method1 <- p.adjust.method1
  attributes(res.dt)$pattern <- pattern
  return(res.dt)
}

#' Transfer p-value of KOs to Z-score
#'
#' @param ko_pvalue data.frame from \code{\link{ko.test}}
#' @inheritParams reporter_score
#'
#' @return ko_stat data.frame
#' @export
#' @details
#' '\strong{mixed}' mode is the original reporter-score method from Patil, K. R. et al. PNAS 2005.
#' In this mode, the reporter score is \strong{undirected}, and the larger the reporter score, the more significant the enrichment, but it cannot indicate the up-and-down regulation information of the pathway! (Liu, L. et al. iMeta 2023.)
#'
#' steps:
#'
#' 1. Use the Wilcoxon rank sum test to obtain the P value of the significance of each KO difference between the two groups (ie \eqn{P_{koi}}, i represents a certain KO);
#'
#' 2. Using an inverse normal distribution, convert the P value of each KO into a Z value (\eqn{Z_{koi}}), the formula:
#'
#' \eqn{Z_{koi}=\theta ^{-1}(1-P_{koi})}
#'
#' 3. 'Upgrade' KO to pathway: \eqn{Z_{koi}}, calculate the Z value of the pathway, the formula:
#'
#' \eqn{Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}}
#'
#' where k means A total of k KOs were annotated to the corresponding pathway;
#'
#' 4. Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of \eqn{Z_{pathway}}, the formula:
#'
#' \eqn{Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k}
#'
#' \eqn{\mu _k} is The mean of the random distribution, \eqn{\sigma _k} is the standard deviation of the random distribution.
#'
#' Instead, '\strong{directed}' mode is a derived version of 'mixed', referenced from \code{https://github.com/wangpeng407/ReporterScore}.
#'
#' This approach is based on the same assumption of many differential analysis methods: the expression of most genes has no significant change.
#'
#' steps:
#'
#' 1. Use the Wilcoxon rank sum test to obtain the P value of the significance of each KO difference between the two groups (ie \eqn{P_{koi}}, i represents a certain KO), and then divide the P value by 2, that is, the range of (0,1] becomes (0,0.5], \eqn{P_{koi}=P_{koi}/2};
#'
#' 2. Using an inverse normal distribution, convert the P value of each KO into a Z value (\eqn{Z_{koi}}), the formula:
#'
#' \eqn{Z_{koi}=\theta ^{-1}(1-P_{koi})}
#'
#' since the above P value is less than 0.5, all Z values will be greater than 0;
#'
#' 3. Considering whether each KO is up-regulated or down-regulated, calculate \eqn{diff\_KO},
#'
#' \eqn{Z_{koi}=-Z_{koi}\ \ \ \ (diff\_KO<0)},
#'
#' so \eqn{Z_{koi}} is greater than 0 Up-regulation, \eqn{Z_{koi}} less than 0 is down-regulation;
#'
#' 4. 'Upgrade' KO to pathway: \eqn{Z_{koi}}, calculate the Z value of the pathway, the formula:
#'
#'  \eqn{Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}}
#'
#'  where k means A total of k KOs were annotated to the corresponding pathway;
#'
#' 5. Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of \eqn{Z_{pathway}}, the formula:
#'
#' \eqn{Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k}
#'
#' \eqn{\mu _k} is The mean of the random distribution, \eqn{\sigma _k} is the standard deviation of the random distribution.
#'
#' The finally obtained \eqn{Z_{adjustedpathway}} is the Reporter score value enriched for each pathway.
#' In this mode, the Reporter score is directed, and a larger positive value represents a significant up-regulation enrichment, and a smaller negative values represent significant down-regulation enrichment.
#'
#' However, the disadvantage of this mode is that when a pathway contains about the same number of significantly up-regulates KOs and significantly down-regulates KOs, the final absolute value of Reporter score may approach 0, becoming a pathway that has not been significantly enriched.
#'
#'
#' @references
#' 1. Patil, K. R. & Nielsen, J. Uncovering transcriptional regulation of metabolism by using metabolic network topology. Proc Natl Acad Sci U S A 102, 2685–2689 (2005).
#' 2. Liu, L., Zhu, R. & Wu, D. Misuse of reporter score in microbial enrichment analysis. iMeta n/a, e95.
#' 3. \code{https://github.com/wangpeng407/ReporterScore}
#' @family GRSA
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' ko_pvalue <- ko.test(KO_abundance, "Group", metadata)
#' ko_stat <- pvalue2zs(ko_pvalue, mode = "directed")
#' }
pvalue2zs <- function(ko_pvalue, mode = c("directed", "mixed")[1], p.adjust.method1 = "none") {
  p.adjust <- type <- Significantly <- NULL
  res.dt <- ko_pvalue

  if ("origin_p.value" %in% colnames(res.dt)) {
    message("detect the origin_p.value, use the origin_p.value")
    res.dt$p.value <- res.dt$origin_p.value
  }
  if (!all(c("p.value") %in% colnames(res.dt))) {
    stop("check if `p.value` in your ko_stat dataframe!")
  }

  stopifnot(mode %in% c("mixed", "directed"))

  if ("diff_mean" %in% colnames(res.dt)) {
    res.dt$sign <- ifelse(res.dt$diff_mean < 0, -1, 1)
    res.dt$type <- ifelse(res.dt$diff_mean < 0, paste0("Depleted"), paste0("Enriched"))
  }

  if ("cor" %in% colnames(res.dt)) {
    res.dt$sign <- ifelse(res.dt$cor < 0, -1, 1)
    if (is.null(attributes(ko_pvalue)$pattern)) {
      res.dt$type <- ifelse(res.dt$cor < 0, paste0("Depleted"), paste0("Enriched"))
    } else {
      res.dt$type <- ifelse(res.dt$cor < 0, paste0("None"), paste0("Positive"))
    }
  }
  if ("b-coefficient" %in% colnames(res.dt)) {
    res.dt$sign <- ifelse(res.dt$`b-coefficient` < 0, -1, 1)
    if (is.null(attributes(ko_pvalue)$pattern)) {
      res.dt$type <- ifelse(res.dt$`b-coefficient` < 0, paste0("Depleted"), paste0("Enriched"))
    } else {
      res.dt$type <- ifelse(res.dt$`b-coefficient` < 0, paste0("None"), paste0("Positive"))
    }
  }

  # mixed不考虑正负号，p.adjust不除以2，考虑的话除以2
  if (mode == "mixed") {
    res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method1)
    if (quantile(res.dt$p.adjust, 0.1, na.rm = TRUE) == 1) {
      warning("most of p.adjust are 1! try use the p.adjust.method1='none'\n")
    }

    # 逆正态分布
    zs <- stats::qnorm(1 - (res.dt$p.adjust))
    # 防止太小的p.adjust产生Inf
    zs <- ifelse(zs > 8.209536, 8.209536, zs)
    zs <- ifelse(zs < -8.209536, -8.209536, zs)
    res.dt$Z_score <- zs
    attributes(res.dt)$mode <- "mixed"
  }
  if (mode == "directed") {
    if (!"type" %in% colnames(res.dt)) {
      stop("directed mode only use for two groups differential analysis or multi-groups correlation analysis.")
    }
    res.dt$origin_p.value <- ko_pvalue$p.value
    if ("p.adjust" %in% colnames(res.dt)) {
      res.dt$origin_p.adjust <- ko_pvalue$p.adjust
    }
    pn_sign <- 2
    res.dt$p.value <- res.dt$p.value / pn_sign
    res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method1)
    if (quantile(res.dt$p.adjust, 0.1, na.rm = TRUE) == 0.5) {
      warning("most of p.adjust are 1! try use the p.adjust.method1='none'\n")
    }

    # 这种做法可能要基于一个前提，就是上下调ko数量基本一致,才能保证正负都是显著差异的，或者分开正负分析？
    # up_down_ratio=table(res.dt%>%dplyr::filter(abs(p.adjust)<=stats::quantile(res.dt$p.adjust,0.05,na.rm=TRUE))%>%dplyr::pull(type))
    # kafang_res=stats::chisq.test(up_down_ratio) pcutils::dabiao('') pcutils::dabiao('Chi-squared test for up and down ko ratio') message('X-squared =
    # ',round(kafang_res$statistic,4), ' p-value = ',round(kafang_res$p.value,4)) #if p-value>0.05，正负一致。 if(kafang_res$p.value<0.05){ message('The
    # overall up-down ratio of ko abundance is unbalanced!\n Continuing to use the directed mode may lead to wrong conclusions') }

    # 逆正态分布
    zs <- stats::qnorm(1 - (res.dt$p.adjust))
    # 防止太小的p.adjust产生Inf
    zs <- ifelse(zs > 8.209536, 8.209536, zs)
    zs <- ifelse(zs < -8.209536, -8.209536, zs)

    # 通过判断上下调给予z-score正负号，让最后的reporter-score正负号作为上下调标志
    res.dt$Z_score <- ifelse(res.dt$sign < 0, -zs, zs)
    attributes(res.dt)$mode <- "directed"
  }

  if (is.null(attributes(res.dt)$mode)) {
    p_th <- 0.05
  } else {
    p_th <- ifelse(attributes(res.dt)$mode == "directed", 0.025, 0.05)
  }
  attributes(res.dt)$pattern <- attributes(ko_pvalue)$pattern

  if ("type" %in% colnames(res.dt) & mode == "directed") {
    res.dt <- dplyr::mutate(res.dt, Significantly = ifelse(p.adjust < p_th, type, "None"))
    if (!is.null(attributes(ko_pvalue)$pattern)) {
      res.dt <- dplyr::mutate(res.dt, Significantly = ifelse(Significantly == "Positive", "Significant", Significantly))
    }
  } else {
    res.dt <- dplyr::mutate(res.dt, Significantly = ifelse(p.adjust < p_th, "Significant", "None"))
  }
  return(res.dt)
}

#' Calculate the mean and standard deviation of the random distribution
#'
#' @param vec numeric
#' @param Knum numeric
#' @param perm numeric
#'
#' @return list
#' @noRd
random_mean_sd <- function(vec, Knum, perm = 1000) {
  # Permutation就是不放回抽样😭，Bootstrap才是有放回

  replace <- FALSE
  temp <- vapply(seq_len(perm), \(i) {
    sum(sample(vec, Knum, replace = replace)) / sqrt(Knum)
  }, numeric(1))
  list(vec = temp, mean_sd = c(mean(temp), stats::sd(temp)))
}

#' Check the input of kodf and modulelist
#'
#' @param type "pathway", "module"
#' @param feature "ko", "compound", "gene"
#' @param gene one of "symbol","id"
#' @param verbose logical
#' @param chr keep chraacter or not
#'
#' @return modulelist
#' @noRd
get_modulelist <- function(type, feature, gene = "symbol", verbose = TRUE, chr = FALSE) {
  if (type %in% c("pathway", "module")) {
    # reference pathway
    type <- match.arg(type, c("pathway", "module"))
    if (feature == "ko") {
      KOlist <- load_KOlist(verbose = verbose)
      modulelist <- KOlist[[type]]
    }
    if (feature == "compound") {
      CPDlist <- load_CPDlist(verbose = verbose)
      modulelist <- CPDlist[[type]]
    }
    if (feature == "gene") {
      stop("If the feature of your table is \"gene\", please sepcify the \"type\" arugment as an org in listed in https://www.genome.jp/kegg/catalog/org_list.html or \"CC\", \"MF\", \"BP\", \"ALL\" for default GOlist")
    }
  } else if (type %in% c("CC", "MF", "BP", "ALL")) {
    if (feature != "gene") {
      stop("\"CC\", \"MF\", \"BP\", \"ALL\" using GO database, which only support feature=\"gene\"")
    }
    GOlist <- load_GOlist(verbose = verbose)
    if (type == "ALL") {
      modulelist <- "lapply(names(GOlist), \\(i) cbind(GOlist[[i]], ONT = i)) %>%
                    do.call(rbind, .)"
      if (!chr) modulelist <- eval(parse(text = modulelist))
    } else {
      modulelist <- paste0("GOlist[['", type, "']]")
      if (!chr) modulelist <- eval(parse(text = modulelist))
    }
  } else {
    # KEGG pathway of other organisms
    modulelist <- custom_modulelist_from_org(type, feature = feature, gene = gene, verbose = verbose)
    if (verbose) {
      message("You choose the feature: '", feature, "', make sure the rownames of your input table are right.")
    }
  }
  return(modulelist)
}

#' Calculate reporter score
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}}
#' @inheritParams reporter_score
#'
#' @return reporter_res data.frame
#' @export
#' @family GRSA
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' ko_pvalue <- ko.test(KO_abundance, "Group", metadata)
#' ko_stat <- pvalue2zs(ko_pvalue, mode = "directed")
#' reporter_s1 <- get_reporter_score(ko_stat, perm = 499)
#' }
get_reporter_score <- function(
    ko_stat, type = c("pathway", "module")[1], feature = "ko",
    threads = 1, modulelist = NULL, perm = 4999,
    verbose = TRUE, p.adjust.method2 = "BH",
    min_exist_KO = 3, max_exist_KO = 600) {
  ReporterScore <- KOlist <- i <- NULL
  type_flag <- FALSE
  t1 <- Sys.time()

  check_kodf_modulelist(ko_stat, type, feature, modulelist, verbose, mode = 2)

  if (is.null(modulelist)) {
    modulelist <- get_modulelist(type, feature, verbose)
    type_flag <- TRUE
  }

  # calculate each pathway
  if (verbose) {
    if (type %in% c("module")) {
      pcutils::dabiao("Calculating each module")
    }
    if (type %in% c("CC", "MF", "BP", "ALL")) {
      pcutils::dabiao("Calculating each GO term, please wait as there are lots of GO terms")
    } else {
      pcutils::dabiao("Calculating each pathway")
    }
  }

  # parallel
  reps <- nrow(modulelist)

  if (is.null(attributes(ko_stat)$mode)) {
    p_th <- 0.05
  } else {
    p_th <- ifelse(attributes(ko_stat)$mode == "directed", 0.025, 0.05)
  }

  clean.KO <- ko_stat$Z_score[!is.na(ko_stat$Z_score)]
  # main function
  loop <- function(i) {
    # 找到在该pathway里的所有ko的zs
    tmp_kos <- strsplit(modulelist$KOs[i], ",")[[1]]

    z <- ko_stat[ko_stat$KO_id %in% tmp_kos, ]
    exist_KO <- nrow(z)

    # significant_KO=sum(z$p.adjust<p_th)
    significant_KO <- sum(z$Significantly != "None")
    sig_up <- sum(z$Significantly == "Enriched")
    sig_down <- sum(z$Significantly == "Depleted")

    # 如果一条通路里压根没找到几个ko，就不应该有reporterscore
    if ((exist_KO < min_exist_KO) | (exist_KO > max_exist_KO)) {
      if (attributes(ko_stat)$mode == "directed") {
        return(c(exist_KO, significant_KO, sig_up, sig_down, NA, NA, NA, NA, NA))
      } else {
        return(c(exist_KO, significant_KO, NA, NA, NA, NA, NA))
      }
    }
    # KOnum <- modulelist$K_num[i] KOnum <- ifelse(length(clean.KO) >= KOnum, KOnum, length(clean.KO))

    KOnum <- exist_KO
    # 以整个输入ko文件作为背景,抽取KOnum应该是exist_KO，而不是所有的KOnum，可以在iMeta文章看到
    mean_sd <- random_mean_sd(clean.KO, KOnum, perm = perm)
    Z_score <- sum(z$Z_score, na.rm = TRUE) / sqrt(KOnum)

    reporter_score <- (Z_score - mean_sd$mean_sd[1]) / mean_sd$mean_sd[2]

    if (attributes(ko_stat)$mode == "directed") {
      if (reporter_score > 0) {
        p.value <- (sum(mean_sd$vec >= Z_score) + 1) / (length(mean_sd$vec) + 1)
      } else {
        p.value <- (sum(mean_sd$vec <= Z_score) + 1) / (length(mean_sd$vec) + 1)
      }
    } else {
      p.value <- (sum(mean_sd$vec >= Z_score) + 1) / (length(mean_sd$vec) + 1)
    }

    if (verbose & (i %% 50 == 0)) {
      message(i, " pathways done.")
    }
    if (attributes(ko_stat)$mode == "directed") {
      return(c(exist_KO, significant_KO, sig_up, sig_down, Z_score, mean_sd$mean_sd, reporter_score, p.value))
    } else {
      return(c(exist_KO, significant_KO, Z_score, mean_sd$mean_sd, reporter_score, p.value))
    }
  }
  {
    if (threads > 1) {
      pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
      if (verbose) {
        pb <- utils::txtProgressBar(max = reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      } else {
        opts <- NULL
      }
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::`%dopar%`(
        foreach::foreach(
          i = seq_len(reps), .options.snow = opts,
          .export = c("random_mean_sd")
        ),
        loop(i)
      )
      snow::stopCluster(cl)
      gc()
    } else {
      res <- lapply(seq_len(reps), loop)
    }
  }
  # simplify method
  res <- do.call(rbind, res)
  if (attributes(ko_stat)$mode == "directed") {
    colnames(res) <- c("Exist_K_num", "Significant_K_num", "Significant_up_num", "Significant_down_num", "Z_score", "BG_Mean", "BG_Sd", "ReporterScore", "p.value")
  } else {
    colnames(res) <- c("Exist_K_num", "Significant_K_num", "Z_score", "BG_Mean", "BG_Sd", "ReporterScore", "p.value")
  }
  res <- as.data.frame(res)
  reporter_res <- data.frame(ID = modulelist$id, Description = modulelist$Description, K_num = modulelist$K_num, res)
  if (type == "ALL") {
    reporter_res$ONT <- modulelist$ONT
  }

  reporter_res$p.adjust <- stats::p.adjust(reporter_res$p.value, p.adjust.method2)
  if (attributes(ko_stat)$mode == "directed") {
    reporter_res <- dplyr::arrange(reporter_res, -abs(ReporterScore))
  } else {
    reporter_res <- dplyr::arrange(reporter_res, -(ReporterScore))
  }
  attributes(reporter_res)$mode <- attributes(ko_stat)$mode
  attributes(reporter_res)$method <- attributes(ko_stat)$method
  attributes(reporter_res)$vs_group <- attributes(ko_stat)$vs_group
  attributes(reporter_res)$pattern <- attributes(ko_stat)$pattern

  attributes(reporter_res)$perm <- perm

  if (type_flag) {
    attributes(reporter_res)$type <- type
    attributes(reporter_res)$feature <- feature
  }

  rownames(reporter_res) <- reporter_res$ID

  t2 <- Sys.time()

  stime <- sprintf("%.3f", t2 - t1)
  resinfo <- paste0("ID number: ", reps, "\n", "Time use: ", stime, attr(stime, "units"), "\n")
  if (verbose) {
    message(resinfo)
  }
  return(reporter_res)
}

#' Check the input of kodf and modulelist
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}}
#' @param type "pathway", "module"
#' @param feature "ko", "compound", "gene"
#' @param modulelist NULL or customized modulelist dataframe, must contain 'id','K_num','KOs','Description' columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#' @param mode 1 or 2
#'
#' @return invisible
#' @noRd
check_kodf_modulelist <- function(ko_stat, type, feature, modulelist, verbose, mode = 1) {
  if (!is.null(attributes(ko_stat)$check)) {
    if (attributes(ko_stat)$check) {
      return(invisible())
    }
  }
  if (mode == 1) {
  } else {
    if (!all(c("KO_id", "Z_score") %in% colnames(ko_stat))) {
      stop("Some wrong with ko_stat, please check if the colnames of ko_stat contains KO_id and Z_score?")
    }
    if (!all(rownames(ko_stat) == ko_stat$KO_id)) {
      stop("Rownames of ko_stat do not match the KO_id in ko_stat")
    }
  }

  if (verbose) pcutils::dabiao("Use feature: ", feature)
  if (verbose) pcutils::dabiao("Checking rownames")
  if (feature == "ko") {
    rowname_check <- grepl("K\\d{5}", rownames(ko_stat))
    if (!all(rowname_check)) {
      if (verbose) message("Some of your ko_stat are not KO id, check the format! (e.g. K00001)\n")
    }
  }
  if (feature == "gene") {
    if (verbose) message("please make sure your input table rows are gene symbol!\n")
  }
  if (feature == "compound") {
    rowname_check <- grepl("C\\d{5}", rownames(ko_stat))
    if (!all(rowname_check)) {
      if (verbose) message("Some of your ko_stat are not 'KEGG' compound id, check the format! (e.g. C00001)\n")
    }
  }

  if (is.null(modulelist)) {
    modulelist <- get_modulelist(type, feature, verbose = FALSE)
  }

  if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) {
    stop("Please check your modulelist format!")
  }

  all_kos <- unique(transform_modulelist(modulelist, 2)[, 2])
  ko_ratio <- sum(rownames(ko_stat) %in% all_kos) / nrow(ko_stat)

  if (ko_ratio == 0) {
    stop("All your ", feature, "s do not exist in the modulelist! Please check your input table and modulelist (or  parameters `feature` and `type`).")
  } else if (ko_ratio < 0.1) {
    warning("Only ", 100 * round(ko_ratio, 4), "% of your ", feature, "s in the modulelist! Please make sure your input table and modulelist (or  parameters `feature` and `type`) are right.")
  } else {
    if (verbose) message(100 * round(ko_ratio, 4), "% of your ", feature, "s in the modulelist!")
  }
}

#' get features in a modulelist
#'
#' @param map_id map_id in modulelist
#' @param ko_stat NULL or ko_stat result from \code{\link{pvalue2zs}}
#' @param modulelist NULL or customized modulelist dataframe, must contain 'id','K_num','KOs','Description' columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#'
#' @family modulelist
#' @aliases get_KOs
#' @export
#' @return KOids, or data.frame with these KOids.
#' @examples
#' get_features(map_id = "map00010")
get_features <- function(map_id = "map00010", ko_stat = NULL, modulelist = NULL) {
  KOlist <- GOlist <- NULL
  if (is.null(modulelist)) {
    message("modulelist is NULL, use default modulelist!")
    KOlist <- load_KOlist()
    if (grepl("map", map_id[1])) {
      modulelist <- KOlist$pathway
    }
    if (grepl("M", map_id[1])) {
      modulelist <- KOlist$module
    }
    if (grepl("GO:", map_id[1])) {
      GOlist <- load_GOlist()
      modulelist <- lapply(names(GOlist), function(i) cbind(GOlist[[i]], ONT = i)) %>%
        do.call(rbind, .)
    }
  }
  if (is.null(modulelist)) {
    stop("check your modulelist.")
  }

  if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) {
    stop("check your modulelist format!")
  }
  kos <- modulelist[which(modulelist$id %in% map_id), "KOs"]
  if (length(kos) > 0) {
    kos <- lapply(kos, strsplit, ",") %>%
      unlist()
  } else {
    kos <- NULL
  }
  if (!"KO_id" %in% colnames(ko_stat)) {
    rownames(ko_stat) -> ko_stat$KO_id
  }
  if (!is.null(ko_stat)) {
    return(ko_stat[ko_stat$KO_id %in% kos, ])
  }
  return(kos)
}


#' Build a custom modulelist
#'
#' @param pathway2ko user input annotation of Pathway to KO mapping, a data.frame of 2 column with pathway and ko.
#' @param pathway2desc user input of Pathway TO Description mapping, a data.frame of 2 column with pathway and description.
#' @param verbose verbose
#'
#' @family modulelist
#' @return a custom modulelist
#' @export
#' @examples
#' mydat <- data.frame(pathway = paste0("PATHWAY", rep(seq_len(2), each = 5)), ko = paste0("K", 1:10))
#' mymodulelist <- custom_modulelist(mydat)
#' print(mymodulelist)
#' transform_modulelist(mymodulelist)
custom_modulelist <- function(pathway2ko, pathway2desc = NULL, verbose = TRUE) {
  Pathway <- NULL
  pathway2ko <- pathway2ko[, c(1, 2)]
  colnames(pathway2ko) <- c("Pathway", "KO")
  pathway2ko_num <- pathway2ko %>% dplyr::count(Pathway, name = "K_num")

  if (is.null(pathway2desc)) {
    pathway2desc <- data.frame(Pathway = pathway2ko_num$Pathway, Description = pathway2ko_num$Pathway)
  } else {
    pathway2desc <- pathway2desc[, c(1, 2)]
  }
  colnames(pathway2desc) <- c("Pathway", "Description")

  pathway2ko_com <- stats::aggregate(KO ~ Pathway, pathway2ko, paste, collapse = ",")
  pathway_list <- Reduce(\(x, y)dplyr::left_join(x = x, y = y, by = "Pathway"), list(pathway2ko_num, pathway2ko_com, pathway2desc))
  colnames(pathway_list) <- c("id", "K_num", "KOs", "Description")
  if (verbose) message("please assgin this custom modulelist to `reporter_score(modulelist=your_modulelist)` to do a custom enrichment.")
  pathway_list$Description <- ifelse(is.na(pathway_list$Description), pathway_list$id, pathway_list$Description)
  pathway_list
}

#' Transform a modulelist to a list
#'
#' @param mymodulelist mymodulelist
#' @param mode 1~2
#'
#' @family modulelist
#' @export
#' @return modulelist
#' @rdname custom_modulelist
transform_modulelist <- function(mymodulelist, mode = 1) {
  if (mode == 1) {
    setNames(strsplit(mymodulelist$KOs, ","), mymodulelist$id) %>% return()
  } else {
    path2ko <- pcutils::explode(mymodulelist[, c("id", "KOs")], column = 2, split = ",")
    if (mode == 2) {
      return(path2ko)
    } else {
      path2ko_mat <- reshape2::acast(data.frame(path2ko, value = 1), KOs ~ id)
      path2ko_mat[is.na(path2ko_mat)] <- 0
      return(path2ko_mat)
    }
  }
}


#' Custom modulelist from a specific organism
#'
#' @param org kegg organism, listed in https://www.genome.jp/kegg/catalog/org_list.html, default, "hsa"
#' @param feature one of "ko", "gene", "compound"
#' @param verbose logical
#' @param gene one of "symbol","id"
#'
#' @return modulelist
#' @export
#'
#' @family modulelist
#' @examples
#' hsa_pathway <- custom_modulelist_from_org(org = "hsa", feature = "gene")
custom_modulelist_from_org <- function(org = "hsa", feature = "ko", gene = "symbol", verbose = TRUE) {
  org_path <- load_org_pathway(org = org, verbose = verbose)
  if (feature == "ko") {
    pathway2ko <- dplyr::distinct_all(org_path$all_org_gene[, c("pathway_id", "KO_id")])
  } else if (feature == "gene") {
    gene <- match.arg(gene, c("symbol", "id"))
    if (gene == "symbol") {
      pathway2ko <- dplyr::distinct_all(org_path$all_org_gene[, c("pathway_id", "gene_symbol")])
    } else if (gene == "id") pathway2ko <- dplyr::distinct_all(org_path$all_org_gene[, c("pathway_id", "kegg_gene_id")])
  } else if (feature == "compound") {
    pathway2ko <- dplyr::distinct_all(org_path$all_org_compound[, c("pathway_id", "kegg_compound_id")])
  } else {
    stop("feature should be one of 'ko', 'gene', 'compound'.")
  }
  org_modulist <- custom_modulelist(pathway2ko = pathway2ko, pathway2desc = org_path$org_pathway[, c("pathID", "org_pathway")], verbose = verbose)
  org_modulist
}

#' Modify the pathway description before plotting
#'
#' @param reporter_res reporter_res
#' @param pattern str, like " - Homo sapiens (human)"
#' @param replacement str, like ""
#'
#' @export
#' @return reporter_res
#' @examples
#' data("reporter_score_res")
#' modify_description(reporter_score_res, pattern = " - Homo sapiens (human)")
modify_description <- function(reporter_res, pattern = " - Homo sapiens (human)", replacement = "") {
  if (inherits(reporter_res, "reporter_score")) {
    reporter_res$reporter_s$Description <- gsub(pattern, replacement, reporter_res$reporter_s$Description, fixed = TRUE)
    reporter_res$modulelist$Description <- gsub(pattern, replacement, reporter_res$modulelist$Description, fixed = TRUE)
    return(reporter_res)
  }
  if (inherits(reporter_res, "rs_by_cm")) {
    rsa_cm_res <- reporter_res
    for (i in grep("Cluster", names(rsa_cm_res), value = TRUE)) {
      rsa_cm_res[[i]]$reporter_s$Description <- gsub(pattern, replacement, rsa_cm_res[[i]]$reporter_s$Description, fixed = TRUE)
    }
    rsa_cm_res$modulelist$Description <- gsub(pattern, replacement, rsa_cm_res$modulelist$Description, fixed = TRUE)
    return(rsa_cm_res)
  }
  if (is.data.frame(reporter_res)) {
    reporter_res$Description <- gsub(pattern, replacement, reporter_res$Description, fixed = TRUE)
    return(reporter_res)
  }
}

#' Upgrade the KO level
#'
#' @param KO_abundance KO_abundance
#' @param level one of 'pathway', 'module', 'level1', 'level2', 'level3', 'module1', 'module2', 'module3'.
#' @param show_name logical
#' @param modulelist NULL or customized modulelist dataframe, must contain 'id','K_num','KOs','Description' columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#'
#' @return data.frame
#' @export
#'
#' @examples
#' data("KO_abundance_test")
#' KO_level1 <- up_level_KO(KO_abundance, level = "level1", show_name = TRUE)
up_level_KO <- function(KO_abundance, level = "pathway", show_name = FALSE, modulelist = NULL, verbose = TRUE) {
  a <- KO_abundance
  a$KO_id <- rownames(a)

  if (level %in% c("pathway", "module")) {
    if (is.null(modulelist)) {
      KOlist <- load_KOlist(verbose = verbose)
      modulelist <- KOlist[[level]]
    }
    if (!all(c("id", "K_num", "KOs", "Description") %in% colnames(modulelist))) {
      stop("check your modulelist format!")
    }
    path2ko <- pcutils::explode(modulelist[, c("id", "KOs")], 2, split = ",")
    path2name <- setNames(modulelist$Description, modulelist$id)
  } else if (level %in% c("level1", "level2", "level3")) {
    KO_htable <- NULL
    KO_htable <- load_KO_htable(verbose = verbose)
    if (level == "level3") {
      path2ko <- KO_htable[, c(paste0(level, "_id"), "KO_id")]
      path2name <- KO_htable[, paste0(level, c("_id", "_name"))] %>%
        dplyr::distinct()
      path2name <- setNames(path2name[, paste0(level, "_name"), drop = TRUE], path2name[, paste0(level, "_id"), drop = TRUE])
    } else {
      path2ko <- KO_htable[, c(paste0(level, "_name"), "KO_id")]
      show_name <- FALSE
    }
  } else if (level %in% c("module1", "module2", "module3")) {
    Module_abundance <- up_level_KO(KO_abundance, level = "module")
    a <- Module_abundance
    a$KO_id <- rownames(a)
    Module_htable <- load_Module_htable(verbose = verbose)
    path2ko <- Module_htable[, c(paste0(level, "_name"), "Module_id")]
    show_name <- FALSE
  } else {
    stop("level should be one of \"pathway\", \"module\", \"level1\", \"level2\", \"level3\", \"module1\", \"module2\", \"module3\".")
  }

  colnames(path2ko) <- c("Pathway", "KO_id")
  path2ko <- dplyr::distinct_all(path2ko)
  dplyr::left_join(a, path2ko, by = "KO_id") -> aa
  aa$Pathway[is.na(aa$Pathway)] <- "Unknown"

  b <- pcutils::hebing(dplyr::select(aa, -c("KO_id", "Pathway")), aa$Pathway, 1, act = "sum")
  if (show_name) {
    rownames(b) <- c(path2name, Unknown = "Unknown")[rownames(b)]
  }
  b
}


#' Transfer gene symbol table to KO table
#' @description
#' You can use `clusterProfiler::bitr()` to transfer your table from other gene_id to gene_symbol.
#'
#' @param genedf ,rowname is gene symbol (e.g. PFKM), colnames is samples
#' @param org kegg organism, listed in 'https://www.genome.jp/kegg/catalog/org_list.html', default, 'hsa'
#'
#' @return kodf
#' @export
#'
#' @examples
#' data("genedf")
#' KOdf <- gene2ko(genedf, org = "hsa")
gene2ko <- function(genedf, org = "hsa") {
  org_path <- load_org_pathway(org = org)

  gene_mp_ko <- dplyr::distinct_all(org_path$all_org_gene[, c("gene_symbol", "KO_id")])

  gene_symbol <- data.frame(gene_symbol = rownames(genedf))
  in_genes <- length(intersect(rownames(genedf), gene_mp_ko$gene_symbol))
  out_genes <- nrow(genedf) - in_genes
  message("\n", in_genes, " genes are in the pathways, while ", out_genes, " genes are not.")
  gene_symbol <- dplyr::left_join(gene_symbol, gene_mp_ko, by = "gene_symbol")
  gene_symbol$KO_id <- ifelse(is.na(gene_symbol$KO_id), "others", gene_symbol$KO_id)

  KOdf <- pcutils::hebing(genedf[gene_symbol$gene_symbol, ], gene_symbol$KO_id, 1, "sum")
  KOdf
}

# 做一个函数进行cm_cluster，找出十分显著的pattern 然后再对每一个pattern做ko.test_with_pattern，即可以得到每个pattern的RSA结果

if (FALSE) {
  pcutils::hebing(KO_abundance, metadata$Group2, 2) -> KO_abundance_Group2
  cm_test_k(KO_abundance_Group2, 0.7)
  cm_res <- c_means(KO_abundance_Group2, 3, 0.7)
  plot(cm_res, filter_membership = 0.8, mode = 1)
  cm_res$centers

  test1 <- reporter_score(KO_abundance, "Group2", metadata, mode = "directed", method = "pearson", pattern = cm_res$centers[1, ])
  plot_report(test1, show_ID = TRUE)
  plot_KOs_in_pathway(test1, map_id = "map03010")
}



#' Test the proper clusters k for c_means
#'
#' @param otu_group grouped otutab
#' @param filter_var filter the highest var
#' @param fast whether do the gap_stat?
#'
#' @return ggplot
#' @export
#'
#' @rdname c_means
cm_test_k <- function(otu_group, filter_var, fast = TRUE) {
  data_scaled <- filter_top_var(otu_group, filter_var)

  # 判断聚类个数
  # 输入文件最好是按你想要的分组合并过的
  lib_ps("factoextra", library = FALSE)
  #-------determining the number of clusters
  # 1 Elbow method
  cp1 <- factoextra::fviz_nbclust(x = data_scaled, FUNcluster = stats::kmeans, method = "wss", verbose = TRUE) +
    labs(subtitle = "Elbow method")
  # 2 Silhouette method
  cp2 <- factoextra::fviz_nbclust(x = data_scaled, FUNcluster = stats::kmeans, method = "silhouette") +
    labs(subtitle = "Silhouette method")
  # 3 Gap statistic
  # nboot = 50 to keep the function speedy.
  # recommended value: nboot= 500 for your analysis.
  # Use verbose = FALSE to hide computing progression.
  cp3 <- NULL
  if (!fast) {
    cp3 <- factoextra::fviz_nbclust(x = data_scaled, FUNcluster = stats::kmeans, nstart = 25, method = "gap_stat", nboot = 50) +
      labs(subtitle = "Gap statistic method")
  }
  return(list(cp1 = cp1, cp2 = cp2, cp3 = cp3))
}

#' Filter the highest var and scale
#'
#' @param otu_group grouped otutab
#' @param filter_var filter the highest var
#' @return data_scaled
#' @noRd
filter_top_var <- function(otu_group, filter_var) {
  # trans
  pcutils::dabiao("Filter top ", (1 - filter_var) * 100, "% var and scale")
  group.var <- apply(otu_group, 1, var)
  otu_group.sel <- otu_group[group.var >= stats::quantile(group.var, filter_var), ] # 挑出变化较大的部分
  weight <- c(apply(otu_group.sel, 1, var))
  data_scaled <- pcutils::trans(otu_group.sel, method = "standardize", MARGIN = 1)
  data_scaled
}

#' C-means cluster
#'
#' @param otu_group standardize data
#' @param k_num cluster number
#' @param filter_var filter the highest var
#'
#' @return ggplot
#' @export
#' @family C_means
#' @examples
#' \donttest{
#' if (requireNamespace("e1071") && requireNamespace("factoextra")) {
#'   data(otutab, package = "pcutils")
#'   pcutils::hebing(otutab, metadata$Group) -> otu_group
#'   cm_test_k(otu_group, filter_var = 0.7)
#'   cm_res <- c_means(otu_group, k_num = 3, filter_var = 0.7)
#'   plot(cm_res, 0.8)
#' }
#' }
c_means <- function(otu_group, k_num, filter_var) {
  lib_ps("e1071", library = FALSE)

  data_scaled <- filter_top_var(otu_group, filter_var)

  #-----Start clustering
  # set.seed(123)
  cm <- e1071::cmeans(data_scaled, center = k_num, iter.max = 500)

  cm_data <- cbind.data.frame(
    Name = row.names(data_scaled), data_scaled,
    Weight = apply(otu_group[rownames(data_scaled), ], 1, stats::var),
    Cluster = cm$cluster,
    Membership = apply(cm$membership, 1, max)
  )
  res <- list(
    data = otu_group, filter_var = filter_var, data_scaled = data_scaled,
    cm_data = cm_data, centers = cm$centers, membership = cm$membership
  )
  class(res) <- "cm_res"
  return(res)
}


#' Reporter score analysis after C-means clustering
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples.
#' @param group The comparison groups (at least two categories) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf. And you can use factor levels to change order.
#' @param metadata sample information data.frame contains group
#' @param k_num if NULL, perform the cm_test_k, else an integer
#' @param filter_var see c_means
#' @param verbose verbose
#' @param method method from \code{\link{reporter_score}}
#' @param ... additional arguments for \code{\link{reporter_score}}
#'
#' @return rs_by_cm
#' @export
#' @aliases GRSA_by_cm
#'
#' @family C_means
#' @examples
#' message("The following example require some time to run:")
#' \donttest{
#' if (requireNamespace("e1071") && requireNamespace("factoextra")) {
#'   data("KO_abundance_test")
#'   rsa_cm_res <- RSA_by_cm(KO_abundance, "Group2", metadata,
#'     k_num = 3,
#'     filter_var = 0.7, method = "pearson", perm = 199
#'   )
#'   extract_cluster(rsa_cm_res, cluster = 1)
#' }
#' }
RSA_by_cm <- function(kodf, group, metadata = NULL, k_num = NULL, filter_var = 0.7, verbose = TRUE, method = "pearson", ...) {
  if (verbose) {
    pcutils::dabiao("Checking group")
  }
  if (!is.null(metadata)) {
    if (length(group) != 1) {
      stop("'group' should be one column name of metadata when metadata exsit!")
    }
    idx <- rownames(metadata) %in% colnames(kodf)
    metadata <- metadata[idx, , drop = FALSE]
    kodf <- kodf[, rownames(metadata), drop = FALSE]
    if (verbose) {
      message(nrow(metadata), " samples are matched for next step.")
    }
    if (length(idx) < 2) {
      stop("too less common samples")
    }
    sampFile <- data.frame(group = metadata[, group], row.names = row.names(metadata))
  } else {
    if (length(group) != ncol(kodf)) {
      stop("'group' length should equal to columns number of kodf when metadata is NULL!")
    }
    sampFile <- data.frame(row.names = colnames(kodf), group = group)
  }

  pcutils::hebing(kodf, sampFile$group, 2) -> kodf_g
  if (is.null(k_num)) {
    cm_test_k_res <- cm_test_k(kodf_g, filter_var = filter_var)
    print(cm_test_k_res)
    while (TRUE) {
      if (!is.null(k_num)) {
        message("wrong format, please choose a proper k_num for clustering.")
      } else {
        message("please choose a proper k_num for clustering.")
      }
      k_num <- readline("k_num: ")
      suppressWarnings({
        k_num <- as.numeric(k_num)
      })
      if (is.numeric(k_num) & (!is.na(k_num))) {
        break
      }
    }
  }
  if (verbose) {
    pcutils::dabiao("Choose k_num: ", k_num, ". Start to cluster")
  }

  cm_res <- c_means(kodf_g, k_num = k_num, filter_var = filter_var)

  if (verbose) {
    pcutils::dabiao("Get ReporterScore for each cluster")
  }
  all_res <- list()
  for (i in seq_len(k_num)) {
    if (verbose) {
      pcutils::dabiao("For cluster ", i)
    }
    pattern1 <- cm_res$centers[i, ]
    tmp_res <- reporter_score(kodf = kodf, group = group, metadata = metadata, mode = "directed", method = method, pattern = pattern1, verbose = verbose, ...)
    all_res[[paste0("Cluster", i)]] <- append(tmp_res[2:3], list(pattern = pattern1))
  }
  all_res <- append(all_res, append(tmp_res[c(1, 4:6)], list(cm_res = cm_res)))
  class(all_res) <- "rs_by_cm"
  all_res
}

#' Extract one cluster from rs_by_cm object
#'
#' @rdname RSA_by_cm
#' @param rsa_cm_res rs_by_cm object
#' @param cluster integer
#'
#' @return reporter_score object
#' @export
extract_cluster <- function(rsa_cm_res, cluster = 1) {
  stopifnot(inherits(rsa_cm_res, "rs_by_cm"))
  ncluster <- sum(grepl("Cluster", names(rsa_cm_res)))
  if (!paste0("Cluster", cluster[1]) %in% names(rsa_cm_res)) {
    stop("check the cluster")
  }
  res <- append(rsa_cm_res[[paste0("Cluster", cluster[1])]], rsa_cm_res[c("kodf", "modulelist", "group", "metadata")])
  class(res) <- "reporter_score"
  res
}
