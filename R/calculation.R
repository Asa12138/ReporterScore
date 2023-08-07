#' Print reporter_score
#'
#' @param x reporter_score
#' @param ... add
#'
#' @return No value
#' @exportS3Method
#' @method print reporter_score
#'
print.reporter_score=function(x,...){
    reporter_score_res=x
    pcutils::dabiao("KO abundance table",print=TRUE)
    cat("With ", nrow(reporter_score_res$kodf)," KOs and ",ncol(reporter_score_res$kodf)," samples.\n")
    pcutils::dabiao("group",print=TRUE)
    vs_group=attributes(reporter_score_res$reporter_s)$vs_group

    if(attributes(reporter_score_res$reporter_s)$mode=="directed")title=paste0(vs_group,collapse = "/ ")
    else title=paste0(vs_group,collapse = "/")
    cat("vs group: ",title,"\n",sep = "")
    pcutils::dabiao("parameter",print=TRUE)
    cat("use mode: ",attributes(reporter_score_res$reporter_s)$mode,
        "; use method: ",attributes(reporter_score_res$reporter_s)$method,sep = "")

}

#' One step to get the reporter score of your KO abundance table.
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples.
#' @param group The comparison groups (at least two categories) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf. And you can use factor levels to change order.
#' @param metadata sample information data.frame contains group
#' @param mode "mixed" or "directed" (only for two groups differential analysis or multi-groups correlation analysis.), see details in \code{\link{pvalue2zs}}.
#' @param verbose logical
#' @param method the type of test. Default is `wilcox.test`. Allowed values include:
#' \itemize{
#' \item \code{\link[stats]{t.test}} (parametric) and \code{\link[stats]{wilcox.test}} (non-parametric). Perform comparison between two groups of samples. If the grouping variable contains more than two levels, then a pairwise comparison is performed.
#' \item \code{\link[stats]{anova}} (parametric) and \code{\link[stats]{kruskal.test}} (non-parametric). Perform one-way ANOVA test comparing multiple groups.
#' \item "pearson", "kendall", or "spearman" (correlation), see \code{\link[stats]{cor}}.}
#' @param threads default 1
#' @param p.adjust.method1 p.adjust.method for `ko.test`, see \code{\link[stats]{p.adjust}}
#' @param p.adjust.method2 p.adjust.method for the correction of ReporterScore, see \code{\link[stats]{p.adjust}}
#' @param type "pathway" or "module" for default KOlist_file.
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param perm permutation number, default: 1000.
#'
#' @return reporter_score object：
#' \item{kodf}{your input KO_abundance table}
#' \item{ko_pvalue}{ko statistics result contains p.value}
#' \item{ko_stat}{ko statistics result contains p.value and z_score}
#' \item{reporter_s}{the reporter score in each pathway}
#' \item{modulelist}{default KOlist or customized modulelist dataframe}
#' \item{group}{The comparison groups in your data}
#' \item{metadata}{sample information dataframe contains group}
#' @export
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' reporter_score_res=reporter_score(KO_abundance,"Group",metadata,mode="directed")
#' reporter_score_res2=reporter_score(KO_abundance,"Group2",metadata,mode="mixed",
#'      method = "kruskal.test",p.adjust.method1 = "none")
#' }
reporter_score=function(kodf,group,metadata=NULL,mode=c("mixed","directed")[1],
                        verbose=TRUE,method="wilcox.test",threads=1,
                        p.adjust.method1='BH',p.adjust.method2='BH',
                        type=c("pathway","module")[1],perm =1000,
                        KOlist_file=NULL,modulelist=NULL){
    KOlist=NULL

    stopifnot(mode%in%c("mixed","directed"))

    if(verbose)pcutils::dabiao("1.KO test")
    ko_pvalue=ko.test(kodf,group,metadata,method = method,threads =threads,
                      p.adjust.method =p.adjust.method1,verbose = verbose)
    if(verbose)pcutils::dabiao("2.Transfer p.value to z-score")
    ko_stat=pvalue2zs(ko_pvalue,mode=mode,p.adjust.method =p.adjust.method1)
    if(verbose)pcutils::dabiao("3.Calculating reporter score")
    reporter_s=get_reporter_score(ko_stat,type = type,threads = threads,
                                  p.adjust.method = p.adjust.method2,
                                  KOlist_file =KOlist_file,modulelist = modulelist,
                                  perm = perm,verbose = verbose)
    if(verbose)pcutils::dabiao("All done")

    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment(),verbose=verbose)
        modulelist=KOlist[[type]]
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")

    res=list(kodf=kodf,ko_pvalue=ko_pvalue,ko_stat=ko_stat,reporter_s=reporter_s,modulelist=modulelist,group=group,metadata=metadata)
    class(res)="reporter_score"
    res
}


#' Differential analysis or Correlation analysis for KO-abundance table
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001), colnames is samples
#' @param group The comparison groups (at least two categories) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf. And you can use factor levels to change order.
#' @param metadata sample information data.frame contains group
#' @param method the type of test. Default is `wilcox.test`. Allowed values include:
#' \itemize{
#' \item \code{\link[stats]{t.test}} (parametric) and \code{\link[stats]{wilcox.test}} (non-parametric). Perform comparison between two groups of samples. If the grouping variable contains more than two levels, then a pairwise comparison is performed.
#' \item \code{\link[stats]{anova}} (parametric) and \code{\link[stats]{kruskal.test}} (non-parametric). Perform one-way ANOVA test comparing multiple groups.
#' \item "pearson", "kendall", or "spearman" (correlation), see \code{\link[stats]{cor}}.}
#' @param threads default 1
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#' @param verbose logical
#'
#' @return ko_pvalue dataframe
#' @export
#'
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' ko_pvalue=ko.test(KO_abundance,"Group",metadata)
#' }
ko.test=function(kodf,group,metadata=NULL,method="wilcox.test",threads=1,p.adjust.method='BH',verbose=TRUE){
    i=NULL
    t1 <- Sys.time()

    if(verbose)pcutils::dabiao("Checking rownames")
    rowname_check=grepl("K\\d{5}",rownames(kodf))
    if(!all(rowname_check))message("Some of your kodf are not KO id, check the format! (e.g. K00001)\n")

    if(verbose)pcutils::dabiao("Checking group")
    if(!is.null(metadata)){
        if(length(group)!=1)stop("'group' should be one column name of metadata when metadata exsit!")
        idx = rownames(metadata) %in% colnames(kodf)
        metadata = metadata[idx, , drop = FALSE]
        kodf = kodf[, rownames(metadata),drop=FALSE]
        if(verbose)message(nrow(metadata)," samples are matched for next step.")
        if(length(idx)<2)stop("too less common samples")
        sampFile = data.frame(group=metadata[, group], row.names = row.names(metadata))
    }
    else {
        if(length(group)!=ncol(kodf))stop("'group' length should equal to columns number of kodf when metadata is NULL!")
        sampFile =data.frame(row.names =colnames(kodf),group=group)
    }

    if(verbose)pcutils::dabiao("Removing all-zero rows: ",sum(rowSums(abs(kodf))==0))
    kodf=kodf[rowSums(abs(kodf))>0,]

    vs_group=levels(factor(sampFile$group))

    if(length(vs_group)==1)stop("'group' should be at least two elements factor")
    if(length(vs_group)>2){
        if(method%in%c("t.test","wilcox.test"))stop("'group' more than two elements, try 'kruskal.test' or 'anova'")}

    #calculate each
    if(verbose)pcutils::dabiao("Calculating each KO")
    if(verbose)pcutils::dabiao("Using method: ",method)
    tkodf=t(kodf)%>%as.data.frame()
    group=sampFile$group
    if(method%in%c("pearson", "kendall", "spearman")){
        if(verbose)message("Using correlation analysis: ",method," the groups will be transform to numeric, note the factor level of group.")
        group2=as.numeric(factor(group))
    }

    res.dt=data.frame("KO_id"=rownames(kodf))

    for (i in vs_group) {
        tmpdf=data.frame(average=apply(kodf[,which(group==i)],1,mean),sd=apply(kodf[,which(group==i)],1,sd))
        colnames(tmpdf)=paste0(colnames(tmpdf),"_",i)
        res.dt=cbind(res.dt,tmpdf)
    }
    if(length(vs_group)==2){
        #update, make sure the control group is first one.
        res.dt$diff_mean=res.dt[,paste0("average_",vs_group[2])]-res.dt[,paste0("average_",vs_group[1])]
    }
    if(method%in%c("pearson", "kendall", "spearman")){
        res.dt$cor=cor(tkodf,group2,method = method)[,1]
    }

    high_group <- apply(res.dt[,paste0("average_",vs_group)], 1, function(a) which(a == max(a))[[1]])
    res.dt$Highest=vs_group[high_group]

    #parallel
    reps=nrow(kodf)

    tkodf=t(kodf)
    #main function
    loop=function(i){
        val <- tkodf[,i]
        if(method=="wilcox.test"){
            pval <- stats::wilcox.test(val~group)$p.value
        }
        if(method=="t.test"){
            pval <- stats::t.test(val~group)$p.value
        }
        if(method=="kruskal.test"){
            pval <- stats::kruskal.test(val~group)$p.value
        }
        if(method=="anova"){
            pval <- stats::lm(val~group) %>% stats::anova(.) %>%.$`Pr(>F)` %>% .[1]
        }
        if(method%in%c("pearson", "kendall", "spearman")){
            pval <- stats::cor.test(val,group2,method = method)$p.value
        }
        if(verbose&(i%%1000==0))message(paste(i,"done."))
        pval
    }

    {
        if(threads>1){
            pcutils::lib_ps("foreach","doSNOW","snow")
            pb <- utils::txtProgressBar(max =reps, style = 3)
            opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
            cl <- snow::makeCluster(threads)
            doSNOW::registerDoSNOW(cl)
            res <- foreach::foreach(i = 1:reps,.options.snow = opts
            ) %dopar% {
                suppressWarnings(loop(i))
            }
            snow::stopCluster(cl)
            gc()
            pcutils::del_ps("doSNOW","snow","foreach")
        }
        else {
            res <-suppressWarnings(lapply(1:reps, loop))
        }}
    #simplify method
    res=do.call(c,res)
    res.dt$p.value=res

    t2 <- Sys.time()
    stime <- sprintf("%.3f", t2 - t1)

    res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)

    resinfo <- paste0('Compared groups: ', paste(vs_group,collapse = ', '), "\n",
                      'Total KO number: ', reps, "\n",
                      'Compare method: ', method, "\n",
                      'Time use: ', stime, attr(stime, 'units'), "\n")
    message(resinfo)

    attributes(res.dt)$vs_group=vs_group
    attributes(res.dt)$method=method
    attributes(res.dt)$p.adjust.method=p.adjust.method
    return(res.dt)
}

#' Transfer p-value of KOs to Z-score
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko.test}} or your statistics test data.frame contains column: `p.value`
#' @param mode "mixed" or "directed" (only for two groups differential analysis or multi-groups correlation analysis.), see details
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#'
#' @return ko_stat dataframe
#' @export
#' @details
#' "\strong{mixed}" mode is the original reporter-score method from Patil, K. R. et al. PNAS 2005.
#' In this mode, the reporter score is \strong{undirected}, and the larger the reporter score, the more significant the enrichment, but it cannot indicate the up-and-down regulation information of the pathway！(Liu, L. et al. iMeta 2023.)
#'
#' steps:
#'
#' 1. Use the Wilcoxon rank sum test to obtain the P value of the significance of each KO difference between the two groups (ie \eqn{P_{koi}}, i represents a certain KO);
#'
#' 2. Using an inverse normal distribution, convert the P value of each KO into a Z value (\eqn{Z_{koi}}), the formula:
#'
#' \eqn{Z_{koi}=\theta ^{-1}(1-P_{koi})}
#'
#' 3. "Upgrade" KO to pathway: \eqn{Z_{koi}}, calculate the Z value of the pathway, the formula:
#'
#' \eqn{Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}}
#'
#' where k means A total of k KOs were annotated to the corresponding pathway;
#'
#' 4. Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of \eqn{Z_{pathway}}, the formula:
#'
#' \eqn{Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k}
#'
#' \eqn{μ_k} is The mean of the random distribution, \eqn{σ_k} is the standard deviation of the random distribution.
#'
#' Instead, "\strong{directed}" mode is a derived version of "mixed", referenced from \code{https://github.com/wangpeng407/ReporterScore}.
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
#' 4. "Upgrade" KO to pathway: \eqn{Z_{koi}}, calculate the Z value of the pathway, the formula:
#'
#'  \eqn{Z_{pathway}=\frac{1}{\sqrt{k}}\sum Z_{koi}}
#'
#'  where k means A total of k KOs were annotated to the corresponding pathway;
#'
#' 5. Evaluate the degree of significance: permutation (permutation) 1000 times, get the random distribution of \eqn{Z_{pathway}}, the formula:
#'
#' \eqn{Z_{adjustedpathway}=(Z_{pathway}-\mu _k)/\sigma _k}
#'
#' \eqn{μ_k} is The mean of the random distribution, \eqn{σ_k} is the standard deviation of the random distribution.
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
#'
#' @examples
#' \donttest{
#' data(KO_abundance_test)
#' ko_pvalue=ko.test(KO_abundance,"Group",metadata)
#' ko_stat=pvalue2zs(ko_pvalue,mode="directed")
#' }
pvalue2zs=function(ko_pvalue,mode=c("mixed","directed")[1],p.adjust.method='BH'){
    p.adjust=type=NULL
    res.dt=ko_pvalue
    if(!all(c("p.value")%in%colnames(res.dt))){stop("check if `p.value` in your ko_stat dataframe!")}

    stopifnot(mode%in%c("mixed","directed"))

    if("diff_mean"%in%colnames(res.dt)){
        res.dt$sign <- ifelse(res.dt$diff_mean < 0, -1, 1)
        res.dt$type <- ifelse(res.dt$diff_mean < 0, paste0('Depleted'), paste0('Enriched'))
    }
    if("cor"%in%colnames(res.dt)){
        res.dt$sign <- ifelse(res.dt$cor < 0, -1, 1)
        res.dt$type <- ifelse(res.dt$cor < 0, paste0('Depleted'), paste0('Enriched'))
    }

    #mixed不考虑正负号，p.adjust不除以2，考虑的话除以2
    if(mode=="mixed"){
        res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)
        #逆正态分布
        zs <- stats::qnorm(1-(res.dt$p.adjust))
        #防止太小的p.adjust产生Inf
        zs <- ifelse( zs > 8.209536, 8.209536, zs)
        zs<- ifelse( zs < -8.209536, -8.209536, zs)
        res.dt$Z_score <- zs
        attributes(res.dt)$mode="mixed"
    }
    if(mode=="directed"){
        if(!"type"%in%colnames(res.dt))stop("directed mode only use for two groups differential analysis or multi-groups correlation analysis.")
        pn_sign=2
        res.dt$p.value=res.dt$p.value/pn_sign
        res.dt$p.adjust <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)

        #这种做法可能要基于一个前提，就是上下调ko数量基本一致,才能保证正负都是显著差异的，或者分开正负分析？
        # up_down_ratio=table(res.dt%>%dplyr::filter(abs(p.adjust)<=stats::quantile(res.dt$p.adjust,0.05,na.rm=TRUE))%>%dplyr::pull(type))
        # kafang_res=stats::chisq.test(up_down_ratio)
        # pcutils::dabiao("")
        # pcutils::dabiao("Chi-squared test for up and down ko ratio")
        # message("X-squared = ",round(kafang_res$statistic,4), "   p-value = ",round(kafang_res$p.value,4))
        #
        # #if p-value>0.05，正负一致。
        # if(kafang_res$p.value<0.05){
        #     message("The overall up-down ratio of ko abundance is unbalanced!\n Continuing to use the directed mode may lead to wrong conclusions")
        # }

        #逆正态分布
        zs <- stats::qnorm(1-(res.dt$p.adjust))
        #防止太小的p.adjust产生Inf
        zs <- ifelse( zs > 8.209536, 8.209536, zs)
        zs <- ifelse( zs < -8.209536, -8.209536, zs)

        #通过判断上下调给予z-score正负号，让最后的reporter-score正负号作为上下调标志
        res.dt$Z_score <- ifelse(res.dt$sign < 0, -zs, zs)
        attributes(res.dt)$mode="directed"
    }

    if(is.null(attributes(res.dt)$mode))p_th=0.05
    else p_th=ifelse(attributes(res.dt)$mode=="directed",0.025,0.05)

    if("type"%in%colnames(res.dt)&mode=="directed")res.dt=dplyr::mutate(res.dt,Significantly=ifelse(p.adjust<p_th,type,"None"))
    else res.dt=dplyr::mutate(res.dt,Significantly=ifelse(p.adjust<p_th,"Significant","None"))

    return(res.dt)
}

random_mean_sd <- function(vec, Knum, perm = 1000){
    #set.seed((Knum + 1))
    #我认为应该repalce=TRUE，否则当vec长度小于Knum时，每次取到的都是同一个结果，sd就会非常小！
    #修改Knum为exist_Knum就不会有这个问题了
    #但是Permutation就是不放回抽样😭，Bootstrap才是有放回
    #建议输入的KO表的ko数量多一些，保证sd正确。
    #replace=(length(vec)<=Knum)

    replace=FALSE
    temp=sapply(1:perm, \(i){sum(sample(vec, Knum,replace = replace))/sqrt(Knum)})
    list(vec=temp,mean_sd=c(mean(temp), stats::sd(temp)))
}

#' Calculate reporter score
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}}
#' @param type "pathway" or "module" for default KOlist_file.
#' @param threads threads
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param perm permutation number, default: 1000
#' @param verbose logical
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#'
#' @return reporter_res dataframe
#' @export
#'
#' @examples
#' \donttest{
#' data(KO_abundance_test)
#' ko_pvalue=ko.test(KO_abundance,"Group",metadata)
#' ko_stat=pvalue2zs(ko_pvalue,mode="directed")
#' reporter_s1=get_reporter_score(ko_stat)
#' }
get_reporter_score=function(ko_stat,type=c("pathway","module")[1],threads=1,KOlist_file=NULL,modulelist=NULL,perm =1000,verbose=TRUE,p.adjust.method="BH"){
    KOlist=i=NULL
    type_flag=FALSE
    t1 <- Sys.time()

    if(verbose)pcutils::dabiao("Checking file")
    if(!all(c("KO_id","Z_score")%in%colnames(ko_stat)))stop("Some wrong with ko_stat")
    rowname_check=grepl("K\\d{5}",ko_stat$KO_id)
    if(!all(rowname_check)){if(verbose)message("Some of your ko_stat are not KO id, check the format! (e.g. K00001)!\n")}

    if(is.null(modulelist)){
        type=match.arg(type,c("pathway","module"))
        load_KOlist(KOlist_file,envir = environment(),verbose=verbose)
        modulelist=KOlist[[type]]
        type_flag=TRUE
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")

    #calculate each pathway
    if(verbose)pcutils::dabiao("Calculating each pathway")

    #parallel
    reps=nrow(modulelist)

    if(is.null(attributes(ko_stat)$mode))p_th=0.05
    else p_th=ifelse(attributes(ko_stat)$mode=="directed",0.025,0.05)

    clean.KO <- ko_stat$Z_score[!is.na(ko_stat$Z_score)]
    #main function
    loop=function(i){
        #找到在该pathway里的所有ko的zs
        tmp_kos=strsplit(modulelist$KOs[i], ',')[[1]]

        z <- ko_stat[ko_stat$KO_id %in% tmp_kos,]
        exist_KO=nrow(z)

        significant_KO=sum(z$p.adjust<p_th)

        #如果一条通路里压根没找到几个ko，就不应该有reporterscore
        if((exist_KO<3)&(exist_KO/modulelist$K_num[i]<0.2))return(c(exist_KO,significant_KO,NA,NA,NA,NA,NA))

        #KOnum <- modulelist$K_num[i]
        #KOnum <- ifelse(length(clean.KO) >= KOnum, KOnum, length(clean.KO))

        KOnum=exist_KO
        #以整个输入ko文件作为背景,抽取KOnum应该是exist_KO，而不是所有的KOnum，可以在iMeta文章看到
        mean_sd <- random_mean_sd(clean.KO, KOnum, perm =perm)
        Z_score=sum(z$Z_score) / sqrt(KOnum)

        reporter_score <- (Z_score - mean_sd$mean_sd[1])/mean_sd$mean_sd[2]
        p.value=sum(mean_sd$vec>Z_score)/length(mean_sd$vec)
        p.value=ifelse(p.value>0.5,1-p.value,p.value)
        if(verbose&(i%%100==0))message(paste(i,"done."))
        c(exist_KO,significant_KO,Z_score,mean_sd$mean_sd,reporter_score,p.value)
    }
    {
    if(threads>1){
        pcutils::lib_ps("foreach","doSNOW","snow")
        pb <- utils::txtProgressBar(max =reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
        cl <- snow::makeCluster(threads)
        doSNOW::registerDoSNOW(cl)
        res <- foreach::foreach(i = 1:reps,.options.snow = opts
        ) %dopar% {
            loop(i)
        }
        snow::stopCluster(cl)
        gc()
        pcutils::del_ps("doSNOW","snow","foreach")
    }
    else {
        res <-lapply(1:reps, loop)
    }}
    #simplify method
    res=do.call(rbind,res)
    colnames(res)=c("Exist_K_num","Significant_K_num","Z_score","BG_Mean","BG_Sd","ReporterScore","p.value")
    res=as.data.frame(res)
    reporter_res <- data.frame(ID = modulelist$id,
                               Description = modulelist$Description,
                               K_num=modulelist$K_num,res)
    reporter_res$p.adjust=stats::p.adjust(reporter_res$p.value,p.adjust.method)
    attributes(reporter_res)$mode=attributes(ko_stat)$mode
    attributes(reporter_res)$method=attributes(ko_stat)$method
    attributes(reporter_res)$vs_group=attributes(ko_stat)$vs_group
    if(type_flag)attributes(reporter_res)$type=type
    rownames(reporter_res)=reporter_res$ID

    t2 <- Sys.time()

    stime <- sprintf("%.3f", t2 - t1)
    resinfo <- paste0('ID number: ', reps, "\n",
                      'Time use: ', stime, attr(stime, 'units'), '\n')
    message(resinfo)
    return(reporter_res)
}

#' get KOs in a modulelist
#'
#' @param map_id map_id in modulelist
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}.
#' @param ko_stat NULL or ko_stat result from \code{\link{pvalue2zs}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#'
#' @export
#' @return koids, or dataframe with these koids
get_KOs=function(map_id="map00010",KOlist_file=NULL,ko_stat=NULL,modulelist=NULL){
    KOlist=NULL
    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        if(grepl("map",map_id[1]))modulelist=KOlist$pathway
        if(grepl("M",map_id[1]))modulelist=KOlist$module
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")
    kos=modulelist[which(modulelist$id%in%map_id),"KOs"]
    if(length(kos)>0)kos=lapply(kos, strsplit,",")%>%unlist()
    else kos=NULL
    if(!is.null(ko_stat)){return(ko_stat[ko_stat$KO_id%in%kos,])}
    return(kos)
}


#' Upgrade the KO level
#'
#' @param KO_abundance KO_abundance
#' @param level one of "pathway", "module", "level1", "level2", "level3", "module1", "module2", "module3".
#' @param show_name logical
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#'
#' @return data.frame
#' @export
#'
#' @examples
#' KO_level1=up_level_KO(KO_abundance,level = "level1",show_name = TRUE)
#' pcutils::stackplot(KO_level1[-which(rownames(KO_level1)=="Unknown"),])
up_level_KO=function(KO_abundance,level="pathway",show_name=FALSE,
                     KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    a=KO_abundance
    a$KO_id=rownames(a)

    if(level%in%c("pathway", "module")){
        KOlist=NULL
        if(is.null(modulelist)){
            load_KOlist(KOlist_file,envir = environment(),verbose=verbose)
            modulelist=KOlist[[level]]
        }
        if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")
        path2ko=pcutils::explode(modulelist[,c("id","KOs")],2,split = ",")
        path2name=setNames(modulelist$Description,modulelist$id)
    }
    else if (level%in%c("level1", "level2", "level3")){
        KO_htable=NULL
        load_KO_htable(envir = environment(),verbose=verbose)
        path2ko=KO_htable[,c(paste0(level,"_id"),"KO_id")]

        path2name=KO_htable[,paste0(level,c("_id","_name"))]%>%dplyr::distinct()
        path2name=setNames(path2name[,paste0(level,"_name"),drop=TRUE],
                           path2name[,paste0(level,"_id"),drop=TRUE])
    }
    else if (level%in%c("module1", "module2", "module3")){
        Module_abundance=up_level_KO(KO_abundance,level = "module")
        a=Module_abundance
        a$KO_id=rownames(a)
        load_Module_htable(envir = environment(),verbose=verbose)
        path2ko=Module_htable[,c(paste0(level,"_name"),"Module_id")]
        show_name=FALSE
    }
    else stop('level should be one of "pathway", "module", "level1", "level2", "level3", "module1", "module2", "module3".')

    colnames(path2ko)=c("Pathway","KO_id")
    path2ko=dplyr::distinct_all(path2ko)
    dplyr::left_join(a,path2ko,by="KO_id")->aa
    aa$Pathway[is.na(aa$Pathway)]="Unknown"

    b=pcutils::hebing(dplyr::select(aa,-c("KO_id","Pathway")),aa$Pathway,1,act = "sum")
    if(show_name)rownames(b)=c(path2name,"Unknown"="Unknown")[rownames(b)]
    b
}
