
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
    vs_group=grep("avg",colnames(reporter_score_res$ko_stat),value = TRUE)
    cat("vs group: ",sub("avg_","",vs_group[1])," vs ",sub("avg_","",vs_group[2]),
        "; use mode: ",attributes(reporter_score_res$reporter_s)$mode,sep = "")
}


#' One step to get the reporter score of your KO abundance table.
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples.
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf. And you can use factor levels to change order.
#' @param metadata sample information dataframe contains group.
#' @param mode "mixed" or "directed", see details in \code{\link{pvalue2zs}}.
#' @param verbose logical
#' @param threads default 1
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#' @param type "pathway" or "module"
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
#' \item{group}{The compare group (two category) in your data}
#' \item{metadata}{sample information dataframe contains group}
#' @export
#' @examples
#' \donttest{
#' reporter_score_res=reporter_score(KO_abundance,"Group",metadata,mode="directed")
#' }
reporter_score=function(kodf,group,metadata=NULL,mode=c("mixed","directed")[1],
                        verbose=TRUE,threads=1,p.adjust.method='BH',
                        type=c("pathway","module")[1],perm =1000,
                        KOlist_file=NULL,modulelist=NULL){
    KOlist=NULL

    stopifnot(mode%in%c("mixed","directed"))
    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        modulelist=KOlist[[type]]
        if(verbose){
            pcutils::dabiao("load KOlist")
            if(!is.null(attributes(KOlist)$"download_time")){
                pcutils::dabiao(paste0("KOlist download time: ",attributes(KOlist)$"download_time"))
                message("If you want to update KOlist, use `update_KO_file()`")
            }
        }
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")


    if(verbose)pcutils::dabiao("1.KO test")
    ko_pvalue=ko.test(kodf,group,metadata,threads =threads,p.adjust.method =p.adjust.method,verbose = verbose)
    if(verbose)pcutils::dabiao("2.Transfer p.value to z-score")
    ko_stat=pvalue2zs(ko_pvalue,mode=mode,p.adjust.method =p.adjust.method)
    if(verbose)pcutils::dabiao("3.Calculating reporter score")
    reporter_s=get_reporter_score(ko_stat,type = type,threads = threads,
                                  KOlist_file =KOlist_file,modulelist = modulelist,perm = perm,verbose = verbose)
    if(verbose)pcutils::dabiao("All done")

    res=list(kodf=kodf,ko_pvalue=ko_pvalue,ko_stat=ko_stat,reporter_s=reporter_s,modulelist=modulelist,group=group,metadata=metadata)
    class(res)="reporter_score"
    res
}


#' Wilcox-test or t.test for KO-abundance table
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001), colnames is samples
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata sample information dataframe contains group
#' @param threads default 1
#' @param verbose logical
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#'
#' @return ko_pvalue dataframe
#' @export
#'
#' @examples
#' \donttest{
#' data("KO_abundance_test")
#' ko_pvalue=ko.test(KO_abundance,"Group",metadata)
#' }
ko.test=function(kodf,group,metadata=NULL,threads=1,p.adjust.method='BH',verbose=TRUE){
    i=NULL
    t1 <- Sys.time()

    if(verbose)pcutils::dabiao("Checking rownames")
    rowname_check=grepl("K\\d{5}",rownames(kodf))
    if(!all(rowname_check))warning("Some of your kodf are not KO id, check the format! (e.g. K00001)\n")

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

    kodf=kodf[rowSums(kodf)>0,]

    if(!nlevels(factor(sampFile$group))==2)stop("'group' should be two elements factor")
    {vs_group=levels(factor(sampFile$group))}
    #calculate each
    if(verbose)pcutils::dabiao("Calculating each KO")

    #parallel
    reps=nrow(kodf)

    tkodf=t(kodf)
    g1=(sampFile$group==vs_group[1])
    g2=(sampFile$group==vs_group[2])
    exact <- (sum(g1) < 50) && (sum(g2) < 50)
    #main function
    loop=function(i){
        val1 <- tkodf[g1,i]
        val2 <- tkodf[g2,i]

        r <- rank(c(val1, val2))
        TIES <- (length(r) != length(unique(r)))

        if((!exact)|TIES){
            pval <- stats::t.test(val1, val2)$p.value
            #resinfo <- paste0(s, ": Ties exists or exact is false in ", kn, ", using t.test insead!")
        }else{
            pval <- stats::wilcox.test(val1, val2)$p.value
        }
        if(verbose&(i%%100==0))message(paste(i,"done."))
        data.frame(colnames(tkodf)[i],mean(val1), stats::sd(val1), mean(val2), stats::sd(val2), mean(val1) - mean(val2), pval)
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
    res.dt <- data.frame(res, stringsAsFactors = FALSE)

    colnames(res.dt) <- c('KO_id',
                          paste0('avg_', vs_group[1]),
                          paste0('sd_', vs_group[1]),
                          paste0('avg_', vs_group[2]),
                          paste0('sd_', vs_group[2]),
                          'diff_mean', 'p.value')
    t2 <- Sys.time()
    stime <- sprintf("%.3f", t2 - t1)
    resinfo <- paste0('Compared groups: ', vs_group[1], ' and ', vs_group[2], "\n",
                      'Total KO number: ', reps, "\n",
                      'Time use: ', stime, attr(stime, 'units'), "\n")

    message(resinfo)
    res.dt$q.value <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)
    return(res.dt)
}

#' Transfer p-value of KOs to Z-score
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko.test}}
#' @param mode "mixed" or "directed", see details
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#'
#' @return ko_stat dataframe
#' @export
#' @details
#' "\strong{mixed}" mode is the original reporter-score method from Patil, K. R. et al. PNAS 2005.
#' In this mode, the reporter score is \strong{Undirected}, and the larger the reporter score, the more significant the enrichment, but it cannot indicate the up-and-down regulation information of the pathway！(Liu, L. et al. iMeta 2023.)
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
    q.value=type=NULL
    res.dt=ko_pvalue
    if(!all(c("p.value")%in%colnames(res.dt))){stop("check if `p.value` in your ko_stat dataframe!")}
    stopifnot(mode%in%c("mixed","directed"))

    if("diff_mean"%in%colnames(res.dt)){
        res.dt$sign <- ifelse(res.dt$diff_mean < 0, -1, 1)
        res.dt$type <- ifelse(res.dt$diff_mean < 0, paste0('Depleted'), paste0('Enriched'))
    }

    #mixed不考虑正负号，q.value不除以2，考虑的话除以2
    if(mode=="mixed"){
        res.dt$q.value <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)
        #逆正态分布
        zs <- stats::qnorm(1-(res.dt$q.value))
        res.dt$Z_score <- ifelse( zs < -8.209536, -8.209536, zs)
        attributes(res.dt)$mode="mixed"
    }
    if(mode=="directed"){
        if(!"diff_mean"%in%colnames(res.dt))stop("directed mode only use for two group and get the `diff_mean`.")
        pn_sign=2
        res.dt$p.value=res.dt$p.value/pn_sign
        res.dt$q.value <- stats::p.adjust(res.dt$p.value, method = p.adjust.method)

        #这种做法可能要基于一个前提，就是上下调ko数量基本一致,才能保证正负都是显著差异的，或者分开正负分析？
        up_down_ratio=table(res.dt%>%dplyr::filter(abs(q.value)<=stats::quantile(res.dt$q.value,0.05,na.rm=TRUE))%>%dplyr::pull(type))
        kafang_res=stats::chisq.test(up_down_ratio)
        pcutils::dabiao("")
        pcutils::dabiao("Chi-squared test for up and down ko ratio")
        message("X-squared = ",round(kafang_res$statistic,4), "   p-value = ",round(kafang_res$p.value,4))
        #if p-value>0.05，正负一致。
        if(kafang_res$p.value<0.05){
            warning("The overall up-down ratio of ko abundance is unbalanced!\n Continuing to use the directed mode may lead to wrong conclusions")
        }
        #逆正态分布
        zs <- stats::qnorm(1-(res.dt$q.value))
        res.dt$Z_score <- ifelse( zs < -8, -8, zs)
        #通过判断上下调给予z-score正负号，让最后的reporter-score正负号作为上下调标志
        res.dt$Z_score <- ifelse(res.dt$sign < 0, -res.dt$Z_score, res.dt$Z_score)
        attributes(res.dt)$mode="directed"
    }
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
#' @param type "pathway" or "module"
#' @param threads threads
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param perm permutation number, default:1000
#' @param verbose logical
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
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
get_reporter_score=function(ko_stat,type=c("pathway","module")[1],threads=1,KOlist_file=NULL,modulelist=NULL,perm =1000,verbose=TRUE){
    KOlist=i=NULL
    type=match.arg(type,c("pathway","module"))
    t1 <- Sys.time()

    if(verbose)pcutils::dabiao("Checking file")
    if(!all(c("KO_id","Z_score")%in%colnames(ko_stat)))stop("Some wrong with ko_stat")
    rowname_check=grepl("K\\d{5}",ko_stat$KO_id)
    if(!all(rowname_check)){if(verbose)warning("Some of your ko_stat are not KO id, check the format! (e.g. K00001)!\n")}

    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        modulelist=KOlist[[type]]
        if(verbose){
            pcutils::dabiao("load KOlist")
            if(!is.null(attributes(KOlist)$"download_time")){
                pcutils::dabiao(paste0("KOlist download time: ",attributes(KOlist)$"download_time"))
                message("If you want to update KOlist, use `update_KO_file()`")
            }
        }
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")

    #calculate each pathway
    if(verbose)pcutils::dabiao("Calculating each pathway")

    #parallel
    reps=nrow(modulelist)

    p_th=ifelse(attributes(ko_stat)$mode=="mixed",0.05,0.025)
    clean.KO <- ko_stat$Z_score[!is.na(ko_stat$Z_score)]
    #main function
    loop=function(i){
        #找到在该pathway里的所有ko的zs
        tmp_kos=strsplit(modulelist$KOs[i], ',')[[1]]

        z <- ko_stat[ko_stat$KO_id %in% tmp_kos,]
        exist_KO=nrow(z)

        significant_KO=sum(z$q.value<p_th)

        #如果一条通路里压根没找到几个ko，就不应该有reporterscore
        if(exist_KO<3)return(c(exist_KO,significant_KO,NA,NA,NA,NA,NA))

        #KOnum <- modulelist$K_num[i]
        #KOnum <- ifelse(length(clean.KO) >= KOnum, KOnum, length(clean.KO))

        KOnum=exist_KO
        #以整个输入ko文件作为背景,抽取KOnum应该是exist_KO，而不是所有的KOnum，可以在iMeta文章看到
        mean_sd <- random_mean_sd(clean.KO, KOnum, perm =perm)
        Z_score=sum(z$Z_score) / sqrt(KOnum)

        reporter_score <- (Z_score - mean_sd$mean_sd[1])/mean_sd$mean_sd[2]
        p.value=sum(mean_sd$vec>Z_score)/length(mean_sd$vec)
        p.value=ifelse(p.value>0.5,1-p.value,p.value)
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

    attributes(reporter_res)$mode=attributes(ko_stat)$mode
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
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param ko_stat NULL or ko_stat result from \code{\link{pvalue2zs}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#'
#' @export
#' @return koids, or dataframe with these koids
get_KOs=function(map_id="map00010",KOlist_file=NULL,ko_stat=NULL,modulelist=NULL){
    KOlist=NULL
    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        if(grepl("map",map_id))modulelist=KOlist$pathway
        if(grepl("M",map_id))modulelist=KOlist$module
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")
    kos=modulelist[which(modulelist$id==map_id),"KOs"]
    kos=strsplit(kos, ',')[[1]]
    if(!is.null(ko_stat)){return(ko_stat[ko_stat$KO_id%in%kos,])}
    return(kos)
}

#' Plot the reporter_res
#'
#' @param reporter_res result of `get_reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param y_text_size y_text_size
#' @param str_width str_width to wrap
#' @param mode 1～2 plot style.
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_report(reporter_score_res,rs_threshold=c(2,-2),y_text_size=10,str_width=40)
plot_report<-function(reporter_res,rs_threshold=1.64,mode=1,y_text_size=13,str_width=50){
    if(inherits(reporter_res,"reporter_score"))reporter_res=reporter_res$reporter_s
    reporter_res=na.omit(reporter_res)
    Group=Description=ReporterScore=Exist_K_num=NULL
    if(length(rs_threshold)==1)rs_threshold=c(rs_threshold,-rs_threshold)
    rs_threshold=sort(rs_threshold)
    if(rs_threshold[2]>max((reporter_res$ReporterScore))){
        rs_threshold[2]=tail(sort((reporter_res$ReporterScore)))[1]%>%round(.,4)
        warning("Too big rs_threshold, change rs_threshold to ", rs_threshold[1])
    }

    if(attributes(reporter_res)$mode=="directed"){
        if(rs_threshold[1]<min((reporter_res$ReporterScore))){
            rs_threshold[1]=head(sort((reporter_res$ReporterScore)))[5]%>%round(.,4)
            warning("Too small rs_threshold, change rs_threshold to", rs_threshold[1])
        }

        reporter_res2 <- reporter_res[(reporter_res$ReporterScore >= rs_threshold[2])|(reporter_res$ReporterScore <= rs_threshold[1]), ]
    }
    else reporter_res2 <- reporter_res[reporter_res$ReporterScore >= rs_threshold[2],]

    reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0, 'P', 'N')
    reporter_res2 <- reporter_res2[stats::complete.cases(reporter_res2), ]

    if(mode==1){
        p=ggplot(reporter_res2, aes(ReporterScore,stats::reorder(Description, ReporterScore), fill = Group)) +
            geom_bar(stat = 'identity', position='dodge')+
            scale_fill_manual(values=c('P'='orange','N'='seagreen'))+
            theme_light()+
            theme(legend.position = "none")
    }
    if(mode==2){
        p=ggplot(reporter_res2, aes(ReporterScore,stats::reorder(Description, ReporterScore),
                                  size=Exist_K_num, fill =Exist_K_num)) +
            geom_point(shape=21)+
            scale_fill_gradient(low = "#FF000033",high = "red",guide = "legend")+theme_light()
    }

    p <-p+labs(y="")+
        geom_vline(xintercept = rs_threshold[2], linetype =2)+
        scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
        theme(
            axis.text.x = element_text(colour='black',size=13),
            axis.text.y = element_text(colour='black',size=y_text_size)
        )

    if(attributes(reporter_res)$mode=="directed")p=p+geom_vline(xintercept = rs_threshold[1], linetype = 2)
    return(p)
}


#' Plot KOs trend in one pathway or module
#'
#' @param map_id the pathway or module id
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}}
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param box_color box and point color, default: c("#e31a1c","#1f78b4")
#' @param line_color line color, default: c("Depleted"="seagreen","Enriched"="orange","None"="grey")
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_KOs_in_pathway(ko_stat = reporter_score_res,map_id="map00780")
plot_KOs_in_pathway=function(ko_stat,map_id="map00780",
                             KOlist_file=NULL,modulelist=NULL,
                             box_color=c("#e31a1c","#1f78b4"),
                             line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey")){
    Group=value=Group2=value2=Group1=value1=type=Significantly=KOlist=q.value=NULL
    pcutils::lib_ps("ggnewscale","reshape2",library = FALSE)

    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat
        modulelist=reporter_res$modulelist
    }

    A=get_KOs(map_id =map_id ,ko_stat = ko_stat)

    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        if(grepl("map",map_id))modulelist=KOlist$pathway
        if(grepl("M",map_id))modulelist=KOlist$module
    }
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")
    Description=modulelist[modulelist$id==map_id,"Description"]

    p_th=ifelse(attributes(ko_stat)$mode=="mixed",0.05,0.025)
    A=dplyr::mutate(A,Significantly=ifelse(q.value<p_th,type,"None"))

    vs_group=grep("avg",colnames(A),value = TRUE)
    box_df=reshape2::melt(A[,c("KO_id",vs_group)],id.vars="KO_id",variable.name ="Group")
    box_df$Group=factor(box_df$Group,levels = rev(vs_group))

    line_df=A[,c("KO_id",vs_group,"Significantly")]
    colnames(line_df)=c("KO_id","value1","value2","Significantly")
    line_df$Group1=vs_group[1];    line_df$Group2=vs_group[2]

    ggplot()+
        geom_boxplot(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        geom_point(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        scale_color_manual(values = box_color)+
        labs(title = paste0("KOs in ",map_id," (",Description,")"),x=NULL,y="Abundance")+
        ggnewscale::new_scale_color()+
        geom_segment(data = line_df,
                     aes(x=Group2,y=value2,xend=Group1,yend=value1,color=Significantly))+
        scale_color_manual(values = line_color)+
        theme_classic(base_size = 13)+theme(axis.text = element_text(color = "black"))
}

#' Plot KOs boxplot
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param box_param parameters pass to \code{\link[pcutils]{group_box}}
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_box(reporter_score_res,"Group",metadata,
#'      select_ko=c("K00059","K00208","K00647","K00652","K00833","K01012"))
#'
plot_KOs_box=function(kodf,group=NULL,metadata=NULL,
                      map_id="map00780",select_ko=NULL,box_param=NULL,
                      KOlist_file = NULL,modulelist = NULL){
    if(inherits(kodf,"reporter_score")){
        reporter_res=kodf
        kodf=reporter_res$kodf
        group=reporter_res$group
        metadata=reporter_res$metadata
        modulelist=reporter_res$modulelist
    }

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,KOlist_file =KOlist_file,modulelist =modulelist )

    tkodf=kodf[]%>%t()%>%as.data.frame()
    cols=which(colnames(tkodf)%in%select_ko)

    if(length(cols)==0)stop("No select KOs! check map_id or select_ko")
    if(length(cols)>36){
        print(("Too many KOs, do you still want to plot?"))
        flag=readline("yes/no(y/n)?")
        if(!tolower(flag)%in%c("yes","y"))return(NULL)
    }

    metadata[,group]=factor(metadata[,group],levels = rev(levels(factor(metadata[,group]))))

    do.call(pcutils::group_box,
            append(list(tab = tkodf[,cols],group = group,metadata = metadata),
                   pcutils::update_param(list(p_value1 = TRUE,trend_line = TRUE),box_param)))+
        theme_classic(base_size = 13)+theme(axis.text = element_text(color = "black"))+
        scale_fill_manual(values = c("#e31a1c","#1f78b4"))+scale_color_manual(values = c("#e31a1c","#1f78b4"))
}

#' Plot KOs heatmap
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param heatmap_param parameters pass to \code{\link[pheatmap]{pheatmap}}
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_heatmap(reporter_score_res,"Group",metadata,map_id="map00780")
#'
plot_KOs_heatmap=function(kodf,group=NULL,metadata=NULL,
                          map_id="map00780",select_ko=NULL,
                          KOlist_file = NULL,modulelist = NULL,heatmap_param=list()){
    pcutils::lib_ps("pheatmap",library = FALSE)
    if(inherits(kodf,"reporter_score")){
        reporter_res=kodf
        kodf=reporter_res$kodf
        group=reporter_res$group
        metadata=reporter_res$metadata
        modulelist=reporter_res$modulelist
    }

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,KOlist_file =KOlist_file,modulelist =modulelist)
    cols=which(rownames(kodf)%in%select_ko)

    metadata[,group]=factor(metadata[,group],levels = rev(levels(factor(metadata[,group]))))
    if(length(cols)==0)stop("No select KOs! check map_id or select_ko")

    do.call(pheatmap::pheatmap,
            pcutils::update_param(list(mat = kodf[cols,],cluster_cols = F,annotation_col = metadata,scale = "row"),
                                  heatmap_param))
}
