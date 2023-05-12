#' Wilcox-test for KO-abundance table,
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples
#' @param group one column name of metadata when metadata exist or a vector length equal to columns number of kodf
#' @param metadata sample information dataframe contains group
#' @param vs_group vs_group should contains group levels
#' @param verbose logical
#' @param threads default 1
#'
#' @return ko_pvalue dataframe
#' @export
#'
#' @examples
#' data(KO_test)
#' ko_pvalue=ko_test(KO_abundance,"Group",Group_tab)
ko_test=function(kodf,group,metadata=NULL,vs_group=NULL,verbose=T,threads=1){
    t1 <- Sys.time()
    if(verbose)dabiao("Checking rownames")
    rowname_check=grepl("K\\d{5}",rownames(kodf))
    if(!all(rowname_check))warning("Some of your kodf are not KO id, check the format! (e.g. K00001)")

    if(verbose)dabiao("Checking group")
    if(!is.null(metadata)){
        if(length(group)!=1)stop("'group' should be one column name of metadata when metadata exsit!")
        idx = rownames(metadata) %in% colnames(kodf)
        metadata = metadata[idx, , drop = F]
        kodf = kodf[, rownames(metadata),drop=F]
        message(nrow(metadata)," samples are matched for next step.")
        sampFile = data.frame(group=metadata[, group], row.names = row.names(metadata))
    }
    else {
        if(length(group)!=ncol(kodf))stop("'group' length should equal to columns number of kodf!")
        sampFile =data.frame(row.names =colnames(kodf),group=group)
    }

    if(!nlevels(factor(sampFile$group))==2)stop()
    if(!is.null(vs_group)&(!setequal(levels(factor(sampFile$group)),vs_group)))stop("vs_group should contains group levels")
    if(is.null(vs_group)){vs_group=levels(factor(sampFile$group))}
    #calculate each
    if(verbose)dabiao("Calculating each KO")

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
            pval <- t.test(val1, val2)$p.value
            #resinfo <- paste0(s, ": Ties exists or exact is false in ", kn, ", using t.test insead!")
            #cat(resinfo, "\n")
        }else{
            pval <- wilcox.test(val1, val2)$p.value
        }
        if(verbose&(i%%100==0))print(paste(i,"done."))
        data.frame(colnames(tkodf)[i],mean(val1), sd(val1), mean(val2), sd(val2), mean(val1) - mean(val2), pval)
    }
    {
    if(threads>1){
        lib_ps("foreach","doSNOW")
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
    }
    else {
        res <-lapply(1:reps, loop)
    }}
    #simplify method
    res=do.call(rbind,res)
    res.dt <- data.frame(res, stringsAsFactors = F)

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

    cat(resinfo)
    return(res.dt)
}

#' Transfer p-value of KOs to Z-score
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko_test}}
#' @param mode "mixed" or "directed", see details
#' @param p.adjust.method p.adjust.method, see \code{\link[stats]{p.adjust}}
#'
#' @return ko_stat dataframe
#' @export
#' @details
#' "\strong{mixed}" mode is the original reporter-score method from Patil, K. R. et al. PNAS 2005.
#' In this mode, the reporter score is \strong{non-directional}, and the larger the reporter score, the more significant the enrichment, but it cannot indicate the up-and-down regulation information of the pathway！(Liu, L. et al. iMeta 2023.)
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
#' Instead, "\strong{directed}" mode is a derived version of "mixed", referenced from \link{https://github.com/wangpeng407/ReporterScore}.
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
#' In this mode, the Reporter score is directional, and a larger positive value represents a significant up-regulation enrichment, and a smaller Negative values represent significant down-regulation enrichment.
#'
#' However, the disadvantage of this mode is that when a pathway contains about the same number of significantly up-regulates KOs and significantly down-regulates KOs, the final absolute value of Reporter score may approach 0, becoming a pathway that has not been significantly enriched.
#'
#'
#' @references
#' 1. Patil, K. R. & Nielsen, J. Uncovering transcriptional regulation of metabolism by using metabolic network topology. Proc Natl Acad Sci U S A 102, 2685–2689 (2005).
#' 2. Liu, L., Zhu, R. & Wu, D. Misuse of reporter score in microbial enrichment analysis. iMeta n/a, e95.
#'
#'
#' @examples
#' data(KO_test)
#' ko_pvalue=ko_test(KO_abundance,"Group",Group_tab)
#' ko_stat=pvalue2zs(ko_pvalue,mode="directed")
pvalue2zs=function(ko_pvalue,mode=c("mixed","directed")[1],p.adjust.method='BH'){
    res.dt=ko_pvalue
    if(!all(c("p.value")%in%colnames(res.dt))){stop("check if p.value in your ko_stat dataframe!")}

    if("diff_mean"%in%colnames(res.dt)){
        res.dt$sign <- ifelse(res.dt$diff_mean < 0, -1, 1)
        res.dt$type <- ifelse(res.dt$diff_mean < 0, paste0('Depleted'), paste0('Enriched'))
    }

    #mixed不考虑正负号，q.value不除以2，考虑的话除以2
    if(mode=="mixed"){
        res.dt$q.value <- p.adjust(res.dt$p.value, method = p.adjust.method)
        #逆正态分布
        zs <- qnorm(1-(res.dt$q.value))
        res.dt$Z_score <- ifelse( zs < -8.209536, -8.209536, zs)
        attributes(res.dt)$mode="mixed"
    }
    if(mode=="directed"){
        if(!"diff_mean"%in%colnames(res.dt))stop("directed mode only use for two group and get the diff_mean.")
        pn_sign=2
        res.dt$p.value=res.dt$p.value/pn_sign
        res.dt$q.value <- p.adjust(res.dt$p.value, method = p.adjust.method)

        #这种做法可能要基于一个前提，就是上下调ko数量基本一致,才能保证正负都是显著差异的，或者分开正负分析？
        up_down_ratio=table(res.dt%>%filter(abs(q.value)<=quantile(res.dt$q.value,0.05))%>%pull(type))
        kafang_res=chisq.test(up_down_ratio)
        dabiao("")
        print(kafang_res)
        #if p-value>0.05，正负一致。
        if(kafang_res$p.value<0.05){
            warning("The overall up-down ratio of ko abundance is unbalanced!\n Continuing to use the directed mode may lead to wrong conclusions")
        }
        #逆正态分布
        zs <- qnorm(1-(res.dt$q.value))
        res.dt$Z_score <- ifelse( zs < -8.209536, -8.209536, zs)
        #通过判断上下调给予z-score正负号，让最后的reporter-score正负号作为上下调标志
        res.dt$Z_score <- ifelse(res.dt$sign < 0, -res.dt$Z_score, res.dt$Z_score)
        attributes(res.dt)$mode="directed"
    }
    return(res.dt)
}

#' Calculate reporter score
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}}
#' @param verbose logical
#' @param mode "pathway" or "module"
#' @param threads threads
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#'
#' @return reporter_res dataframe
#' @export
#'
#' @examples
#' data(KO_test)
#' ko_pvalue=ko_test(KO_abundance,"Group",Group_tab)
#' ko_stat=pvalue2zs(ko_pvalue,mode="directed")
#' reporter_s=get_reporter_score(ko_stat)
get_reporter_score=function(ko_stat,mode=c("pathway","module")[1],verbose=T,threads=1,KOlist_file=NULL){
    mode=match.arg(mode,c("pathway","module"))
    t1 <- Sys.time()
    if(verbose)dabiao("Checking file")
    if(!all(c("KO_id","")%in%colnames(ko_stat)))
    rowname_check=grepl("K\\d{5}",ko_stat$KO_id)
    if(!all(rowname_check))warning("Some of your ko_stat are not KO id, check the format! (e.g. K00001)!")

    if(is.null(KOlist_file)){
        KOlist_file=system.file("data","new_KOlist.rda",package = "ReporterScore")
        if(!file.exists(KOlist_file))KOlist_file=system.file("data","KOlist.rda",package = "ReporterScore")
    }
    if(file.exists(KOlist_file))load(KOlist_file)

    if(verbose){
        dabiao("load KOlist")
        if(!is.null(attributes(KOlist)$"download_time")){
            dabiao(paste0("KOlist download time: ",attributes(KOlist)$"download_time"))
            message("If you want to update KOlist, use `update_KO_file()`")
        }
    }

    modulelist=KOlist[[mode]]
    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist format!")

    random_mean_sd <- function(vec, Knum, perm = 1000){
        set.seed((Knum + 1))
        #我认为应该repalce=T，否则当vec长度小于Knum时，每次取到的都是同一个结果，sd就会非常小！
        #但是Permutation就是不放回抽样😭，Bootstrap才是有放回
        #建议输入的KO表数量多一些，保证sd正确。
        replace=(length(vec)<=Knum)
        temp=sapply(1:perm, \(i){sum(sample(vec, Knum,replace = replace))/sqrt(Knum)})
        c(mean(temp), sd(temp))
    }

    #calculate each pathway
    if(verbose)dabiao("Calculating each pathway")

    #parallel
    reps=nrow(modulelist)

    #main function
    loop=function(i){
        #找到在该pathway里的所有ko的zs
        z <- ko_stat$Z_score[ko_stat$KO_id %in% strsplit(modulelist$KOs[i], ',')[[1]]]
        KOnum <- modulelist$K_num[i]
        clean.KO <- ko_stat$Z_score[!is.na(ko_stat$Z_score)]
        KOnum <- ifelse(length(clean.KO) >= KOnum, KOnum, length(clean.KO))

        #以整个输入ko文件作为背景
        mean_sd <- random_mean_sd(clean.KO, KOnum, 1000)

        reporter_score <- (sum(z) / sqrt(KOnum) - mean_sd[1])/mean_sd[2]
        reporter_score
    }
    {
    if(threads>1){
        lib_ps("foreach","doSNOW")
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
    }
    else {
        res <-lapply(1:reps, loop)
    }}
    #simplify method
    res=do.call(c,res)

    reporter_res <- data.frame(ID = modulelist$id,
                               ReporterScore = res,
                               Description = modulelist$Description,
                               K_num=modulelist$K_num)

    attributes(reporter_res)$mode=attributes(ko_stat)$mode
    t2 <- Sys.time()
    stime <- sprintf("%.3f", t2 - t1)
    resinfo <- paste0('ID number: ', reps, "\n",
                      'Time use: ', stime, attr(stime, 'units'), '\n')
    cat(resinfo)
    return(reporter_res)
}


#' Plot the reporter_res
#'
#' @param reporter_res result of `get_reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param color color vector c('#e31a1c','#47B0D9')
#' @param y_text_size y_text_size
#' @param str_width str_width to wrap
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(KO_test)
#' ko_pvalue=ko_test(KO_abundance,"Group",Group_tab)
#' ko_stat=pvalue2zs(ko_pvalue,mode="directed")
#' reporter_s=get_reporter_score(ko_stat)
#' plot_report(reporter_s,rs_threshold=c(2,-7),y_text_size=10,str_width=40)
plot_report<-function(reporter_res,rs_threshold=1.64,mode=1,y_text_size=13,str_width=50){
    if(length(rs_threshold)==1)rs_threshold=c(rs_threshold,-rs_threshold)
    if(rs_threshold[1]>max((reporter_res$ReporterScore))){
        rs_threshold[1]=tail(sort((reporter_res$ReporterScore)))[1]
        warning("Too big rs_threshold, change rs_threshold to", rs_threshold)
    }
    if(rs_threshold[1]<min((reporter_res$ReporterScore))){
        rs_threshold[1]=head(sort((reporter_res$ReporterScore)))[1]
        warning("Too small rs_threshold, change rs_threshold to", rs_threshold)
    }

    if(attributes(reporter_res)$mode=="directed"){
        reporter_res2 <- reporter_res[(reporter_res$ReporterScore >= rs_threshold[1])|(reporter_res$ReporterScore <= rs_threshold[2]), ]
    }
    else reporter_res2 <- reporter_res[reporter_res$ReporterScore >= rs_threshold, ]

    reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0, 'P', 'N')
    reporter_res2 <- reporter_res2[complete.cases(reporter_res2), ]

    if(mode==1){
        p=ggplot(reporter_res2, aes(reorder(Description, ReporterScore), ReporterScore, fill = Group)) +
            geom_bar(stat = 'identity', position='dodge')+
            scale_fill_manual(values=c('P'='#e31a1c','N'='#47B0D9'))
    }
    if(mode==2){
        p=ggplot(reporter_res2, aes(reorder(Description, ReporterScore),
                                  ReporterScore,size=K_num, fill =K_num)) +
            geom_point(shape=21)+
            scale_fill_gradient(low = "#FF000033",high = "red",guide = "legend")
    }

    p <-p+
        geom_hline(yintercept = rs_threshold[1], linetype =2)+
        coord_flip()+
        scale_x_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
        theme_light()+
        theme(
            legend.position = "none",
            axis.title.y=element_blank(),
            axis.text.x = element_text(colour='black',size=13),
            axis.text.y = element_text(colour='black',size=y_text_size)
        )

    if(attributes(reporter_res)$mode=="directed")p=p+geom_hline(yintercept = rs_threshold[2], linetype = 2)
    return(p)
}



#' Update the KO2Pathway and KO2Module files
#'
#' @param pack_dir the ReporterScore package location, detect automatically.
#'
#' @export
#' @description
#' Download links:
#'
#' \link{https://rest.kegg.jp/link/pathway/ko}
#' \link{https://rest.kegg.jp/list/pathway}
#' \link{https://rest.kegg.jp/link/module/ko}
#' \link{https://rest.kegg.jp/list/module}
#'
#' you can also download yourself the use \code{\link{make_KO_list}} to get a KOlist object
#' @examples
#' update_KO_file()
#'
update_KO_file=function(pack_dir=NULL){
    if(is.null(pack_dir))pack_dir=paste0(.libPaths()[1],"/ReporterScore")
    if(!dir.exists(pack_dir))stop("check your package location: ",pack_dir)

    dd="ReporterScore_temp_download"
    dir.create(dd)

    dabiao("Trying to download files from https://rest.kegg.jp/ ")

    download.file("https://rest.kegg.jp/link/pathway/ko",destfile = paste0(dd,"/pathway-KO.list"),method = "curl")
    download.file("https://rest.kegg.jp/list/pathway",destfile = paste0(dd,"/pathway.desc.list"),method = "curl")

    download.file("https://rest.kegg.jp/link/module/ko",destfile = paste0(dd,"/module-KO.list"),method = "curl")
    download.file("https://rest.kegg.jp/list/module",destfile = paste0(dd,"/module.desc.list"),method = "curl")
    make_KO_list(dd,paste0(pack_dir,"/data/new_KOlist.rda"))
    dabiao(paste0("Update done at ",Sys.time()))
}

#' Prepare KOlist
#'
#' @param dir directory of your download file
#' @param output default, dir/KOlist.rda
#'
#' @export
#'
make_KO_list=function(dir,output=NULL){
    dd=dir
    pathway2ko=read.table(paste0(dd,"/pathway-KO.list"),sep = "\t",col.names =c("KO","Pathway"))
    pathway2ko=filter(pathway2ko,grepl("map",Pathway))
    pathway2ko$KO=sub("^ko:","",pathway2ko$KO);pathway2ko$Pathway=sub("^path:","",pathway2ko$Pathway)
    pathway2ko_num=pathway2ko%>%count(Pathway,name = "K_num")
    pathway_desc=read.table(paste0(dd,"/pathway.desc.list"),sep = "\t",col.names =c("Pathway","Description"))
    pathway2ko_com=aggregate(KO~Pathway,pathway2ko,paste, collapse = ",")
    pathway_list=Reduce(left_join, list(pathway2ko_num,pathway2ko_com,pathway_desc))
    colnames(pathway_list)=c("id","K_num","KOs","Description")

    module2ko=read.table(paste0(dd,"//module-KO.list"),sep = "\t",col.names =c("KO","Module"))
    module2ko$KO=sub("^ko:","",module2ko$KO);module2ko$Module=sub("^md:","",module2ko$Module)
    module2ko_num=module2ko%>%count(Module,name = "K_num")
    module_desc=read.table(paste0(dd,"/module.desc.list"),sep = "\t",col.names =c("Module","Description"))
    module2ko_com=aggregate(KO~Module,module2ko,paste, collapse = ",")
    module_list=Reduce(left_join, list(module2ko_num,module2ko_com,module_desc))
    colnames(module_list)=c("id","K_num","KOs","Description")

    KOlist=list("pathway"=pathway_list,"module"=module_list)

    attributes(KOlist)$download_time=file.info(paste0(dd,"/pathway-KO.list"))$mtime
    attributes(KOlist)$build_time=Sys.time()
    if(is.null(output))output=paste0(dd,"/KOlist.rda")
    save(KOlist,file =output)
}
