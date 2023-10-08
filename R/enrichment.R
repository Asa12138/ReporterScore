cid2keggid=function(id=c("5281127"),
                        from="PubChem CID",
                        to=c("KEGG"),check_api=FALSE){
    #install_github("dgrapov/CTSgetR")

    lib_ps("CTSgetR","httr",library = F)

    #Make sure CTS API is available
    if(check_api){
        httr::GET('https://cts.fiehnlab.ucdavis.edu/services') %>%
            httr::http_status(.) %>%
            {if( .$category != 'Success'){stop('Oops looks like https://cts.fiehnlab.ucdavis.edu/services is down!') }}
    }
    # db_name<-'ctsgetr.sqlite'
    # CTSgetR::init_CTSgetR_db(db_name)
    # CTSgetR::db_stats()

    suppressMessages({res=lapply(id, \(i)CTSgetR::CTSgetR(i,from,to))%>%do.call(rbind,.)})%>%suppressWarnings()
    #KEGGREST::keggConv("compound","pubchem")
    res
}


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
#'
#' @return A data frame containing the enrichment results.
#' @export
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' enrich_res=KO_enrich(reporter_score_res)
#' plot(enrich_res)
#' }
KO_enrich=function(ko_stat,padj_threshold=0.05,
                   logFC_threshold=NULL,add_mini=NULL,p.adjust.method='BH',
                   type=c("pathway","module")[1],feature="ko",
                   modulelist=NULL,verbose=TRUE){
    KO_enrich_internal(ko_stat,padj_threshold,
                       logFC_threshold,add_mini,p.adjust.method,
                       type,feature,
                       modulelist,verbose,mode=1)
}

KO_enrich_internal=function(ko_stat,padj_threshold=0.05,
                   logFC_threshold=NULL,add_mini=NULL,p.adjust.method='BH',
                   type=c("pathway","module")[1],feature="ko",
                   modulelist=NULL,verbose=TRUE,mode=1,weight="logFC"){

    KO_id=p.adjust=NULL
    pcutils::lib_ps("clusterProfiler",library = F)
    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat
        modulelist=reporter_res$modulelist
        if(is.character(modulelist)){
            load_GOlist(envir = environment())
            modulelist=eval(parse(text = modulelist))
        }
        type=attributes(reporter_res$reporter_s)$type
    }
    res.dt=ko_stat

    KOlist=NULL
    if(is.null(modulelist)){
        modulelist=get_modulelist(type,feature,verbose)
    }
    if(!all(c("id","K_num","KOs","Description")%in%colnames(modulelist)))stop("check your modulelist format!")


    if(!all(c("KO_id","p.adjust")%in%colnames(res.dt))){stop("check if p.adjust in your ko_stat dataframe!")}
    if("origin_p.adjust"%in%colnames(res.dt)){
        message("detect the origin_p.adjust, use the origin_p.adjust.")
        res.dt$p.adjust=res.dt$origin_p.adjust
    }

    #modulist内的KO才考虑
    all_KOs=lapply(modulelist$KOs,\(i)strsplit(i,",")[[1]])%>%do.call(c,.)
    res.dt=dplyr::filter(res.dt,KO_id%in%all_KOs)

    if(!is.null(logFC_threshold)){
        vs_group=grep("average",colnames(res.dt),value = T)
        if(length(vs_group)!=2)stop("logFC only available for two groups")
        tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])

        if (is.null(add_mini))
            add_mini = min(tmp[tmp > 0]) * 0.05
        res.dt$logFC=log2((res.dt[,vs_group[2]]+add_mini)/(res.dt[,vs_group[1]]+add_mini))
        sig_KO=dplyr::filter(res.dt,p.adjust<padj_threshold,abs(logFC)>logFC_threshold)%>%dplyr::pull(KO_id)
    }
    else sig_KO=dplyr::filter(res.dt,p.adjust<padj_threshold)%>%dplyr::pull(KO_id)%>%unique()

    if(length(sig_KO)<1)return(NULL)

    path2ko=pcutils::explode(modulelist[,c("id","KOs")],2,split = ",")
    path2name=modulelist[,c("id","Description")]
    #set background
    #这个跟指定universe的结果一致
    {path2ko=dplyr::filter(path2ko,KOs%in%res.dt$KO_id)}

    if(mode==1){
        e <- clusterProfiler::enricher(gene = sig_KO,TERM2GENE = path2ko,TERM2NAME = path2name,
                                       pAdjustMethod = p.adjust.method, pvalueCutoff = 1, qvalueCutoff = 1)

        if(verbose)pcutils::dabiao("`clusterProfiler::enricher` done")

        GO_res=as.data.frame(e)
        GO_res=rename(GO_res,"p.value"="pvalue","Significant_K_num"="Count")
        GO_res=GO_res[,c(1:6,9)]
        GO_res$Exist_K_num=pcutils::strsplit2(GO_res$BgRatio,"/")[,1]%>%as.numeric()
        GO_res$Exist_K_num=as.integer(GO_res$Exist_K_num)
        GO_res$Significant_K_num=as.integer(GO_res$Significant_K_num)
        class(GO_res)<-c("enrich_res",class(GO_res))
        attributes(GO_res)$method="enricher"
        attributes(GO_res)$type=type
        return(GO_res)
    }
    else if (mode==2){
        sig_K_num = length(sig_KO)
        nosig_K_num = nrow(res.dt)-sig_K_num

        reps=nrow(modulelist)

        lapply(1:reps,\(i){
            tmp_kos=strsplit(modulelist$KOs[i], ',')[[1]]

            z <- res.dt[res.dt$KO_id %in% tmp_kos,]%>%dplyr::distinct(KO_id,.keep_all = T)
            exist_KO=nrow(z)
            significant_KO=sum(z$KO_id%in%sig_KO)

            p_value = stats::fisher.test(matrix(c(significant_KO,exist_KO-significant_KO,
                                                  sig_K_num-significant_KO,nosig_K_num-(exist_KO-significant_KO)), 2, 2, byrow = T),
                                         alternative = "greater")

            c(exist_KO,significant_KO,p_value$p.value)
        })%>%do.call(rbind,.)->res
        colnames(res)=c("Exist_K_num","Significant_K_num","p.value")

        if(verbose)pcutils::dabiao("`fisher.test` done")

        fisher_res <- data.frame(ID = modulelist$id,
                                 Description = modulelist$Description,
                                 K_num=modulelist$K_num,res)

        fisher_res=dplyr::filter(fisher_res,Exist_K_num>0)%>%dplyr::arrange(p.value)
        fisher_res$p.adjust <- stats::p.adjust(fisher_res$p.value, method = p.adjust.method)

        fisher_res$Exist_K_num=as.integer(fisher_res$Exist_K_num)
        fisher_res$Significant_K_num=as.integer(fisher_res$Significant_K_num)
        class(fisher_res)<-c("enrich_res",class(fisher_res))
        attributes(fisher_res)$method="fisher.test"
        attributes(fisher_res)$type=type
        return(fisher_res)
    }
    else if(mode==3){
        if(!weight%in%colnames(res.dt)){
            vs_group=grep("average",colnames(res.dt),value = T)
            if(length(vs_group)!=2)stop("logFC only available for two groups")
            tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])

            if (is.null(add_mini))
                add_mini = min(tmp[tmp > 0]) * 0.05
            res.dt$logFC=log2((res.dt[,vs_group[2]]+add_mini)/(res.dt[,vs_group[1]]+add_mini))
            res.dt_sort <- res.dt %>% dplyr::arrange(dplyr::desc(logFC))
            kos <- res.dt_sort$logFC
            names(kos) <- res.dt_sort$KO_id
        }
        else{
            kos=res.dt[,weight,drop=T]
            names(kos) <- res.dt$KO_id
            kos=sort(kos,decreasing = T)
        }
        e <- clusterProfiler::GSEA(kos, TERM2GENE = path2ko, TERM2NAME = path2name,verbose=F,pvalueCutoff = 1)
        if(verbose)pcutils::dabiao("`clusterProfiler::GSEA` done")
        attributes(e)$type=type
        return(e)
    }
}

#' @rdname KO_enrich
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' fisher_res=KO_fisher(reporter_score_res)
#' plot(fisher_res)
#' }
KO_fisher=function(ko_stat,padj_threshold=0.05,
                   logFC_threshold=NULL,add_mini=NULL,p.adjust.method='BH',
                   type=c("pathway","module")[1],feature="ko",
                   modulelist=NULL,verbose=TRUE){
    KO_enrich_internal(ko_stat,padj_threshold,
                       logFC_threshold,add_mini,p.adjust.method,
                       type,feature,
                       modulelist,verbose,mode=2)
}


#' as enrich_res object
#'
#' @param gsea_res gsea_res from KO_gsea
#'
#' @return enrich_res object
#' @export
#'
as.enrich_res=function(gsea_res){
    res=gsea_res@result
    class(res)<-c("enrich_res",class(res))
    res$Exist_K_num=res$setSize
    attributes(res)$type=attributes(gsea_res)$type
    res
}

#' Plot enrich_res
#'
#' @param enrich_res enrich_res object
#' @param mode plot style: 1~2
#' @param str_width default: 50
#' @param ... add
#' @param padj_threshold p.adjust threshold
#' @param facet_level facet plot if the type is "pathway" or "module"
#' @param facet_str_width str width for facet label
#' @param facet_anno annotation table for facet, two columns, first is level summary, second is pathway id.
#'
#' @return ggplot
#' @export
#'
plot_enrich_res<-function(enrich_res,mode=1,str_width=50,padj_threshold=0.05,
                          facet_level=FALSE,facet_anno=NULL,facet_str_width=15,...){
    Description=Significant_K_num=NULL
    flag=FALSE
    if(inherits(enrich_res,"enrich_res"))GO=enrich_res
    else if(is.list(enrich_res)&all(sapply(enrich_res,\(i)inherits(i,"enrich_res")))){
        multi_enrich_res=enrich_res
        if(is.null(names(multi_enrich_res)))names(multi_enrich_res)=paste0("Res",seq_along(multi_enrich_res))
        GO=lapply(names(multi_enrich_res),
                            \(i){data.frame(multi_enrich_res[[i]][,c("ID","Description","Exist_K_num","Significant_K_num","p.adjust")],Group=i,row.names = NULL)})%>%
            do.call(rbind,.)
        attributes(GO)$type=attributes(multi_enrich_res[[1]])$type
        flag=TRUE
    }
    else GO=enrich_res

    GO=dplyr::filter(GO,p.adjust<=padj_threshold)
    if(nrow(GO)<1)stop("No pathway p.adjst less than ",padj_threshold)

    if(facet_level){
        tmpdf=get_facet_anno(GO,facet_anno)
        GO=dplyr::left_join(GO,tmpdf,by=c("ID"))
    }

    if(flag){
        #经典图
        if(mode==1){
            p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=-log(p.adjust), fill=Group))+
                geom_bar(stat = "identity",width=0.7, position='dodge')+#柱子宽度
                scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
                theme_bw()
        }
        if(mode==2){
            p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=-log(p.adjust), fill=Group, size=Significant_K_num))+
                geom_point(shape=21,  position=position_dodge(width = 0.7))+
                scale_size(range = c(3,8))+
                scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
                theme_bw()
        }
    }
    else {
        #经典图
        if(mode==1){
            p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Significant_K_num/Exist_K_num, fill=p.adjust))+
                geom_bar(stat = "identity",width=0.7)+#柱子宽度
                scale_fill_gradient(low = "red",high ="blue",limits=c(0,padj_threshold))+#颜色自己可以换
                labs(x = "Feature ratio")+
                scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
                theme_bw()
        }
        if(mode==2){
            p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Significant_K_num/Exist_K_num,
                                  fill=p.adjust, size=Significant_K_num))+
                geom_point(shape=21)+
                scale_size(range = c(3,8),)+
                scale_fill_gradient(low = "red",high ="blue",limits=c(0,padj_threshold))+#颜色自己可以换
                labs(x = "Feature ratio")+
                scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
                theme_bw()
        }
    }
    if(facet_level)p=p+facet_grid(facet_level~.,scales = "free_y",space = "free",labeller = label_wrap_gen(facet_str_width))+
        theme(strip.text.y = element_text(angle = 0))

    if(length(grep("^GO",GO$ID))>0)p=p+labs(title = "GO Enrichment",y = "GO Term")
    if(length(grep("^hsa|^map|^M",GO$ID))>0)p=p+labs(title = "KEGG Pathways Enrichment",y = "Pathway")
    if(length(grep("^WP",GO$ID))>0)p=p+labs(title = "WiKi Pathways Enrichment",y = "Pathway")
    p
}

#' Plot enrich_res
#'
#' @param x enrich_res object
#' @param mode plot style: 1~2
#' @param str_width default: 50
#' @param ... add
#' @param padj_threshold p.adjust threshold
#' @param facet_level facet plot if the type is "pathway" or "module"
#' @param facet_str_width str width for facet label
#' @param facet_anno annotation table for facet, two columns, first is level summary, second is pathway id.
#'
#' @return ggplot
#' @exportS3Method
#' @method plot enrich_res
plot.enrich_res<-function(x,mode=1,str_width=50,padj_threshold=0.05,
                          facet_level=FALSE,facet_anno=NULL,facet_str_width=15,...){
    plot_enrich_res(x,mode,str_width,padj_threshold,
                    facet_level,facet_anno,facet_str_width,...)
}

#' Perform KO gene set enrichment analysis
#'
#' @param weight the metric used for ranking, default: logFC
#' @rdname KO_enrich
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' gsea_res=KO_gsea(reporter_score_res)
#' enrichplot::gseaplot(gsea_res,geneSetID = gsea_res@result$ID[1])
#' plot(as.enrich_res(gsea_res@result))
#' }
KO_gsea=function(ko_stat,weight="logFC",add_mini=NULL,
                 padj_threshold=1,p.adjust.method='BH',
                 type=c("pathway","module")[1],feature="ko",
                 modulelist=NULL,verbose=TRUE){
    KO_enrich_internal(ko_stat,padj_threshold,
                       logFC_threshold=NULL,add_mini,p.adjust.method,
                       type,feature,
                       modulelist,verbose,mode=3,weight)
}

G_gsea=function(){

}
