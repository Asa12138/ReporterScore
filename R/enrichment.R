
#' Perform fisher.test for enrichment
#'
#' @rdname KO_enrich
#'
#' @return data.frame
#' @export
#'
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' ko_pvalue=reporter_score_res$ko_pvalue
#' fisher_res=KO_fisher(ko_pvalue)
#' plot(fisher_res)
#' }
KO_fisher=function(ko_pvalue,padj_threshold=0.05,p.adjust.method="BH",type=c("pathway","module")[1],
                   KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    res.dt=ko_pvalue
    if(!all(c("p.adjust")%in%colnames(res.dt))){stop("check if p.adjust in your ko_stat dataframe!")}

    KOlist=NULL
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


    sig_K_num = length(unique(res.dt[res.dt$p.adjust<padj_threshold,]$KO_id))
    nosig_K_num = nrow(res.dt)-sig_K_num

    reps=nrow(modulelist)

    lapply(1:reps,\(i){
        tmp_kos=strsplit(modulelist$KOs[i], ',')[[1]]

        z <- res.dt[res.dt$KO_id %in% tmp_kos,]
        exist_KO=nrow(z)
        significant_KO=sum(z$p.adjust<padj_threshold)

        p_value = stats::fisher.test(matrix(c(significant_KO,exist_KO-significant_KO,
                             sig_K_num-significant_KO,nosig_K_num-(exist_KO-significant_KO)), 2, 2, byrow = T), alternative = "greater")

        c(exist_KO,significant_KO,p_value$p.value)
    })%>%do.call(rbind,.)->res
    colnames(res)=c("Exist_K_num","Significant_K_num","p.value")

    if(verbose)pcutils::dabiao("`fisher.test` done")

    fisher_res <- data.frame(ID = modulelist$id,
                               Description = modulelist$Description,
                               K_num=modulelist$K_num,res)
    fisher_res$p.adjust <- stats::p.adjust(fisher_res$p.value, method = p.adjust.method)
    attributes(fisher_res)$method="fisher.test"
    class(fisher_res)<-c("enrich_res",class(fisher_res))
    fisher_res
}


#' Perform KO enrichment analysis
#'
#' This function performs KO enrichment analysis using the clusterProfiler package.
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link[ReporterScore]{ko.test}}.
#' @param padj_threshold p.adjust threshold to determine whether significant or not.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @param type "pathway" or "module" for default KOlist_file.
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#'
#' @return A data frame containing the enrichment results.
#' @export
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' ko_pvalue=reporter_score_res$ko_pvalue
#' enrich_res=KO_enrich(ko_pvalue)
#' plot(enrich_res)
#' }
KO_enrich=function(ko_pvalue,padj_threshold=0.05,p.adjust.method='BH',type=c("pathway","module")[1],
                   KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    KO_id=p.adjust=NULL
    pcutils::lib_ps("clusterProfiler",library = F)
    res.dt=ko_pvalue
    if(!all(c("p.adjust")%in%colnames(res.dt))){stop("check if p.adjust in your ko_stat dataframe!")}

    KOlist=NULL
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

    path2ko=pcutils::explode(modulelist[,c("id","KOs")],2,split = ",")
    path2name=modulelist[,c("id","Description")]
    #set background
    {path2ko=dplyr::filter(path2ko,KOs%in%res.dt$KO_id)}

    sig_KO=dplyr::filter(res.dt,p.adjust<padj_threshold)%>%dplyr::pull(KO_id)

    e <- clusterProfiler::enricher(sig_KO,TERM2GENE = path2ko,TERM2NAME = path2name,
                  pAdjustMethod = p.adjust.method, pvalueCutoff = 1, qvalueCutoff = 1)

    if(verbose)pcutils::dabiao("`clusterProfiler::enricher` done")

    GO_res=as.data.frame(e)
    GO_res=rename(GO_res,"p.value"="pvalue","Significant_K_num"="Count")
    GO_res=GO_res[,c(1:6,9)]
    class(GO_res)<-c("enrich_res",class(GO_res))
    GO_res
}


#' Plot enrich_res
#'
#' @param x enrich_res object
#' @param mode plot style: 1~2
#' @param str_width default: 50
#' @param ... add
#' @param padj_threshold p.adjust threshold
#'
#' @return ggplot
#' @exportS3Method
#' @method plot enrich_res
#'
plot.enrich_res<-function(x,...,mode=1,str_width=50,padj_threshold=0.05){
    Description=Significant_K_num=NULL
    GO=x
    GO=dplyr::filter(GO,p.adjust<=padj_threshold)
    if(nrow(GO)<1)stop("No pathway p.adjst less than ",padj_threshold)

    #经典图
    if(mode==1){
        p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Significant_K_num, fill=p.adjust))+
            geom_bar(stat = "identity",width=0.7)+#柱子宽度
            scale_fill_gradient(low = "red",high ="blue",limits=c(0,0.05))+#颜色自己可以换
            labs(x = "Gene numbers")+
            scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
            theme_bw()
    }
    if(mode==2){
        p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Significant_K_num, fill=p.adjust,size=Significant_K_num))+
            geom_point(shape=21)+
            scale_size(guide = guide_none())+
            scale_fill_gradient(low = "red",high ="blue",limits=c(0,0.05))+#颜色自己可以换
            labs(x = "Gene numbers")+
            scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
            theme_bw()
    }

    if(length(grep("^GO",GO$ID))>0)p=p+labs(title = "GO Enrichment",y = "GO Term")
    if(length(grep("^hsa|^map|^M",GO$ID))>0)p=p+labs(title = "KEGG Pathways Enrichment",y = "Pathway")
    if(length(grep("^WP",GO$ID))>0)p=p+labs(title = "WiKi Pathways Enrichment",y = "Pathway")
    p
}


#' Perform KO gene set enrichment analysis
#'
#' This function performs KO enrichment analysis using the clusterProfiler package.
#'
#' @rdname KO_enrich
#' @param add_mini add_mini when calculate the logFC. e.g (10+0.1)/(0+0.1), default 0.05*min(avg_abundance)
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' ko_pvalue=reporter_score_res$ko_pvalue
#' gsea_res=KO_gsea(ko_pvalue)
#' enrichplot::gseaplot(gsea_res,geneSetID = gsea_res@result$ID[1])
#' }
KO_gsea=function(ko_pvalue,add_mini=NULL,padj_threshold=0.05,p.adjust.method='BH',type=c("pathway","module")[1],
                 KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    FC=p.adjust=NULL
    pcutils::lib_ps("clusterProfiler",library = F)

    KOlist=NULL
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

    vs_group=grep("average",colnames(ko_pvalue),value = T)
    if(length(vs_group)>2)stop("GESA only available for two groups")
    res.dt=ko_pvalue
    if(!all(c("p.adjust")%in%colnames(res.dt))){stop("check if p.adjust in your ko_stat dataframe!")}
    res.dt=dplyr::filter(res.dt,p.adjust<padj_threshold)

    tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])
    if (is.null(add_mini))
        add_mini = min(tmp[tmp > 0]) * 0.05
    res.dt$FC=(res.dt[,vs_group[1]]+add_mini)/(res.dt[,vs_group[2]]+add_mini)

    res.dt_sort <- res.dt %>% dplyr::arrange(dplyr::desc(FC))

    kos <- log2(res.dt_sort$FC)
    names(kos) <- res.dt_sort$KO_id

    path2ko=pcutils::explode(KOlist$pathway[,c(1,3)],2)
    path2name=KOlist$pathway[,c(1,4)]
    #set background
    {path2ko=dplyr::filter(path2ko,KOs%in%res.dt$KO_id)}

    e <- clusterProfiler::GSEA(kos, TERM2GENE = path2ko, TERM2NAME = path2name,verbose=F,pvalueCutoff = 1)
    e
}
