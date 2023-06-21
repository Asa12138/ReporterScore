
#' Perform KO enrichment analysis
#'
#' This function performs KO enrichment analysis using the clusterProfiler package.
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link[ReporterScore]{ko.test}}.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @param type "pathway" or "module"
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
KO_enrich=function(ko_pvalue,p.adjust.method='BH',type=c("pathway","module")[1],
                   KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    KO_id=q.value=NULL
    pcutils::lib_ps("clusterProfiler",library = F)
    res.dt=ko_pvalue
    if(!all(c("q.value")%in%colnames(res.dt))){stop("check if q.value in your ko_stat dataframe!")}

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

    sig_KO=filter(res.dt,q.value<0.05)%>%pull(KO_id)

    e <- clusterProfiler::enricher(sig_KO,TERM2GENE = path2ko,TERM2NAME = path2name,
                  pAdjustMethod = p.adjust.method, pvalueCutoff = 0.05, qvalueCutoff = 0.05)

    GO_res=as.data.frame(e)
    class(GO_res)<-c("enrich_res",class(GO_res))
    GO_res
}

#' Plot enrich_res
#'
#' @param x enrich_res object
#' @param mode plot style
#' @param str_width default: 50
#' @param ... add
#'
#' @return ggplot
#' @exportS3Method
#' @method plot enrich_res
#'
plot.enrich_res<-function(x,...,mode=1,str_width=50){
    Description=Count=NULL
    GO=x
    #经典图
    if(mode==1){
        p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Count, fill=p.adjust))+
            geom_bar(stat = "identity",width=0.7)+####柱子宽度
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
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko.test}}.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @param type "pathway" or "module"
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param verbose logical
#'
#' @export
#'
#' @examples
#' \donttest{
#' data("reporter_score_res")
#' ko_pvalue=reporter_score_res$ko_pvalue
#' enrich_res2=KO_gsea(ko_pvalue)
#' enrichplot::gseaplot(enrich_res2,geneSetID = enrich_res2@result$ID[1])
#' }
KO_gsea=function(ko_pvalue,p.adjust.method='BH',type=c("pathway","module")[1],
                 KOlist_file=NULL,modulelist=NULL,verbose=TRUE){
    FC=q.value=NULL
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

    vs_group=grep("avg",colnames(ko_pvalue),value = T)
    res.dt=ko_pvalue
    if(!all(c("q.value")%in%colnames(res.dt))){stop("check if q.value in your ko_stat dataframe!")}
    res.dt=filter(res.dt,q.value<0.05)
    res.dt$FC=res.dt[,vs_group[1]]/res.dt[,vs_group[2]]
    res.dt_sort <- res.dt %>% arrange(desc(FC))

    kos <- log2(res.dt_sort$FC)
    names(kos) <- res.dt_sort$KO_id

    path2ko=pcutils::explode(KOlist$pathway[,c(1,3)],2)
    path2name=KOlist$pathway[,c(1,4)]

    e <- clusterProfiler::GSEA(kos, TERM2GENE = path2ko, TERM2NAME = path2name,verbose=F, pvalueCutoff = 1)
    e
}
