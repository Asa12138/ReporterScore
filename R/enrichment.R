
#' Perform KO enrichment analysis
#'
#' This function performs KO enrichment analysis using the clusterProfiler package.
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko_test}}.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#' @return A data frame containing the enrichment results.
#' @export
#' @examples
#' a=KO_enrich(ko_pvalue)
#' plot(a)
#'
KO_enrich=function(ko_pvalue,p.adjust.method='BH'){
    lib_ps("clusterProfiler")
    res.dt=ko_pvalue
    if(!all(c("q.value")%in%colnames(res.dt))){stop("check if q.value in your ko_stat dataframe!")}

    load_KOlist()
    path2ko=explode(KOlist$pathway[,c(1,3)],2)
    path2name=KOlist$pathway[,c(1,4)]

    sig_KO=filter(res.dt,q.value<0.05)%>%pull(KO_id)

    e <- enricher(sig_KO,TERM2GENE = path2ko,TERM2NAME = path2name,
                  pAdjustMethod = p.adjust.method, pvalueCutoff = 0.05, qvalueCutoff = 0.05)

    GO_res=as.data.frame(e)
    class(GO_res)<-c("enrich_res",class(GO_res))
    GO_res
}

#' Plot enrich_res
#'
#' @param GO enrich_res object
#' @param mode mode
#' @param str_width default: 50
#'
#' @return ggplot
#' @exportS3Method
#'
plot.enrich_res<-function(GO,mode=1,str_width=50){
    #经典图
    if(mode==1){p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Count, fill=p.adjust))+
        geom_bar(stat = "identity",width=0.7)+####柱子宽度
        scale_fill_gradient(low = "red",high ="blue",limits=c(0,0.05))+#颜色自己可以换
        labs(x = "Gene numbers")+
        scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
        theme_bw()}

    if(mode==2){p=ggplot(data=GO, aes(y=reorder(Description, -p.adjust),x=Count, fill=ONTOLOGY))+
        geom_bar(stat = "identity",width=0.7)+####柱子宽度
        scale_fill_manual(values = c("#66C3A5", "#8DA1CB", "#FD8D62")) + ###颜色
        labs(x = "Gene numbers")+
        scale_y_discrete(labels = \(x)stringr::str_wrap(x, width = str_width))+
        theme_bw()}

    if(length(grep("^GO",GO$ID))>0)p=p+labs(title = "GO Enrichment",y = "GO Term")
    if(length(grep("^hsa|^map|^M",GO$ID))>0)p=p+labs(title = "KEGG Pathways Enrichment",y = "Pathway")
    if(length(grep("^WP",GO$ID))>0)p=p+labs(title = "WiKi Pathways Enrichment",y = "Pathway")
    p
}


#' Perform KO gene set enrichment analysis
#'
#' This function performs KO enrichment analysis using the clusterProfiler package.
#'
#' @param ko_pvalue ko_pvalue dataframe from \code{\link{ko_test}}.
#' @param p.adjust.method The method used for p-value adjustment (default: "BH").
#'
#' @export
#'
#' @examples
#' a=KO_gsea(ko_pvalue)
#' gseaplot(e,geneSetID = e@result$ID[1])
KO_gsea=function(ko_pvalue,p.adjust.method='BH'){
    lib_ps("clusterProfiler")
    vs_group=grep("avg",colnames(ko_pvalue),value = T)
    res.dt=ko_pvalue
    if(!all(c("q.value")%in%colnames(res.dt))){stop("check if q.value in your ko_stat dataframe!")}
    res.dt=filter(res.dt,q.value<0.05)
    res.dt$FC=res.dt[,vs_group[1]]/res.dt[,vs_group[2]]
    res.dt_sort <- res.dt %>% arrange(desc(FC))

    kos <- log2(res.dt_sort$FC)
    names(kos) <- res.dt_sort$KO_id

    load_KOlist()
    path2ko=explode(KOlist$pathway[,c(1,3)],2)
    path2name=KOlist$pathway[,c(1,4)]

    e <- GSEA(kos, TERM2GENE = path2ko, TERM2NAME = path2name,
              verbose=F, pvalueCutoff = 1)
    e
}
