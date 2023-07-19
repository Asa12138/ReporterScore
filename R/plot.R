reporter_color=c("#e31a1c","#1f78b4","#b15928","#b2df8a","#810f7c")

#' Plot the reporter_res
#'
#' @param reporter_res result of `get_reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param y_text_size y_text_size
#' @param str_width str_width to wrap
#' @param mode 1～2 plot style.
#' @param Pathway_description show KO description rather than KO id.
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_report(reporter_score_res,rs_threshold=c(2,-2),y_text_size=10,str_width=40)
plot_report<-function(reporter_res,rs_threshold=1.64,mode=1,y_text_size=13,str_width=50,Pathway_description=TRUE){
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
        reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0, 'P', 'N')
        breaks=c(scales::breaks_extended(3)(range(reporter_res2$ReporterScore)),rs_threshold)
    }
    else {
        reporter_res2 <- reporter_res[reporter_res$ReporterScore >= rs_threshold[2],]
        reporter_res2$Group="S"
        breaks=c(scales::breaks_extended(3)(range(reporter_res2$ReporterScore)),rs_threshold[2])
    }
    reporter_res2 <- reporter_res2[stats::complete.cases(reporter_res2), ]
    if(!Pathway_description)reporter_res2$Description=reporter_res2$ID

    if(mode==1){
        p=ggplot(reporter_res2, aes(ReporterScore,stats::reorder(Description, ReporterScore), fill = Group)) +
            geom_bar(stat = 'identity', position='dodge')+
            scale_fill_manual(values=c('P'='orange','N'='seagreen',"S"='red2'))+
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
        scale_x_continuous(breaks = breaks)+
        theme(
            axis.text.x = element_text(colour='black',size=13),
            axis.text.y = element_text(colour='black',size=y_text_size)
        )

    if(attributes(reporter_res)$mode=="directed")p=p+geom_vline(xintercept = rs_threshold[1], linetype = 2)
    if(length(attributes(reporter_res)$vs_group)==2)p=p+labs(title = paste(attributes(reporter_res)$vs_group,collapse = " vs "))
    else p=p+labs(title = paste(attributes(reporter_res)$vs_group,collapse = "/ "))
    return(p)
}

#' Plot KOs trend in one pathway or module
#'
#' @param map_id the pathway or module id
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param box_color box and point color, default: c("#e31a1c","#1f78b4")
#' @param line_color line color, default: c("Depleted"="seagreen","Enriched"="orange","None"="grey")
#' @param show_number show the numbers.
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_KOs_in_pathway(ko_stat = reporter_score_res,map_id="map00860")
plot_KOs_in_pathway=function(ko_stat,map_id="map00780",
                             KOlist_file=NULL,modulelist=NULL,
                             box_color=reporter_color,show_number=TRUE,
                             line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")){
    Group=value=Group2=value2=Group1=value1=type=Significantly=KOlist=p.adjust=NULL
    pcutils::lib_ps("ggnewscale","reshape2",library = FALSE)
    flag=FALSE
    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat
        modulelist=reporter_res$modulelist
        flag=TRUE
        RS=reporter_res$reporter_s[reporter_res$reporter_s$ID==map_id,"ReporterScore"]
    }

    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        if(grepl("map",map_id))modulelist=KOlist$pathway
        if(grepl("M",map_id))modulelist=KOlist$module
    }

    A=get_KOs(map_id =map_id,ko_stat = ko_stat,KOlist_file =KOlist_file,modulelist =modulelist )
    if(nrow(A)<1)return(NULL)

    if(any(colnames(modulelist)!=c("id","K_num","KOs","Description")))stop("check your KOlist or modulelist format!")
    Description=modulelist[modulelist$id==map_id,"Description"]

    vs_group=grep("average",colnames(A),value = TRUE)
    box_df=reshape2::melt(A[,c("KO_id",vs_group)],id.vars="KO_id",variable.name ="Group")
    box_df$Group=factor(box_df$Group,levels = rev(vs_group))
    line_df=data.frame()

    for (i in 1:(length(vs_group)-1)) {
        tmp=A[,c("KO_id",vs_group[i:(i+1)],"Significantly")]
        colnames(tmp)=c("KO_id","value1","value2","Significantly")
        tmp$Group1=vs_group[i]; tmp$Group2=vs_group[i+1]
        line_df=rbind(line_df,tmp)
    }
    if(show_number){
        num=dplyr::count(line_df,Significantly)%>%dplyr::mutate(label=paste0(Significantly,": ",n))
        line_df$Significantly=setNames(num$label,num$Significantly)[line_df$Significantly]
        line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
        names(line_color)=setNames(num$label,num$Significantly)[names(line_color)]
    }

    p=ggplot()+
        geom_boxplot(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        geom_point(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        scale_color_manual(values = pcutils::get_cols(nlevels(box_df$Group),box_color))+
        labs(title = paste0("KOs in ",map_id," (",Description,")"),x=NULL,y="Abundance")+
        ggnewscale::new_scale_color()+
        geom_segment(data = line_df,
                     aes(x=Group2,y=value2,xend=Group1,yend=value1,color=Significantly))+
        scale_color_manual(values = line_color)+
        theme_classic(base_size = 13)+theme(axis.text = element_text(color = "black"))
    if(flag)p=p+labs(subtitle = paste0("ReporterScore: ",round(RS,3)))
    p
}

#' Plot KOs boxplot
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples. or result of `get_reporter_score`
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param box_param parameters pass to \code{\link[pcutils]{group_box}}
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param KO_description show KO description rather than KO id.
#' @param str_width str_width to wrap
#' @param only_sig only show the significant KOs
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_box(reporter_score_res,"Group",metadata,
#'      select_ko=c("K00059","K00208","K00647","K00652","K00833","K01012"))
#' plot_KOs_box(reporter_score_res,"Group",metadata,select_ko="K00059",KO_description=TRUE)
plot_KOs_box=function(kodf,group=NULL,metadata=NULL,
                      map_id="map00780",select_ko=NULL,only_sig=FALSE,
                      box_param=NULL,
                      KOlist_file = NULL,modulelist = NULL,
                      KO_description=FALSE,str_width=50){
    flag=FALSE
    if(inherits(kodf,"reporter_score")){
        reporter_res=kodf
        kodf=reporter_res$kodf
        group=reporter_res$group
        metadata=reporter_res$metadata
        modulelist=reporter_res$modulelist
        flag=TRUE
    }

    metadata[,group]=factor(metadata[,group],levels = rev(levels(factor(metadata[,group]))))

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,KOlist_file =KOlist_file,modulelist =modulelist )
    if(only_sig&flag){
        sig_names=reporter_res$ko_stat%>%dplyr::filter(Significantly!="None")%>%rownames()
        select_ko=intersect(select_ko,sig_names)
    }
    tkodf=kodf%>%t()%>%as.data.frame()
    cols=which(colnames(tkodf)%in%select_ko)

    if(length(cols)==0)stop("No select KOs! check map_id or select_ko")
    if(length(cols)>36){
        print(("Too many KOs, do you still want to plot?"))
        flag=readline("yes/no(y/n)?")
        if(!tolower(flag)%in%c("yes","y"))return(NULL)
    }

    plotdat=tkodf[,cols,drop=FALSE]

    if(KO_description){
        load_KOlist(envir = environment())
        colnames(plotdat)=ko_desc[match(colnames(plotdat),ko_desc$KO),"Description"]%>%stringr::str_wrap(., width = str_width)
    }

    do.call(pcutils::group_box,
            append(list(tab = plotdat,group = group,metadata = metadata),
                   pcutils::update_param(list(p_value1 = TRUE,trend_line = TRUE),box_param)))+
        theme_classic(base_size = 13)+theme(axis.text = element_text(color = "black"))+
        scale_fill_manual(values = pcutils::get_cols(nlevels(metadata[,group]),reporter_color))+
        scale_color_manual(values = pcutils::get_cols(nlevels(metadata[,group]),reporter_color))
}

#' Plot KOs heatmap
#'
#' @param kodf KO_abundance table, rowname is ko id (e.g. K00001),colnames is samples. or result of `get_reporter_score`
#' @param group The compare group (two category) in your data, one column name of metadata when metadata exist or a vector whose length equal to columns number of kodf.
#' @param metadata metadata
#' @param map_id the pathway or module id
#' @param select_ko select which ko
#' @param heatmap_param parameters pass to \code{\link[pheatmap]{pheatmap}}
#' @param KOlist_file default NULL, use the internal file. Or you can upload your .rda file from \code{\link{make_KO_list}}
#' @param KO_description show KO description rather than KO id.
#' @param str_width str_width to wrap
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param only_sig only show the significant KOs
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_heatmap(reporter_score_res,map_id="map00780")
plot_KOs_heatmap=function(kodf,group=NULL,metadata=NULL,
                          map_id="map00780",select_ko=NULL,only_sig=FALSE,
                          KOlist_file = NULL,modulelist = NULL,
                          KO_description=FALSE,str_width=50,
                          heatmap_param=list()){
    pcutils::lib_ps("pheatmap",library = FALSE)
    flag=FALSE
    if(inherits(kodf,"reporter_score")){
        reporter_res=kodf
        kodf=reporter_res$kodf
        group=reporter_res$group
        metadata=reporter_res$metadata
        modulelist=reporter_res$modulelist
        flag=TRUE
    }

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,KOlist_file =KOlist_file,modulelist =modulelist)

    if(only_sig&flag){
        sig_names=reporter_res$ko_stat%>%dplyr::filter(Significantly!="None")%>%rownames()
        select_ko=intersect(select_ko,sig_names)
    }

    cols=which(rownames(kodf)%in%select_ko)

    metadata[,group]=factor(metadata[,group],levels = rev(levels(factor(metadata[,group]))))
    if(length(cols)==0)stop("No select KOs! check map_id or select_ko")
    plotdat=kodf[cols,,drop=FALSE]

    annotation_colors=list(pcutils::get_cols(nlevels(factor(metadata[,group])),reporter_color))
    names(annotation_colors)=group
    names(annotation_colors[[group]])=levels(factor(metadata[,group]))

    if(KO_description){
        load_KOlist(envir = environment())
        rownames(plotdat)=ko_desc[match(rownames(plotdat),ko_desc$KO),"Description"]%>%stringr::str_wrap(., width = str_width)
    }

    do.call(pheatmap::pheatmap,
            pcutils::update_param(list(mat = plotdat,cluster_cols = F,
                                       color = pcutils::get_cols(100,c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")),
                                       annotation_colors=annotation_colors,
                                       annotation_col = metadata[group],scale = "row"),
                                  heatmap_param))
}

#' Plot KOs heatmap
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param map_id the pathway or module id
#' @param kos_color default, c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
#' @param near_pathway show the near_pathway if any KOs exist.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param ... additional arguments for \code{\link[MetaNet]{c_net_plot}}
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_network(reporter_score_res,map_id="map05230")
#' plot_KOs_network(reporter_score_res,map_id="map00780",near_pathway=TRUE)
plot_KOs_network=function(ko_stat,map_id="map00780",
                          near_pathway=FALSE,
                          modulelist=NULL,
                          kos_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2","Pathway"="#80b1d3"),
                          ...){
    pcutils::lib_ps("ggnewscale","reshape2","MetaNet",library = FALSE)

    flag=FALSE
    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat[,c("KO_id","Significantly")]
        modulelist=reporter_res$modulelist
        flag=TRUE
        reporter_s=reporter_res$reporter_s[,c("ID","ReporterScore")]
    }
    else stop("Need reporter_score object")

    if(is.null(modulelist)){
        load_KOlist(KOlist_file,envir = environment())
        if(grepl("map",map_id))modulelist=KOlist$pathway
        if(grepl("M",map_id))modulelist=KOlist$module
    }

    id2ko=modulelist%>%dplyr::select(id,KOs)%>%pcutils::explode(2)

    select_ko=get_KOs(map_id=map_id,modulelist =modulelist)

    select_ko=intersect(select_ko,ko_stat$KO_id)
    if(length(select_ko)==0)stop("No select KOs! check map_id or select_ko")

    if(!near_pathway){
        id2ko=dplyr::filter(id2ko,id%in%map_id)
    }

    id2ko=dplyr::filter(id2ko,KOs%in%select_ko)
    colnames(id2ko)=c("Pathway","KOs")

    ko_net=MetaNet::twocol_edgelist(id2ko)
    ko_net=MetaNet::c_net_set(ko_net,ko_stat,vertex_class = "Significantly")
    igraph::graph.attributes(ko_net)$n_type="ko_net"
    plot(ko_net,vertex.color=kos_color,...)
}
