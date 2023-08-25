reporter_color=c("#e31a1c","#1f78b4","#b15928","#b2df8a","#810f7c","#FFF021")

#' Plot the reporter_res
#'
#' @param reporter_res result of `get_reporter_score`
#' @param rs_threshold plot threshold vector, default:1.64
#' @param y_text_size y_text_size
#' @param str_width str_width to wrap
#' @param mode 1～2 plot style.
#' @param Pathway_description show KO description rather than KO id.
#' @param facet_level facet plot if the type is "pathway" or "module"
#' @param facet_str_width str width for facet label
#' @param facet_anno annotation table for facet, two columns, first is level summary, second is pathway id.
#' @param show_ID show pathway id
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data("reporter_score_res")
#' plot_report(reporter_score_res,rs_threshold=c(2,-2),y_text_size=10,str_width=40)
plot_report<-function(reporter_res,rs_threshold=1.64,mode=1,y_text_size=13,str_width=50,show_ID=FALSE,
                      Pathway_description=TRUE,facet_level=FALSE,facet_anno=NULL,facet_str_width=15){
    if(inherits(reporter_res,"reporter_score"))reporter_res=reporter_res$reporter_s

    flag=FALSE
    if(inherits(reporter_res,"rs_by_cm")){
        rsa_cm_res=reporter_res
        ncluster=sum(grepl("Cluster",names(rsa_cm_res)))
        clusters_name=grep("Cluster",names(rsa_cm_res),value = T)
        reporter_res=lapply(clusters_name,
               \(i){data.frame(rsa_cm_res[[i]]$reporter_s,Cluster=i,row.names = NULL)})%>%
            do.call(rbind,.)
        attributes(reporter_res)=pcutils::update_param(attributes(rsa_cm_res[["Cluster1"]]$reporter_s),attributes(reporter_res))
        flag=TRUE
    }

    reporter_res=na.omit(reporter_res)
    Group=Description=ReporterScore=Exist_K_num=NULL
    if(length(rs_threshold)==1)rs_threshold=c(rs_threshold,-rs_threshold)

    rs_threshold=sort(rs_threshold)
    if(rs_threshold[2]>max((reporter_res$ReporterScore))){
        rs_threshold[2]=tail(sort((reporter_res$ReporterScore)))[1]%>%round(.,4)
        warning("Too big rs_threshold, change rs_threshold to ", rs_threshold[1],"\n")
    }
    vs_group=attributes(reporter_res)$vs_group
    if((attributes(reporter_res)$mode=="directed")&is.null(attributes(reporter_res)$pattern)){
        if(rs_threshold[1]<min((reporter_res$ReporterScore))){
            rs_threshold[1]=head(sort((reporter_res$ReporterScore)))[5]%>%round(.,4)
            warning("Too small rs_threshold, change rs_threshold to", rs_threshold[1],"\n")
        }
        reporter_res2 <- reporter_res[(reporter_res$ReporterScore >= rs_threshold[2])|(reporter_res$ReporterScore <= rs_threshold[1]), ]

        if(length(vs_group)==2){
            reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0,
                                          paste0("Enrich in ",vs_group[2]),
                                          paste0("Enrich in ",vs_group[1]))
            cols1=setNames(c('P'='orange','N'='seagreen'), paste0("Enrich in ",vs_group[2:1]))
            title=paste0(vs_group,collapse = "/ ")
        }
        else {
            reporter_res2$Group <- ifelse(reporter_res2$ReporterScore > 0,"Increase","Decrease")
            cols1=setNames(c('P'='orange','N'='seagreen'), c("Increase","Decrease"))
            title=paste0(vs_group,collapse = "/ ")
        }
        breaks=c(scales::breaks_extended(3)(range(reporter_res2$ReporterScore)),rs_threshold)
    }
    else {
        reporter_res2 <- reporter_res[reporter_res$ReporterScore >= rs_threshold[2],]
        reporter_res2$Group="Significant"
        cols1=c("Significant"='red2')
        title=paste0(vs_group,collapse = "/")
        breaks=c(scales::breaks_extended(3)(range(reporter_res2$ReporterScore)),rs_threshold[2])
    }
    if(flag){
        reporter_res2$Group=reporter_res2$Cluster
        cols1=setNames(pcutils::get_cols(ncluster),clusters_name)
    }

    if(facet_level){
        if(!is.null(facet_anno)){
            tmpdf=facet_anno
            colnames(tmpdf)=c("facet_level","ID")
            reporter_res2=dplyr::left_join(reporter_res2,tmpdf,by=c("ID"))
        }
        else{
            if(is.null(attributes(reporter_res)$type)){
                warning("No attributes(reporter_res)$type found.")
                facet_level=FALSE
            }
            else if(attributes(reporter_res)$type=="pathway"){
                load_Pathway_htable(envir = environment())
                tmpdf=Pathway_htable[,c("level1_name","Pathway_id")]
                colnames(tmpdf)=c("facet_level","ID")
                reporter_res2=dplyr::left_join(reporter_res2,tmpdf,by=c("ID"))
            }
            else if(attributes(reporter_res)$type=="module"){
                load_Module_htable(envir = environment())
                tmpdf=Module_htable[c("module2_name","Module_id")]
                colnames(tmpdf)=c("facet_level","ID")
                reporter_res2=dplyr::left_join(reporter_res2,tmpdf,by=c("ID"))
            }
            else if(attributes(reporter_res)$type=="ALL"){
                tmpdf=reporter_res[c("ONT","ID")]
                colnames(tmpdf)=c("facet_level","ID")
                reporter_res2=dplyr::left_join(reporter_res2,tmpdf,by=c("ID"))
            }
            else {
                load_Pathway_htable(envir = environment())
                tmpdf=Pathway_htable[,c("level1_name","Pathway_id")]
                colnames(tmpdf)=c("facet_level","ID")
                tmpdf$ID=gsub("map",attributes(reporter_res)$type,tmpdf$ID)
                reporter_res2=dplyr::left_join(reporter_res2,tmpdf,by=c("ID"))
            }
        }
    }

    reporter_res2 <- reporter_res2[stats::complete.cases(reporter_res2), ]

    if(show_ID)reporter_res2$Description=paste0(reporter_res2$ID,": ",reporter_res2$Description)
    if(!Pathway_description)reporter_res2$Description=reporter_res2$ID

    if(mode==1){
        if(flag){
            reporter_res2$Description=factor(reporter_res2$Description,
                                             levels = arrange(reporter_res2,Group,ReporterScore)%>%pull(Description))
            p=ggplot(reporter_res2, aes(ReporterScore,Description, fill = Group))
            }
        else p=ggplot(reporter_res2, aes(ReporterScore,stats::reorder(Description, ReporterScore), fill = Group))
        p=p+geom_bar(stat = 'identity', position='dodge')+
            scale_fill_manual(values=cols1)+
            theme_light()
    }
    if(mode==2){
        p=ggplot(reporter_res2, aes(ReporterScore,stats::reorder(Description, ReporterScore),
                                    size=Exist_K_num, fill =Exist_K_num)) +
            geom_point(shape=21)+
            scale_fill_gradient(low = "#FF000033",high = "red",guide = "legend")+theme_light()
    }

    p <-p+labs(y="")+
        geom_vline(xintercept = rs_threshold[2], linetype =2)+
        scale_y_discrete(labels = label_wrap_gen(width = str_width))+
        scale_x_continuous(breaks = breaks)+
        theme(
            axis.text.x = element_text(colour='black',size=13),
            axis.text.y = element_text(colour='black',size=y_text_size)
        )
    if(facet_level)p=p+facet_grid(facet_level~.,scales = "free_y",space = "free",labeller = label_wrap_gen(facet_str_width))+
        theme(strip.text.y = element_text(angle = 0))
    if(attributes(reporter_res)$mode=="directed"&is.null(attributes(reporter_res)$pattern))p=p+geom_vline(xintercept = rs_threshold[1], linetype = 2)
    #if(length(attributes(reporter_res)$vs_group)==2)p=p+labs(title = paste(attributes(reporter_res)$vs_group,collapse = " vs "))
    p=p+labs(title = title)
    return(p)
}


#' Plot KOs trend in one pathway or module
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param map_id the pathway or module id
#' @param select_ko select which ko
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
                             modulelist=NULL,select_ko=NULL,
                             box_color=reporter_color,show_number=TRUE,
                             line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")){
    Group=value=Group2=value2=Group1=value1=type=Significantly=KOlist=p.adjust=NULL
    if(is.null(names(line_color)))names(line_color)=c("Depleted","Enriched","None","Significant")[seq_along(line_color)]
    pcutils::lib_ps("ggnewscale","reshape2",library = FALSE)
    flag=FALSE
    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat
        modulelist=reporter_res$modulelist
        if(is.character(modulelist)){
            load_GOlist(envir = environment())
            modulelist=eval(parse(text = modulelist))
        }
        flag=TRUE
        RS=reporter_res$reporter_s[reporter_res$reporter_s$ID==map_id,"ReporterScore"]
    }

    if(is.null(select_ko))A=get_KOs(map_id =map_id,ko_stat = ko_stat,modulelist =modulelist)
    else A=ko_stat[ko_stat$KO_id%in%select_ko,]

    if(nrow(A)<1)return(NULL)

    if(!all(c("id","K_num","KOs","Description")%in%colnames(modulelist)))stop("check your modulelist format!")
    Description=modulelist[modulelist$id==map_id,"Description"]

    vs_group=attributes(ko_stat)$vs_group
    colnames(A)[colnames(A)%in%paste0("average_",vs_group)]=vs_group
    box_df=reshape2::melt(A[,c("KO_id",vs_group)],id.vars="KO_id",variable.name ="Group")
    box_df$Group=factor(box_df$Group,levels = (vs_group))
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
        #line_color=c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
        names(line_color)=setNames(num$label,num$Significantly)[names(line_color)]
    }

    p=ggplot()+
        geom_boxplot(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        geom_point(data =box_df,aes(x=Group,y=value,color=Group),show.legend = FALSE)+
        scale_color_manual(values = pcutils::get_cols(nlevels(box_df$Group),box_color))+
        labs(title = ifelse(is.null(select_ko),paste0("KOs in ",map_id," (",Description,")"),"Selected KOs"),
             x=NULL,y="Abundance")+
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
                      modulelist = NULL,
                      KO_description=FALSE,str_width=50){
    flag=FALSE
    if(inherits(kodf,"reporter_score")){
        reporter_res=kodf
        kodf=reporter_res$kodf
        group=reporter_res$group
        metadata=reporter_res$metadata
        modulelist=reporter_res$modulelist
        if(is.character(modulelist)){
            load_GOlist(envir = environment())
            modulelist=eval(parse(text = modulelist))
        }
        flag=TRUE
    }

    metadata[,group]=factor(metadata[,group],levels = levels(factor(metadata[,group])))

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,modulelist =modulelist)
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
        load_ko_desc(envir = environment())
        newname=ko_desc[match(rownames(plotdat),ko_desc$KO_id),"KO_name",drop=T]
        if(all(is.na(newname)))warning("No description for KO found, are you sure rownames of kodf are KOs?")
        rownames(plotdat)=ifelse(is.na(newname),rownames(plotdat),newname)%>%stringr::str_wrap(., width = str_width)
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
#' @param KO_description show KO description rather than KO id.
#' @param str_width str_width to wrap
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param only_sig only show the significant KOs
#' @param columns change columns
#'
#' @export
#' @examples
#' data("reporter_score_res")
#' plot_KOs_heatmap(reporter_score_res,map_id="map00780")
plot_KOs_heatmap=function(kodf,group=NULL,metadata=NULL,
                          map_id="map00780",select_ko=NULL,
                          only_sig=FALSE,columns=NULL,
                          modulelist = NULL,
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
        if(is.character(modulelist)){
            load_GOlist(envir = environment())
            modulelist=eval(parse(text = modulelist))
        }
        flag=TRUE
    }

    if(is.null(select_ko))select_ko=get_KOs(map_id =map_id,modulelist =modulelist)

    if(only_sig&flag){
        sig_names=reporter_res$ko_stat%>%dplyr::filter(Significantly!="None")%>%rownames()
        select_ko=intersect(select_ko,sig_names)
    }

    cols=which(rownames(kodf)%in%select_ko)

    metadata[,group]=factor(metadata[,group],levels = levels(factor(metadata[,group])))
    if(length(cols)==0)stop("No select KOs! check map_id or select_ko")
    plotdat=kodf[cols,,drop=FALSE]
    if(!is.null(columns))plotdat=plotdat[,columns]

    annotation_colors=list(pcutils::get_cols(nlevels(factor(metadata[,group])),reporter_color))
    names(annotation_colors)=group
    names(annotation_colors[[group]])=levels(factor(metadata[,group]))

    if(KO_description){
        load_ko_desc(envir = environment())
        newname=ko_desc[match(rownames(plotdat),ko_desc$KO_id),"KO_name",drop=T]
        if(all(is.na(newname)))warning("No description for KO found, are you sure rownames of kodf are KOs?")
        rownames(plotdat)=ifelse(is.na(newname),rownames(plotdat),newname)%>%stringr::str_wrap(., width = str_width)
    }

    do.call(pheatmap::pheatmap,
            pcutils::update_param(list(mat = plotdat,cluster_cols = F,
                                       color = pcutils::get_cols(100,c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")),
                                       annotation_colors=annotation_colors,
                                       annotation_col = metadata[group],scale = "row"),
                                  heatmap_param))
}

#' Plot KOs network
#'
#' @param ko_stat ko_stat result from \code{\link{pvalue2zs}} or result of `get_reporter_score`
#' @param map_id the pathway or module id
#' @param kos_color default, c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
#' @param near_pathway show the near_pathway if any KOs exist.
#' @param modulelist NULL or customized modulelist dataframe, must contain "id","K_num","KOs","Description" columns. Take the `KOlist` as example, use \code{\link{custom_modulelist}}.
#' @param ... additional arguments for \code{\link[MetaNet]{c_net_plot}}
#' @param pathway_label show pathway_label?
#' @param kos_label show kos_label?
#' @param mark_module mark the modules?
#' @param mark_color mark colors, default, c("Depleted"="seagreen","Enriched"="orange","None"="grey","Significant"="red2")
#' @param return_net return the network
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
                          pathway_label=TRUE,kos_label=TRUE,
                          mark_module=FALSE,mark_color=NULL,
                          return_net=FALSE,
                          ...){
    pcutils::lib_ps("ggnewscale","reshape2","MetaNet",library = FALSE)

    if(inherits(ko_stat,"reporter_score")){
        reporter_res=ko_stat
        ko_stat=reporter_res$ko_stat[,c("KO_id","Significantly")]
        modulelist=reporter_res$modulelist
        if(is.character(modulelist)){
            load_GOlist(envir = environment())
            modulelist=eval(parse(text = modulelist))
        }
        reporter_s=reporter_res$reporter_s
    }
    else stop("Need reporter_score object")

    if(is.null(modulelist)){
        load_KOlist(envir = environment())
        if(grepl("map",map_id))modulelist=KOlist$pathway
        if(grepl("M",map_id))modulelist=KOlist$module
        if(grepl("GO:",map_id[1])){
            load_GOlist(envir = environment())
            modulelist=lapply(names(GOlist),function(i)cbind(GOlist[[i]],ONT=i))%>%do.call(rbind,.)
        }
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
    igraph::vertex.attributes(ko_net)[["color"]]=MetaNet::tidai(igraph::vertex.attributes(ko_net)[["v_class"]],kos_color)
    tmp_v=MetaNet::get_v(ko_net)
    if(!pathway_label)tmp_v$label=ifelse(tmp_v$v_group=="Pathway",NA,tmp_v$label)
    if(!kos_label)tmp_v$label=ifelse(tmp_v$v_group=="KOs",NA,tmp_v$label)
    if(mark_module){
        ko_net_m=MetaNet::modu_dect(ko_net,method = "cluster_walktrap")
        if(return_net)return(ko_net_m)
        tmp_v2=MetaNet::get_v(ko_net_m)

        tmp_v2=dplyr::left_join(tmp_v2,reporter_s,by=c("name"="ID"))
        modules=dplyr::group_by(tmp_v2,module)%>%dplyr::summarise(RS=mean(ReporterScore,na.rm = T))%>%as.data.frame()
        if(is.null(mark_color))mark_color=kos_color
        if(attributes(reporter_res$ko_stat)$mode=="directed"){
            modules$color=ifelse(modules$RS>1.64,kos_color["Enriched"],
                                 ifelse(modules$RS<(-1.64),kos_color["Depleted"],kos_color["None"]))
        }
        else modules$color=ifelse(modules$RS>1.64,kos_color["Significant"],kos_color["None"])
        modules_col=setNames(modules$color,modules$module)
        plot(ko_net_m,vertex.color=kos_color,vertex.label=tmp_v$label,mark_module=T,mark_color=modules_col,...)
    }
    else {
        if(return_net)return(ko_net)
        plot(ko_net,vertex.color=kos_color,vertex.label=tmp_v$label,...)
    }
}

#' plot KO_htable levels
#'
#' @param select_ko select kos
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data("KO_abundance_test")
#' plot_KO_htable(rownames(KO_abundance))
plot_KO_htable=function(select_ko=NULL){
    load_KO_htable(envir = environment())
    if(!is.null(select_ko))KO_htable=dplyr::filter(KO_htable,KO_id%in%select_ko)

    dplyr::count(KO_htable,level1_name,level2_name)->a
    a$level2_name=factor(a$level2_name,levels = a$level2_name)
    patt=setNames(pcutils::get_cols(length(unique(a$level1_name))),unique(a$level1_name))

    ggplot(a)+
        geom_col(aes(x=n,y=level2_name,fill=level1_name))+
        scale_fill_manual(values = patt)+
        geom_text(aes(x=n,y=level2_name,label=n),
                  nudge_x = max(a$n)*1/200,hjust=0)+
        labs(y=NULL,x="Number of KOs",title = "KEGG pathway annotation")+
        theme_classic()+
        theme(axis.text.y = element_text(color = patt[a$level1_name]))+
        scale_x_continuous(expand = c(0,0),limits = c(0,max(a$n)*1.1))
}

#' Plot c_means result
#'
#' @param rsa_cm_res a cm_res object
#' @param ... additional
#' @param filter_membership filter membership 0~1.
#' @param mode 1~2
#' @param show.clust.cent show cluster center?
#' @param show_num show number of each cluster?
#'
#' @export
#'
plot_c_means=function(rsa_cm_res,filter_membership,mode=1,show.clust.cent=TRUE,show_num=TRUE,...){
    stopifnot(inherits(rsa_cm_res,"rs_by_cm"))
    plot(rsa_cm_res$cm_res,filter_membership = filter_membership,show.clust.cent=show.clust.cent,mode=mode,show_num=show_num,...)
}
