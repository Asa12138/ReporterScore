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
