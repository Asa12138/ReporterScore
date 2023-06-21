#' Update the KO2Pathway and KO2Module files
#'
#' @param download_dir the ReporterScore user dir location, detect automatically by `tools::R_user_dir("ReporterScore")`.
#'
#' @export
#' @description
#' Download links:
#'
#' \code{https://rest.kegg.jp/link/pathway/ko}
#' \code{https://rest.kegg.jp/list/pathway}
#' \code{https://rest.kegg.jp/link/module/ko}
#' \code{https://rest.kegg.jp/list/module}
#' \code{https://rest.kegg.jp/list/ko}
#'
#' you can also download yourself the use \code{\link{make_KO_list}} to get a KOlist object
update_KO_file=function(download_dir=NULL){
    #if(is.null(pack_dir))pack_dir=paste0(.libPaths()[1],"/ReporterScore")
    pack_dir=tools::R_user_dir("ReporterScore")
    if(!dir.exists(pack_dir))stop("check your package user dir location: ",pack_dir)

    if(is.null(pack_dir))dd="ReporterScore_temp_download"
    else dd=download_dir

    if(!dir.exists(dd))dir.create(dd)

    pcutils::dabiao("Trying to download files from https://rest.kegg.jp/ ")

    utils::download.file("https://rest.kegg.jp/link/pathway/ko",destfile = paste0(dd,"/pathway-KO.list"),method = "curl")
    utils::download.file("https://rest.kegg.jp/list/pathway",destfile = paste0(dd,"/pathway.desc.list"),method = "curl")

    utils::download.file("https://rest.kegg.jp/link/module/ko",destfile = paste0(dd,"/module-KO.list"),method = "curl")
    utils::download.file("https://rest.kegg.jp/list/module",destfile = paste0(dd,"/module.desc.list"),method = "curl")

    utils::download.file("https://rest.kegg.jp/list/ko",destfile = paste0(dd,"/ko.desc.list"),method = "curl")

    make_KO_list(dd,paste0(pack_dir,"/new_KOlist.rda"))

    pcutils::dabiao(paste0("Update done at ",Sys.time()))
}

#' Prepare KOlist
#'
#' @param download_dir directory of your download file
#' @param output default, download_dir/KOlist.rda
#'
#' @export
make_KO_list=function(download_dir,output=NULL){
    Pathway=Module=NULL
    dd=download_dir
    pathway2ko=utils::read.table(paste0(dd,"/pathway-KO.list"),sep = "\t",col.names =c("KO","Pathway"))
    pathway2ko=dplyr::filter(pathway2ko,grepl("map",Pathway))
    pathway2ko$KO=sub("^ko:","",pathway2ko$KO);pathway2ko$Pathway=sub("^path:","",pathway2ko$Pathway)
    pathway2ko_num=pathway2ko%>%dplyr::count(Pathway,name = "K_num")
    pathway_desc=utils::read.table(paste0(dd,"/pathway.desc.list"),sep = "\t",col.names =c("Pathway","Description"))
    pathway2ko_com=stats::aggregate(KO~Pathway,pathway2ko,paste, collapse = ",")
    pathway_list=Reduce(dplyr::left_join, list(pathway2ko_num,pathway2ko_com,pathway_desc))
    colnames(pathway_list)=c("id","K_num","KOs","Description")

    module2ko=utils::read.table(paste0(dd,"/module-KO.list"),sep = "\t",col.names =c("KO","Module"))
    module2ko$KO=sub("^ko:","",module2ko$KO);module2ko$Module=sub("^md:","",module2ko$Module)
    module2ko_num=module2ko%>%dplyr::count(Module,name = "K_num")
    module_desc=utils::read.table(paste0(dd,"/module.desc.list"),sep = "\t",col.names =c("Module","Description"))
    module2ko_com=stats::aggregate(KO~Module,module2ko,paste, collapse = ",")
    module_list=Reduce(dplyr::left_join, list(module2ko_num,module2ko_com,module_desc))
    colnames(module_list)=c("id","K_num","KOs","Description")

    KOlist=list("pathway"=pathway_list,"module"=module_list)

    attributes(KOlist)$download_time=file.info(paste0(dd,"/pathway-KO.list"))$mtime
    attributes(KOlist)$build_time=Sys.time()

    ko_desc=utils::read.table(paste0(dd,"//ko.desc.list"),sep = "\t",col.names =c("KO","Description"),quote = "")

    if(is.null(output))output=paste0(dd,"/KOlist.rda")
    save(KOlist,ko_desc,file =output)
}


#' Build a customize moduleist
#'
#' @param pathway2ko user input annotation of Pathway to KO mapping, a data.frame of 2 column with pathway and ko.
#' @param pathway2desc user input of Pathway TO Description mapping, a data.frame of 2 column with pathway and description.
#'
#' @return a customize modulelist
#' @export
custom_modulelist=function(pathway2ko,pathway2desc=NULL){
    Pathway=NULL
    pathway2ko=pathway2ko[,1:2]
    colnames(pathway2ko)=c("Pathway","KO")
    pathway2ko_num=pathway2ko%>%dplyr::count(Pathway,name = "K_num")

    if(is.null(pathway2desc)){
        pathway2desc=data.frame(Pathway=pathway2ko_num$Pathway,Description=pathway2ko_num$Pathway)
    }
    else pathway2desc=pathway2desc[,1:2]
    colnames(pathway2desc)=c("Pathway","Description")

    pathway2ko_com=stats::aggregate(KO~Pathway,pathway2ko,paste, collapse = ",")
    pathway_list=Reduce(dplyr::left_join, list(pathway2ko_num,pathway2ko_com,pathway2desc))
    colnames(pathway_list)=c("id","K_num","KOs","Description")
    pathway_list
}


#' Load the KOlist
#'
#' @param KOlist_file NULL or your `KOlist.rda` results from `make_KO_list`
#' @param envir `.GlobalEnv` (default) or `environment()`
#'
#' @export
#' @return KOlist in `.GlobalEnv`
load_KOlist=function(KOlist_file=NULL,envir=.GlobalEnv){
    if(is.null(KOlist_file)){
        #KOlist_file=system.file("data","new_KOlist.rda",package = "ReporterScore")
        KOlist_file=file.path(tools::R_user_dir("ReporterScore"),"new_KOlist.rda")
    }
    if(file.exists(KOlist_file))load(KOlist_file,envir = envir)
    else data("KOlist",package = "ReporterScore",envir = envir)
}
