#' Update the KO2Pathway and KO2Module files from KEGG
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
    if(!dir.exists(pack_dir))dir.create(pack_dir,recursive = TRUE)

    if(is.null(download_dir))dd="ReporterScore_temp_download"
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

#' Prepare KOlist from KEGG
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

#' Build a custom moduleist
#'
#' @param pathway2ko user input annotation of Pathway to KO mapping, a data.frame of 2 column with pathway and ko.
#' @param pathway2desc user input of Pathway TO Description mapping, a data.frame of 2 column with pathway and description.
#'
#' @return a custom modulelist
#' @export
#' @examples
#' mydat=data.frame(pathway=paste0("PATHWAY",rep(1:2,each=5)),ko=paste0("K",1:10))
#' mymodulelist=custom_modulelist(mydat)
#' mymodulelist
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

#' Load the KOlist (from KEGG)
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

#' Load the KO_htable (from KEGG)
#'
#' @param envir `.GlobalEnv` (default) or `environment()`
#'
#' @export
#' @return KO_htable in `.GlobalEnv`
load_KO_htable=function(envir=.GlobalEnv){
    KO_htable_file=file.path(tools::R_user_dir("ReporterScore"),"new_KO_htable.rda")
    if(file.exists(KO_htable_file))load(KO_htable_file,envir = envir)
    else data("KO_htable",package = "ReporterScore",envir = envir)
}

if(F){
    update_GO_file=function(){
        GO_list=clusterProfiler:::get_GOTERM()
        GO_list=GO_list[,-1]
    }

    gene_to_go=data.frame(gene=rep(c("dnaJ","hspR"),each=3),go=c("GO:0005575","GO:0005618","GO:0005623",
                                                                 "GO:0005618","GO:0005623","GO:0005506"))
    term2gene1 <- gene_to_go[, c(2, 1)]
    #为直接注释补充为间接注释
    term2gene <- clusterProfiler::buildGOmap(term2gene1)
    #将GoId转换为GoTerm
    go2term <- clusterProfiler::go2term(term2gene$GO)
    gene1=c("dnaJ")
    df <- enricher(gene = gene1, TERM2GENE = term2gene, TERM2NAME = go2term, pvalueCutoff = 1, qvalueCutoff = 1)
    #GO注释不是像KEGG一样，KO map到Pathway，而是直接gene map到 GOterm
}


#' update_KO_htable from KEGG
#'
#' @param file ko00001.keg from https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext
#' @param download_dir download_dir
#'
#' @return KO_htable
#' @export
#'
update_KO_htable=function(file=NULL,download_dir=NULL){
    pcutils::lib_ps("readr")
    #if(is.null(pack_dir))pack_dir=paste0(.libPaths()[1],"/ReporterScore")
    pack_dir=tools::R_user_dir("ReporterScore")
    if(!dir.exists(pack_dir))dir.create(pack_dir,recursive = TRUE)

    if(is.null(download_dir))dd="ReporterScore_temp_download"
    else dd=download_dir

    if(!dir.exists(dd))dir.create(dd)

    if(is.null(file)){
        pcutils::dabiao("Trying to download files from https://rest.kegg.jp/ ")
        utils::download.file("https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext",
                             destfile = paste0(dd,"/ko00001.keg"),method = "curl")
        file=paste0(dd,"/ko00001.keg")
    }
    else {
        if(!grepl("ko00001.keg",file))stop("file should be path of ko00001.keg.")
    }
    ko00001_htext2df(file,paste0(dd,"/ko00001.tsv"))
    readr::read_delim(paste0(dd,"/ko00001.tsv"),delim = "\t",
                      col_names = c("level1_id","level1_name","level2_id","level2_name",
                                    "level3_id","level3_name","KO_id","KO_name"))->KO_htable
    # readr::read_delim(paste0(dd,"/ko00002.tsv"),delim = "\t",
    #                   col_names = c("level1_name","level2_name",
    #                                 "level3_name","Module_id","Module_name"))->M_htable
    attributes(KO_htable)$download_time=file.info(file)$mtime
    attributes(KO_htable)$build_time=Sys.time()
    save(KO_htable,file = paste0(pack_dir,"/new_KO_htable.rda"))
    pcutils::dabiao(paste0("Update done at ",Sys.time()))
}

ko00001_htext2df=function(input_file, output_file,header_num = 0) {
    # Read the input file
    input <- readLines(input_file)

    # Open the output file
    output <- file(output_file, "w")

    count <- list()
    while (header_num > 0) {
        input <- input[-1]
        header_num <- header_num - 1
    }

    # Process the input file
    for (line in input) {
        if (grepl("^A", line)) {
            a <- line
            a <- sub(" ", "\t", a)
            next
        }
        if (grepl("^B  ", line)) {
            b <- line
            b <- sub("B  ", "", b)
            b <- sub(" ", "\t", b)
            next
        }
        if (grepl("^C    ", line)) {
            c <- line
            c <- sub("C    ", "", c)
            c <- sub(" ", "\t", c)
            next
        }
        if (grepl("^D      ", line)) {
            d <- line
            d <- sub("D      ", "", d)
            d <- sub("  ", "\t", d)
            writeLines(paste(a, b, c, d, sep = "\t"), output)
        }
    }

    # Close the output file
    close(output)
}

ko00002_htext2df=function(input_file, output_file,header_num = 0) {
    # Read the input file
    input <- readLines(input_file)

    # Open the output file
    output <- file(output_file, "w")

    count <- list()
    while (header_num > 0) {
        input <- input[-1]
        header_num <- header_num - 1
    }

    # Process the input file
    for (line in input) {
        if (grepl("^A", line)) {
            a <- line
            a <- sub("A", "", a)
            a <- gsub("</?b>", "", a)
            next
        }
        if (grepl("^B  ", line)) {
            b <- line
            b <- sub("B  ", "", b)
            b <- gsub("</?b>", "", b)
            next
        }
        if (grepl("^C    ", line)) {
            c <- line
            c <- sub("C    ", "", c)
            next
        }
        if (grepl("^D      ", line)) {
            d <- line
            d <- sub("D      ", "", d)
            d <- sub("  ", "\t", d)
            writeLines(paste(a, b, c, d, sep = "\t"), output)
        }
    }

    # Close the output file
    close(output)
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
