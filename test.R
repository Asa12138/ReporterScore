#使用iMeta文章中的数据看看结果

data=read.csv("~/Documents/git/ReporterScore-main/demo/dataset/Demo_dataset.tsv",sep = "\t",row.names = 1)
data=data[,-ncol(data)]
Group_tab=data.frame(row.names = colnames(data),Group=gsub("_.*","",colnames(data)))

ko_pvalue=ko_test(data,"Group",Group_tab)
ko_stat=pvalue2zs(ko_pvalue,mode="directed")
#ko_stat=pvalue2zs(ko_pvalue,mode="mixed")
reporter_s=get_reporter_score(ko_stat)
plot_KOs_in_pathway("map01230",ko_stat)

ko_pvalue$q.value=p.adjust(ko_pvalue$p.value,method = "BH")
list(a=filter(ko_pvalue,q.value<0.05)%>%pull(KO_id),b=diffdata$ID)%>%pcutils::venn()

report_feature=read.csv('~/Documents/git/ReporterScore-main/demo/results/reporterfeature.csv')

#clusterProfiler
library('clusterProfiler')
diffdata <- read.csv('~/Documents/git/ReporterScore-main/demo/results/differential KOs.tsv', sep='\t',row.names = 1)
diffdata <- diffdata[which(diffdata$P.Value<0.05),]
head(diffdata)
dim(diffdata)

path2name <- read.csv('~/Documents/git/ReporterScore-main/demo/dataset/KEGG_pathways.csv',row.names = 1)
head(path2name)
path2name <- path2name[, c('ID', 'Name')]
path2ko <- read.csv('~/Documents/git/ReporterScore-main/demo/dataset/KEGG_path2ko.csv')
path2ko <- path2ko[, c('Pathway', 'KO')]

e <- enricher(diffdata$ID,TERM2GENE = path2ko,
              pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.05)

barplot(e, showCategory = 30)


