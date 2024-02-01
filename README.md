
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![](https://img.shields.io/badge/doi-10.1101/2023.10.13.562235-yellow.svg)](https://doi.org/10.1101/2023.10.13.562235)
[![](https://img.shields.io/badge/blog-@asa-blue.svg)](https://asa-blog.netlify.app/)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ReporterScore)](https://cran.r-project.org/package=ReporterScore)
[![](http://cranlogs.r-pkg.org/badges/last-month/ReporterScore)](https://cran.r-project.org/package=ReporterScore)
[![](https://www.r-pkg.org/badges/version/ReporterScore?color=green)](https://cran.r-project.org/package=ReporterScore)
[![](https://img.shields.io/badge/devel%20version-0.1.3-green.svg)](https://github.com/Asa12138/ReporterScore)
<!-- badges: end -->

# ReporterScore

Generalized Reporter Score-based Enrichment Analysis for Omics Data

<img src="README_files/1-workflow.png" width="1985" />

## Citation

To cite ReporterScore in publications use:

C. Peng, Q. Chen, S. Tan, X. Shen, C. Jiang, Generalized Reporter
Score-based Enrichment Analysis for Omics Data. bioRxiv \[Preprint\]
(2023). <https://doi.org/10.1101/2023.10.13.562235>.

## Install

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Asa12138/pcutils")
devtools::install_github("Asa12138/ReporterScore")

library(ReporterScore)
```

## Usage

### 1. Inputdata (KO abundance table and metadata)

For transcriptomic, scRNA-seq, and related gene-based omics data of a
specific species, a complete gene abundance table can be used. For
metagenomic and metatranscriptomic data, which involve many different
species, a KO abundance table can be used, generated using Blast,
Diamond, or KEGG official mapper software to align the reads or contigs
to the KEGG or the EggNOG database For metabolomic data, an annotated
compound abundance table can be used, but the standardization of
compound IDs (e.g., convert compound IDs to C numbers in the KEGG
database) is required.

An example code tailored for a KO abundance table is as follows:

``` r
data("KO_abundance_test")
head(KO_abundance[, 1:6])
##                WT1         WT2         WT3         WT4         WT5         WT6
## K03169 0.002653545 0.005096380 0.002033923 0.000722349 0.003468322 0.001483028
## K07133 0.000308237 0.000280458 0.000596527 0.000859854 0.000308719 0.000878098
## K03088 0.002147068 0.002030742 0.003797459 0.004161979 0.002076596 0.003091182
## K03530 0.003788366 0.000239298 0.000445817 0.000557271 0.000222969 0.000529624
## K06147 0.000785654 0.001213630 0.001312569 0.001662740 0.002387006 0.001725797
## K05349 0.001816325 0.002813642 0.003274701 0.001089906 0.002371921 0.001795214
```

And you should also offer a experimental metadata (rows are samples,
columns are groups).

``` r
head(metadata)
##     Group Group2
## WT1    WT     G3
## WT2    WT     G3
## WT3    WT     G3
## WT4    WT     G3
## WT5    WT     G3
## WT6    WT     G1
```

### 2. Pathway database

The `ReporterScore` package has built-in KEGG pathway, module, gene,
compound, and GO databases and also allows customized databases, making
it compatible with feature abundance tables from diverse omics data.

1.  `ReporterScore` has built-in KEGG pathway-KO and module-KO databases
    (2023-08 version) for KO abundance table. You can use
    `load_KOlist()` to have a look and use `update_KO_file()` to update
    these databases (by KEGG API) as using latest database is very
    important.

2.  `ReporterScore` has built-in KEGG pathway-compound and
    module-compound databases (2023-08 version) for compound abundance
    table. You can use `load_CPDlist()` to have a look and use
    `update_KO_file()` to update these databases (by KEGG API).

3.  `ReporterScore` has built-in pathway-ko, pathway-gene, and
    pathway-compound databases of human (hsa) and mouse (mmu) for
    ko/gene/compound abundance table. You can use
    `custom_modulelist_from_org()` to have a look. Use
    `update_org_pathway()` to update these databases and download other
    organism databases (by KEGG API).

4.  `ReporterScore` has built-in GO-gene database, You can use
    `load_GOlist()` to have a look and use `update_GOlist()` to update
    these databases (by KEGG API).

5.  You can just customize your own pathway databases (gene set of
    interest) by using `custom_modulelist()`.

``` r
# 1. KEGG pathway-KO and module-KO databases
KOlist <- load_KOlist()
head(KOlist$pathway)

# 2. KEGG pathway-compound and module-compound databases
CPDlist <- load_CPDlist()
head(CPDlist$pathway)

# 3. human (hsa) pathway-ko/gene/compound databases
hsa_pathway_gene <- custom_modulelist_from_org(
    org = "hsa",
    feature = c("ko", "gene", "compound")[2]
)
head(hsa_pathway_gene)

# 4. GO-gene database
GOlist <- load_GOlist()
head(GOlist$BP)

# 5. customize your own pathway databases
?custom_modulelist()
```

### 3. One step enrichment

Use function `reporter_score` can get the reporter score result by one
step.

there are some important arguments for analysis:

- **mode**: “mixed” or “directed” (only for two groups differential
  analysis or multi-groups correlation analysis.), see details in
  `pvalue2zs`.
- **method**: the type of test. Default is `wilcox.test`:
  - `t.test` (parametric) and `wilcox.test` (non-parametric). Perform
    comparison between two groups of samples. If the grouping variable
    contains more than two levels, then a pairwise comparison is
    performed.
  - `anova` (parametric) and `kruskal.test` (non-parametric). Perform
    one-way ANOVA test comparing multiple groups.
  - “pearson”, “kendall”, or “spearman” (correlation), see `cor`.
- **type**/**modulelist**: choose the pathway database, “pathway”,
  “module”, for default database, or use a customized modulelist.
- **feature**: one of “ko”, “gene”, “compound”.

The first level will be set as the **control group**, you can change the
factor level to change your comparison.

For example, we want to compare two groups ‘WT-OE’, and use the
“directed” mode as we just want know which pathways are enriched or
depleted in **OE group**:

#### KO-pathway

``` r
cat("Comparison: ", levels(factor(metadata$Group)))
## Comparison:  WT OE

# for microbiome!!!
reporter_res <- reporter_score(KO_abundance, "Group", metadata,
    mode = "directed",
    method = "wilcox.test", perm = 999
)
## ================================Use feature: ko=================================
## ===============================Checking rownames================================
## Some of your ko_stat are not KO id, check the format! (e.g. K00001)
## 52.7% of your kos in the modulelist!
## 30 samples are matched for next step.
## ===========================Removing all-zero rows: 0============================
## ===================================1.KO test====================================
## =================================Checking group=================================
## 30 samples are matched for next step.
## ===========================Removing all-zero rows: 0============================
## ==============================Calculating each KO===============================
## ===========================Using method: wilcox.test============================
## 1000 features done.
## 2000 features done.
## 3000 features done.
## 4000 features done.
## 
## Compared groups: WT, OE
## Total KO number: 4535
## Compare method: wilcox.test
## Time use: 1.177
## =========================2.Transfer p.value to z-score==========================
## ==========================3.Calculating reporter score==========================
## ==================================load KOlist===================================
## ===================KOlist download time: 2023-08-14 16:00:52====================
## If you want to update KOlist, use `update_KO_file()`
## ============================Calculating each pathway============================
## 100 pathways done.
## 400 pathways done.
## ID number: 481
## Time use: 1.604
## ====================================All done====================================
```

The result is a “reporter_score” object:

| elements     | description                                       |
|--------------|---------------------------------------------------|
| `kodf`       | your input KO_abundance table                     |
| `ko_stat`    | ko statistics result contains p.value and z_score |
| `reporter_s` | the reporter score in each pathway                |
| `modulelist` | default KOlist or customized modulelist dataframe |
| `group`      | The comparison groups in your data                |
| `metadata`   | sample information dataframe contains group       |

#### Gene-pathway

When you use the gene abundance table of a specific species
(e.g. human), remember to set the `feature` and `type`!!! Or give the
database through `modulelist`.

``` r
data("genedf")

# Set the `feature` and `type`!
reporter_res_gene <- reporter_score(genedf, "Group", metadata,
    mode = "directed",
    feature = "gene", type = "hsa",
    method = "wilcox.test", perm = 999
)

# Or give the database through `modulelist`
hsa_pathway_gene <- custom_modulelist_from_org(
    org = "hsa",
    feature = "gene"
)

reporter_res_gene <- reporter_score(genedf, "Group", metadata,
    mode = "directed",
    modulelist = hsa_pathway_gene,
    method = "wilcox.test", perm = 999
)
```

#### Compound-pathway

``` r
reporter_res_gene <- reporter_score(chem_df, "Group", metadata,
    mode = "directed",
    feature = "compound", type = "hsa",
    method = "wilcox.test", perm = 999
)
```

### 4. Visualization

Plot the most significantly enriched pathways:

``` r
# View(reporter_res$reporter_s)
plot_report_bar(reporter_res, rs_threshold = c(-2.5, 2.5), facet_level = TRUE)
## ==============================load Pathway_htable===============================
## ===============Pathway_htable download time: 2024-01-12 00:52:39================
## If you want to update Pathway_htable, use `update_htable(type='pathway')`
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

Plot the most significantly enriched pathways (circle packing):

``` r
plot_report_circle_packing(reporter_res, rs_threshold = c(-2.5, 2.5))
## ==============================load Pathway_htable===============================
## ===============Pathway_htable download time: 2024-01-12 00:52:39================
## If you want to update Pathway_htable, use `update_htable(type='pathway')`
## Non-leaf weights ignored
## Scale for fill is already present.
## Adding another scale for fill, which will replace the existing scale.
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

When we focus on one pathway, e.g. “map00780”:

``` r
plot_KOs_in_pathway(reporter_res, map_id = "map00780")
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Or show the distribution of Z-scores

``` r
plot_KOs_distribution(reporter_res, map_id = "map00780")
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Or display as a network:

``` r
plot_KOs_network(reporter_res,
    map_id = c("map00780", "map00785", "map00900"),
    main = "", mark_module = TRUE
)
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

And we also look at the KOs abundance in a pathway:

``` r
plot_KOs_box(reporter_res, map_id = "map00780", only_sig = TRUE)
## `geom_smooth()` using formula = 'y ~ x'
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Or display as a heatmap:

``` r
plot_KOs_heatmap(reporter_res,
    map_id = "map00780", only_sig = TRUE,
    heatmap_param = list(cutree_rows = 2)
)
```

![](README_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Or plot the KEGG pathway:

``` r
plot_KEGG_map(reporter_res$ko_stat, map_id = "map00780", color_var = "Z_score")
```

<img src="README_files/ko00780.Z_score.png" width="586" />

### Example for multi-group or longitudinal

If our experimental design is more than two groups or longitudinal, we
can choose multi-groups comparison (or correlation):

``` r
cat("Comparison: ", levels(factor(metadata$Group2)))
## Comparison:  G1 G2 G3

reporter_res2 <- reporter_score(KO_abundance, "Group2", metadata,
    mode = "directed",
    method = "spearman", p.adjust.method1 = "none", perm = 999
)
## ================================Use feature: ko=================================
## ===============================Checking rownames================================
## Some of your ko_stat are not KO id, check the format! (e.g. K00001)
## 52.7% of your kos in the modulelist!
## 30 samples are matched for next step.
## ===========================Removing all-zero rows: 0============================
## ===================================1.KO test====================================
## =================================Checking group=================================
## 30 samples are matched for next step.
## ===========================Removing all-zero rows: 0============================
## ==============================Calculating each KO===============================
## =============================Using method: spearman=============================
## Using correlation analysis: spearman, the groups will be transform to numeric, note the factor feature of group.
## 1000 features done.
## 2000 features done.
## 3000 features done.
## 4000 features done.
## 
## Compared groups: G1, G2, G3
## Total KO number: 4535
## Compare method: spearman
## Time use: 0.510
## =========================2.Transfer p.value to z-score==========================
## ==========================3.Calculating reporter score==========================
## ==================================load KOlist===================================
## ===================KOlist download time: 2023-08-14 16:00:52====================
## If you want to update KOlist, use `update_KO_file()`
## ============================Calculating each pathway============================
## 100 pathways done.
## 400 pathways done.
## ID number: 481
## Time use: 1.578
## ====================================All done====================================

plot_KOs_in_pathway(reporter_res2, map_id = "map02060") + scale_y_log10()
```

![](README_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Example for specified pattern

For example, groups “G1”, “G2”, and “G3” can be set as 1, 10, and 100 if
an exponentially increasing trend is expected.

We use 1,5,1 to found pathways with the down-up-down pattern

``` r
reporter_res3 <- reporter_score(KO_abundance, "Group2", metadata,
    mode = "directed", perm = 999,
    method = "pearson", pattern = c("G1" = 1, "G2" = 5, "G3" = 1)
)
plot_report_bar(reporter_res3, rs_threshold = 3, show_ID = TRUE)
```

![](README_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
plot_KOs_in_pathway(reporter_res3, map_id = "map00860")
```

![](README_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

To explore potential patterns within the data, clustering methods, such
as C-means clustering, can be used.

``` r
rsa_cm_res <- RSA_by_cm(KO_abundance, "Group2", metadata,
    method = "pearson",
    k_num = 3, perm = 999
)
# show the patterns
plot_c_means(rsa_cm_res, filter_membership = 0.7)
```

![](README_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r

plot_report_bar(rsa_cm_res, rs_threshold = 2.5, y_text_size = 10)
```

![](README_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

## Details

### Step by step

The one step function `reporter_score` consists of three parts：

1.  `ko.test`: this function help to calculate *p-value* for
    KO_abundance by various built-in methods such as differential
    analysis (t.test, wilcox.test, kruskal.test, anova) or correlation
    analysis (pearson, spearman, kendall). **You can also calculate this
    *p-value* for KO_abundance by other methods** like “DESeq2”,
    “Edger”, “Limma”, “ALDEX”, “ANCOM” and do a p.adjust yourself, then
    skip `ko.test` step go to step2…
2.  `pvalue2zs`: this function transfers p-value of KOs to Z-score
    (select mode: “mixed” or “directed”).
3.  `get_reporter_score` this function calculate reporter score of each
    pathways in a specific database. You can use a custom database here.

So that you can get reporter score step by step.

### Other commonly used enrichment methods

`ReporterScore` also provides other enrichment methods like
`KO_fisher`(fisher.test), `KO_enrich`(fisher.test, from
`clusterProfiler`), `KO_gsea` (GSEA, from `clusterProfiler`), The input
data is from `reporter_score`, and also supports custom databases, so
you can easily compare the results of various enrichment methods and
conduct a comprehensive analysis:

``` r
# View(reporter_res2$reporter_s)
# reporter_score
filter(reporter_res$reporter_s, abs(ReporterScore) > 1.64, p.adjust < 0.05) %>% pull(ID) -> RS
# fisher
fisher_res <- KO_fisher(reporter_res)
filter(fisher_res, p.adjust < 0.05) %>% pull(ID) -> Fisher
# enricher
enrich_res <- KO_enrich(reporter_res)
filter(enrich_res, p.adjust < 0.05) %>% pull(ID) -> clusterProfiler
# GESA
set.seed(1234)
gsea_res <- KO_gsea(reporter_res, weight = "Z_score")
filter(data.frame(gsea_res), p.adjust < 0.05) %>% pull(ID) -> GSEA

venn_res <- list(GRSA = RS, Fisher = Fisher, CP = clusterProfiler, GSEA = GSEA)
library(pcutils)
venn(venn_res)
```

![](README_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
venn(venn_res, "network")
```

![](README_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

## Other features

### uplevel the KOs

[KEGG BRITE](https://www.genome.jp/kegg/brite.html) is a collection of
hierarchical classification systems capturing functional hierarchies of
various biological objects, especially those represented as KEGG
objects.

We collected k00001 KEGG Orthology (KO) table so that you can summaries
each levels abundance. Use `load_KO_htable` to get KO_htable and use
`update_KO_htable` to update. Use `up_level_KO` can upgrade to specific
level in one of “pathway”, “module”, “level1”, “level2”, “level3”,
“module1”, “module2”, “module3”.

``` r
KO_htable <- load_KO_htable()
## =================================load KO_htable=================================
## ==================KO_htable download time: 2024-01-12 00:49:03==================
## If you want to update KO_htable, use `update_htable(type='ko')`
head(KO_htable)
##   level1_name             level2_name level3_id                  level3_name
## 1  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
## 2  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
## 3  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
## 4  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
## 5  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
## 6  Metabolism Carbohydrate metabolism  map00010 Glycolysis / Gluconeogenesis
##    KO_id                                                    KO_name
## 1 K00844                                HK; hexokinase [EC:2.7.1.1]
## 2 K12407                              GCK; glucokinase [EC:2.7.1.2]
## 3 K00845                              glk; glucokinase [EC:2.7.1.2]
## 4 K25026                              glk; glucokinase [EC:2.7.1.2]
## 5 K01810       GPI, pgi; glucose-6-phosphate isomerase [EC:5.3.1.9]
## 6 K06859 pgi1; glucose-6-phosphate isomerase, archaeal [EC:5.3.1.9]
plot_htable(type = "ko")
## =================================load KO_htable=================================
## ==================KO_htable download time: 2024-01-12 00:49:03==================
## If you want to update KO_htable, use `update_htable(type='ko')`
```

![](README_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
KO_level1 <- up_level_KO(KO_abundance, level = "level1", show_name = TRUE)
## =================================load KO_htable=================================
## ==================KO_htable download time: 2024-01-12 00:49:03==================
## If you want to update KO_htable, use `update_htable(type='ko')`
pcutils::stackplot(KO_level1[-which(rownames(KO_level1) == "Unknown"), ]) +
    ggsci::scale_fill_d3() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

### CARD for ARGs

For convenience, I also included the CARD database from
<https://card.mcmaster.ca/download/0/broadstreet-v3.2.8.tar.bz2>.

``` r
CARDinfo <- load_CARDinfo()
## =================================load CARDinfo==================================
## ==================CARDinfo download time: 2024-01-12 01:12:11===================
## If you want to update CARDinfo, use `update_GOlist()`
head(CARDinfo$ARO_index)
##         ARO Accession CVTERM ID Model Sequence ID Model ID
## 3005099   ARO:3005099     43314              6143     3831
## 3002523   ARO:3002523     38923              8144     1781
## 3002524   ARO:3002524     38924                85      746
## 3002525   ARO:3002525     38925              4719     1246
## 3002526   ARO:3002526     38926               228     1415
## 3002527   ARO:3002527     38927              5510     2832
##                                                     Model Name
## 3005099 23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)
## 3002523                                             AAC(2')-Ia
## 3002524                                             AAC(2')-Ib
## 3002525                                             AAC(2')-Ic
## 3002526                                             AAC(2')-Id
## 3002527                                             AAC(2')-Ie
##                                                       ARO Name
## 3005099 23S rRNA (adenine(2058)-N(6))-methyltransferase Erm(A)
## 3002523                                             AAC(2')-Ia
## 3002524                                             AAC(2')-Ib
## 3002525                                             AAC(2')-Ic
## 3002526                                             AAC(2')-Id
## 3002527                                             AAC(2')-Ie
##         Protein Accession DNA Accession                         AMR Gene Family
## 3005099        AAB60941.1    AF002716.1 Erm 23S ribosomal RNA methyltransferase
## 3002523        AAA03550.1      L06156.2                                 AAC(2')
## 3002524        AAC44793.1      U41471.1                                 AAC(2')
## 3002525        CCP42991.1    AL123456.3                                 AAC(2')
## 3002526        AAB41701.1      U72743.1                                 AAC(2')
## 3002527        CAC32082.1    AL583926.1                                 AAC(2')
##                                                                   Drug Class
## 3005099 lincosamide antibiotic;macrolide antibiotic;streptogramin antibiotic
## 3002523                                            aminoglycoside antibiotic
## 3002524                                            aminoglycoside antibiotic
## 3002525                                            aminoglycoside antibiotic
## 3002526                                            aminoglycoside antibiotic
## 3002527                                            aminoglycoside antibiotic
##                 Resistance Mechanism CARD Short Name
## 3005099 antibiotic target alteration  Spyo_ErmA_MLSb
## 3002523      antibiotic inactivation      AAC(2')-Ia
## 3002524      antibiotic inactivation      AAC(2')-Ib
## 3002525      antibiotic inactivation      AAC(2')-Ic
## 3002526      antibiotic inactivation      AAC(2')-Id
## 3002527      antibiotic inactivation      AAC(2')-Ie
```

# Reference

1.  Patil, K. R. & Nielsen, J. Uncovering transcriptional regulation of
    metabolism by using metabolic network topology. Proc Natl Acad Sci U
    S A 102, 2685–2689 (2005).

2.  L. Liu, R. Zhu, D. Wu, Misuse of reporter score in microbial
    enrichment analysis. iMeta. 2, e95 (2023).

3.  <https://github.com/wangpeng407/ReporterScore>
