# ReporterScore v0.1.8 Notes

## Fixed

- fixed the filepath problem ':' in Windows, issue [#6](https://github.com/Asa12138/ReporterScore/issues/6). <2024-09-14, Sat>

# ReporterScore v0.1.7 Notes

## Added

- add README_zh_CN.md <2024-08-20, Tue>
- add `export_report_table` function <2024-08-20, Tue>

# ReporterScore v0.1.6 Notes

## Fixed

- fixed the `ko.test` for one sample group <2024-06-25, Tue>

## Added

- added `kos_description` argument for `plot_features_network()` <2024-06-25, Tue>

# ReporterScore v0.1.5 Notes

## Added

- added citation information <2024-06-03, Mon>

# ReporterScore v0.1.4 Notes

## Fixed

- fixed `up_level_KO` as level3_name would be duplicated <2024-04-10, Wed>
- fixed `load_CPDlist()` as no return <2024-04-08, Mon>
- fixed `get_modulelist()` as the third argument is `gene` <2024-04-07, Sun>

## Added

- added `combine_rs_res()` for combining the results of 'step by step GRSA' to `reporterscore` object <2024-04-08, Mon>
- method can be 'none' in `ko.test()`, and return a NA pvalue <2024-04-08, Mon>

# ReporterScore v0.1.3 Notes

## Added

- add the `KO_gsva()`, `KO_sea()`, `KO_safe()` and `KO_padog()`. <2024-01-19, Fri>

# ReporterScore v0.1.2 Notes

## Added

- add the `plot_KOs_network()` as `MetaNet` is available on CRAN <2024-01-15>

# ReporterScore v0.1.1 Notes

## Fixed

- revise and re-submit to CRAN <2024-01-12>
- available on CRAN <https://CRAN.R-project.org/package=ReporterScore> <2024-01-12>

# ReporterScore v0.1.0 Notes

## Fixed

- fixed R CMD check and Bioc check <2023-12-29>
- submit to CRAN <2024-01-10, Wed>

