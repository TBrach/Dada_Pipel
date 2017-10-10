### wilcoxon version median

```{r, warning = FALSE, message = FALSE, echo = FALSE}
# TbTresult_list_median <- evaluate_TbTmatrixes_wilcoxTest(TbTmatrixes_list, physeq = ps_filt, group_var = group_var, p.adjust.method = "fdr", type = "median")
# measure_matrixes <- lapply(TbTresult_list_median, `[[`, 2)
# TbTresult_list_median <- lapply(TbTresult_list_median, `[[`, 1)
# TbTresult_list_median_show <- lapply(TbTresult_list_median, function(df) {
#         select(df, Family, Genus, Species, 2:19, 22)
# })
```

```{r, warning = FALSE, message = FALSE, echo = FALSE, results = "asis"}
# suppressWarnings(head_values <- sapply(TbTresult_list_median, function(df){
#         max(which(df$p_val_adj < 0.05))
# }))
# head_values <- head_values + 5
# head_values[head_values == -Inf] <- 10
# head_values[head_values > 35] <- 35
# for (i in 1:length(TbTresult_list_median_show)){
#         print(knitr::kable(head(TbTresult_list_median_show[[i]], head_values[i]), caption = names(TbTresult_list_median_show)[i]))
# }
```

### wilcoxon version sum

```{r, warning = FALSE, message = FALSE, echo = FALSE}
# TbTresult_list_sum <- evaluate_TbTmatrixes_wilcoxTest(TbTmatrixes_list, physeq = ps_filt, group_var = group_var, p.adjust.method = "fdr", type = "sum")
# measure_matrixes <- lapply(TbTresult_list_median, `[[`, 2)
# TbTresult_list_median <- lapply(TbTresult_list_median, `[[`, 1)
# TbTresult_list_sum_show <- lapply(TbTresult_list_sum, function(df) {
#         select(df, Family, Genus, Species, 2:19, 22)
# })
# ```
# 
# ```{r, warning = FALSE, message = FALSE, echo = FALSE, results = "asis"}
# suppressWarnings(head_values <- sapply(TbTresult_list_sum, function(df){
#         max(which(df$p_val_adj < 0.05))
# }))
# head_values <- head_values + 5
# head_values[head_values == -Inf] <- 10
# head_values[head_values > 35] <- 35
# for (i in 1:length(TbTresult_list_sum_show)){
#         print(knitr::kable(head(TbTresult_list_sum_show[[i]], head_values[i]), caption = names(TbTresult_list_sum_show)[i]))
# }
```

### taxon vs taxon tile plots TAKES TIME, better only for few taxa

```{r, warning = FALSE, message = FALSE, echo = FALSE}
# TbT_tiles <- create_TbT_TilePlots(TbTmatrixes_list, physeq = ps_filt, group_var = group_var)
# TbT_tiles_show <- lapply(TbT_tiles, `[[`, 2)
```

```{r, warning = FALSE, message = FALSE, echo = FALSE, results = "asis"}
# do.call("grid.arrange", c(TbT_tiles_show, ncol = 1))
```