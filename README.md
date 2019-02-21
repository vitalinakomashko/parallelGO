
<!-- README.md is generated from README.Rmd. Please edit that file -->
parallelGO <img src="man/figures/logo.png" align="right" height=140/>
=====================================================================

The goal of `parallelGO` is to provide parallel implementation (using [foreach](https://cran.r-project.org/web/packages/foreach/) package) for gene ontology (GO) enrichment analysis given labeled sets of genes. The project was inspired by the need to do a fast evaluation of biological significance of identified clusters (also known as modules or sets) within gene networks.

Installation
------------

You can install the released version of parallelGO from [github](https://github.com) with:

``` r
remotes::install_github("vitalinakomashko/parallelGO")
```

How to use
----------

If you have a file with human gene symbols where each symbol is assigned to a gene set you can run the analysis by invoking the following code:

``` r
out <- wrap_and_go("path_to_my_file", col_names = TRUE, delim = ",", 
                   id = "symbol", species = "human")
```

Mouse species and ensembl identifiers are also supported.

**IMPORTANT** Please prepare your file so that gene identifiers are in the first column and set labels are in the second column.

`wrap_and_go` will read the file, perform some minor cleaning, map identifiers to ENTREZ gene identifiers and run GO enrichment analysis in parallel. Please read documentation for `wrap_and_go` to learn about additional parameters.

If you don't have a file, but rather a data frame (`dat`) in memory, you can perform your analyses step by step in the same manner as it as performed by `wrap_and_go`:

``` r
# remove duplicated rows if present:
dat_clean <- deduplicate_rows(dat)
# map to ENTREZ gene identifiers:
dat_mapped <- map_genes(dat_clean, id = "symbol", species = "human")
# extract the universe for hypergeometric test:
universe <- unique(dat_mapped$entrez)
# remove sets with a small number of genes, OPTIONAL:
dat_large_sets <- remove_small_sets(dat_mapped, min_set_size = 30)
# run GO analysis in parallel
res <- run_parallel_go(dat_large_sets, species = "human", universe = universe)
```

Example dataset
---------------

The package also provides a small example dataset for illustration purposes

``` r
data("human_symbol")
head(human_symbol)
#>        id set_label
#> 1   STPG1        24
#> 2 CACNA1G        23
#> 3  MAP3K9        24
#> 4    MDH1       156
#> 5   SYT13       159
#> 6  GABRA1        24
str(human_symbol)
#> 'data.frame':    280 obs. of  2 variables:
#>  $ id       : chr  "STPG1" "CACNA1G" "MAP3K9" "MDH1" ...
#>  $ set_label: int  24 23 24 156 159 24 158 157 158 156 ...
```
