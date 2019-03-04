
<!-- README.md is generated from README.Rmd. Please edit that file -->
parallelGO <img src="man/figures/logo.png" align="right" height=140/>
=====================================================================

The goal of `parallelGO` is to provide parallel implementation (using [foreach](https://cran.r-project.org/web/packages/foreach/) package) for gene ontology (GO) enrichment analysis given labeled sets of genes. The project was inspired by the need to do a fast evaluation of biological significance of identified clusters (also known as modules or sets) within gene networks.

Installation
------------

You can install the released version of parallelGO from [github](https://github.com) with:

``` r
install.packages("remotes")
remotes::install_github("vitalinakomashko/parallelGO")
```

This package requires two Bioconductor packages: [AnnotationDbi](https://bioconductor.org/packages/release/bioc/html/AnnotationDbi.html) and [GOstats](https://bioconductor.org/packages/release/bioc/html/GOstats.html). If you don't have these packages installed you can either install them by following instructions on the pages for these packages or, prior to installation of parallelGO from github, run this line in R:

``` r
setRepositories(ind = 1:3)
```

This will set the repositories to CRAN, BioC software and BioC annotation. This idea was found [here](https://stackoverflow.com/a/20479243/1655368).

How to use
----------

If you have a file with human gene symbols where each symbol is assigned to a gene set you can run the analysis by invoking the following code:

``` r
library(parallelGO)
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

In addition, we provide a function for filtering p-values (raw or adjusted) which can be called after either scenario described above:

``` r
res_filtered <- filter_pvalues(res, cutoff = 0.05)
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

Benchmarking using the example dataset
--------------------------------------

We benchmarked performance using the code above and the sample dataset.

macOS: MacBook Pro 2.7 GHz Intel Core i7, 16 Gb 2133 MHz LPDDR3.

Serial execution:

``` r
system.time(res <- run_parallel_go(dat_large_sets, 
                                   species = "human", 
                                   universe = universe, 
                                   run_parallel = FALSE))
```

    Parameter run_parallel is FALSE. Computation will be run sequentially.
       user  system elapsed 
    236.094  46.896 284.544

Parallel execution using 4 cores:

``` r
system.time(res <- run_parallel_go(dat_large_sets, 
                                   species = "human", 
                                   universe = universe))
```

    Parameter 'cores' is not provided. Getting the number of available cores
    with foreach::getDoParWorkers(). GO enrichment will be run in parallel on
    4 cores using doParallelMC backend.
       user  system elapsed 
    243.602  55.488 108.968
