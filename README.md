
# TAPE

<!-- badges: start -->
<!-- badges: end -->

The goal of TAPE is to ...

## Installation

You can install the development version of TAPE from [Github](https://github.com/styvon/TAPE) with:

``` r
devtools::install_github("styvon/TAPE",ref="main")
```

## Example


``` r
library(TAPE)
## input
file_cova = system.file("extdata", "example_cova.txt", package = "TAPE")
file_kmat = system.file("extdata", "example_kmat.rds", package = "TAPE")
file_geno_grm = tools::file_path_sans_ext(system.file("extdata", "example_geno_grm.bim", package = "TAPE"))
file_geno_test = tools::file_path_sans_ext(system.file("extdata", "example_geno_test.bgen", package = "TAPE"))
file_bgi_test = system.file("extdata", "example_geno_test.bgi", package = "TAPE")
file_idsingeno = system.file("extdata", "example_idsingeno.txt", package = "TAPE")
file_output_s2 = paste0(dirname(file_cova),"/output_s2.txt") # change this to your path for output file

## step 1
data = read.table(file_cova,header=T)
kmat = readRDS(file_kmat)

obj_null <- TAPE_Null_Model(y ~ sex+age, data = data, K=kmat, KgenFile=file_geno_grm, idstoIncludeFile=file_idsingeno, tau=rep(0,3),fixtau=rep(0,3),verbose=T)

## step 2
n_variants_tested = TAPEtest(null_object=obj_null, genfile=file_geno_test, samplefile=file_idsingeno, outfile=file_output_s2, genfile_format="bgen", bgi_file="1")
```

