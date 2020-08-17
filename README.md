# TMEscore
--------------------------------------


TME infiltration patterns were determined and systematically correlated with TME cell phenotypes, genomic traits, and patient clinicopathological features to establish the "TMEscore" [1]

TMEscore is a R package to do Tumor microenvironment analysis. Main advantages:
- Provides functionality to calculate Tumor microenvironment (TME) score (PCA or z-score)
- Functions to visulize TME data.
- Identify TME relevant mutaions.

The package is not yet on CRAN. You can install from Github:

``` r
# install.packages("devtools")
devtools::install_github("DongqiangZeng0808/TMEscore")
```
Main documentation is on the `tmescore` function in the package:

``` r
library('TMEscore')
help('tmescore')
```

Example

``` r
tmescore<-tmescore(eset = eset_stad, pdata = pdata_stad, colum_of_sample = "ID",classify = T)
head(tmescore)
              ID subtype   time status TMEscoreA TMEscoreB  TMEscore TMEscore_binary
284 TCGA-RD-A8N2    <NA> 118.00      0 -7.306563  13.24346 -20.55003             Low
95  TCGA-BR-A4IV      GS  28.97      1 -6.743132  12.61978 -19.36292             Low
66  TCGA-BR-8371      GS  11.97      1 -7.024702  12.56123 -19.58593             Low
69  TCGA-BR-8380      GS     NA      1 -5.855567  12.97473 -18.83030             Low
101 TCGA-BR-A4J9      GS   0.47      0 -6.643521  11.98279 -18.62631             Low
82  TCGA-BR-8592      GS   6.37      1 -5.143173  12.47427 -17.61744             Low
```

References
----------
Manuscript is under review.

E-mail any questions to dongqiangzeng0808@gmail.com
