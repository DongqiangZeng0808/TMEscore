TMEscore
=====================================

TME infiltration patterns were determined and systematically correlated with TME cell phenotypes, genomic traits, and patient clinicopathological features to establish the "TMEscore" [1]

TMEscore is a R package to do Tumor microenvironment analysis. Main advantages:
- Provides functionality to calculate Tumor microenvironment (TME) score (PCA, ssGSEA, and z-score)
- Functions to visulize TME data.
- Deconvolution of cellular composition of tumor.
- Identify TME relevant mutaions.

The package is not yet on CRAN.
You can install from Github:

``` r
# install.packages("devtools")
devtools::install_github("DongqiangZeng0808/TMEscore")
```
Main documentation is on the `tmescore` function in the package:

``` r
help('tmescore', 'TMEscore')
```

References
----------
[1] Zeng DQ#, Li MY#, Zhou R, Zhang JW, Sun HY, Shi M, Bin JP, Liao YL, Rao JJ, Liao WJ*, Tumor microenvironment characterization in gastric cancer identifies prognostic and imunotherapeutically relevant gene signatures. 2019, Cancer Immunology Research. 
DOI: 10.1158/2326-6066.CIR-18-0436
PMID: *30842092_
.. _30842092: https://www.ncbi.nlm.nih.gov/pubmed/30842092
