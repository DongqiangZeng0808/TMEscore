#TMEscore
=====================================
=====================================

TME infiltration patterns were determined and systematically correlated with TME cell phenotypes, genomic traits, and patient clinicopathological features to establish the "TMEscore" [1]

TMEscore is a R package to do Tumor microenvironment analysis. Main advantages:
- Provides functionality to calculate Tumor microenvironment (TME) score (PCA, ssGSEA, and z-score)
- Functions to visulize TME data.
- Deconvolution of cellular composition of tumor.
- Identify TME relevant mutaions.

=======
TME infiltration patterns were determined and systematically correlated with TME cell phenotypes, genomic traits, and patient clinicopathological features to establish the "TMEscore" [1]

TMEscore is a R package to do Tumor microenvironment analysis. Main advantages:
- Provides functionality to calculate Tumor microenvironment (TME) score (PCA, ssGSEA, and z-score)
- Functions to visulize TME data.
- Deconvolution of cellular composition of tumor.
- Identify TME relevant mutaions.

The package is not yet on CRAN. You can install from Github:

``` r
# install.packages("devtools")
devtools::install_github("DongqiangZeng0808/TMEscore")
```
Main documentation is on the `tmescore` function in the package:

``` r
help('tmescore')
```

Example

``` r
tmescore_acrg<-tmescore(eset = eset_acrg, pdata = pdata, method = "pca")
head(tmescore_acrg)
      ID      time    status km.cluster.IP TMEscoreA   TMEscoreB    TMEscore
1 GSM1523727  88.73      0   TMEclusterC -0.60373345 -0.9236988597  0.31996541
2 GSM1523728  88.23      0   TMEclusterC  0.71439724 -0.2958801840  1.01027742
3 GSM1523729  88.23      0   TMEclusterA -0.49132521 -0.0116628526 -0.47966236
4 GSM1523744 105.70      0   TMEclusterC  0.09007943 -0.7960630547  0.88614249
5 GSM1523745 105.53      0   TMEclusterB -0.06481076  0.0005292262 -0.06533999
6 GSM1523746  25.50      1   TMEclusterB -0.49198077  0.2359644650 -0.72794524
```

References
----------
[1] Zeng DQ#, Li MY#, Zhou R, Zhang JW, Sun HY, Shi M, Bin JP, Liao YL, Rao JJ, Liao WJ*, Tumor microenvironment characterization in gastric cancer identifies prognostic and imunotherapeutically relevant gene signatures. 2019, Cancer Immunology Research. 
DOI: 10.1158/2326-6066.CIR-18-0436
PMID: [30842092](https://www.ncbi.nlm.nih.gov/pubmed/30842092)
