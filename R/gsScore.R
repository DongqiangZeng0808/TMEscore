


##' Calculate score across genes and samples
##'
##' This wrapper function combines filtering out genes with low reads in a number of samples (recommended for limma:voom) with normalization
##' @param gm normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
##' @param summarizationFunction character vector defining whether signature scores shall be based on principal component 1 ("PC", default) or z-scores (other value)
##' @return numeric vector of length ncol(gm); a score summarizing the rows of gm
##' @author Dorothee Nickles
##' @export
gsScore <- function(gm, summarizationFunction="PC") {
  if (summarizationFunction == "PC") {
    pc <- prcomp(t(gm),
                 retx=TRUE)
    gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
  } else {
    gss <- colMeans(gm)
  }
  return(gss)
}
