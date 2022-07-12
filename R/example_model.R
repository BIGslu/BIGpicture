#' kmFit example model results
#'
#' @format List of 6:
#' Results for linear mixed effects (lme) or linear mixed effects model with kinship relatedness (lmerel) for `~ virus*asthma + (1|ptID)`
#' \enumerate{
#' \item \strong{lme/lmerel} A data frame with 1000 rows and 10 columns
#' \describe{
#'   \item{model}{character. Model name.}
#'   \item{gene}{character. ENSEMBL gene ID.}
#'   \item{variable}{character. Variable name.}
#'   \item{estimate}{numeric. Model estimate, i.e. log fold change.}
#'   \item{pval}{numeric. Significance estimate.}
#'   \item{FDR}{numeric. FDR corrected signficant estimate.}
#'   \item{hgnc_symbol}{character. Current HUGO/HGNC gene symbol.}
#'   \item{Previous symbols}{comma separated character. Previous HUGO/HGNC gene symbols, if any.}
#'   \item{Alias symbols}{comma separated character. Alias HUGO/HGNC gene symbols, if any.}
#'   \item{gene_biotype}{character. Gene biotype. All protein coding in this data set.}
#'       }
#'
#' \item \strong{lme.contrast/lmerel.contrast} A data frame with 8000 rows and 15 columns
#' \describe{
#'   \item{model}{character. Model name.}
#'   \item{gene}{character. ENSEMBL gene ID.}
#'   \item{variable}{character. Variable name.}
#'   \item{contrast_ref}{character. Model name.}
#'   \item{contrast_lvl}{character. Model name.}
#'   \item{estimate}{numeric. Model estimate, i.e. log fold change.}
#'   \item{pval}{numeric. Significance estimate.}
#'   \item{FDR}{numeric. FDR corrected signficant estimate.}
#'   \item{std.error}{numeric. Dtandard error.}
#'   \item{df}{numeric. Degrees of freedom.}
#'   \item{statistic}{numeric. T statistic.}
#'   \item{hgnc_symbol}{character. Current HUGO/HGNC gene symbol.}
#'   \item{Previous symbols}{comma separated character. Previous HUGO/HGNC gene symbols, if any.}
#'   \item{Alias symbols}{comma separated character. Alias HUGO/HGNC gene symbols, if any.}
#'   \item{gene_biotype}{character. Gene biotype. All protein coding in this data set.}
#'       }
#'
#' \item \strong{lme.fit/lmerel.fit} A data frame with 1000 rows and 11 columns
#' \describe{
#'   \item{model}{character. Model name.}
#'   \item{gene}{character. ENSEMBL gene ID.}
#'   \item{sigma}{numeric. Estimated standard deviation of the errors, ie residual standard deviation.}
#'   \item{AIC}{numeric. Akaike information criterion.}
#'   \item{BIC}{numeric. Bayesian information criterion.}
#'   \item{Rsq}{numeric. R-squared.}
#'   \item{adj_Rsq}{numeric. Adjusted R-squared.}
#'   \item{hgnc_symbol}{character. Current HUGO/HGNC gene symbol.}
#'   \item{Previous symbols}{comma separated character. Previous HUGO/HGNC gene symbols, if any.}
#'   \item{Alias symbols}{comma separated character. Alias HUGO/HGNC gene symbols, if any.}
#'   \item{gene_biotype}{character. Gene biotype. All protein coding in this data set.}
#'       }
#' }
#' @source \url{https://github.com/altman-lab/P259_pDC_public}
#' @references Dill-McFarland et al. 2021. Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. In revision
#' @description A dataset containing linear mixed effects model results with and without genetic kinship correction. Raw data were obtained from the kimma package.
#' @docType data
#' @name example_model
#' @keywords datasets
"example_model"
