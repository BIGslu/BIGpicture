#' kmFit example model results
#'
#' @format List of 2:
#' \enumerate{
#' \item \strong{lme} A data frame with 1000 rows and 7 columns
#' \describe{
#'   \item{model}{character. Model name.}
#'   \item{gene}{character. ENSEMBL gene ID.}
#'   \item{variable}{character. Variable name.}
#'   \item{estimate}{numeric. Model estimate, log fold change.}
#'   \item{pval}{numeric. Significance estimate.}
#'   \item{sigma}{numeric. Model fit, mean residual error.}
#'   \item{FDR}{numeric. FDR corrected signficant estimate.}
#'       }
#'
#' \item \strong{lmekin} A data frame with 2000 rows and 7 columns
#' \describe{
#'   \item{model}{character. Model name.}
#'   \item{gene}{character. ENSEMBL gene ID.}
#'   \item{variable}{character. Variable name.}
#'   \item{estimate}{numeric. Model estimate, log fold change.}
#'   \item{pval}{numeric. Significance estimate.}
#'   \item{sigma}{numeric. Model fit, mean residual error.}
#'   \item{FDR}{numeric. FDR corrected signficant estimate.}
#'       }
#' }
#' @source \url{https://github.com/altman-lab/P259_pDC_public}
#' @references Dill-McFarland et al. 2021. Eosinophil-mediated suppression and Anti-IL-5 enhancement of plasmacytoid dendritic cell interferon responses in asthma. J Allergy Clin Immunol. In revision
#' @description A dataset containing linear mixed effects model results with and without genetic kinship correction. Raw data were obtained from the kimma package.
#' @docType data
#' @name example_model
#' @keywords datasets
"example_model"
