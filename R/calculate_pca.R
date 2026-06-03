#' Calculate PCA. Input either voom object or counts and metadata tables
#' @param dat edgeR DGEList, or limma EList object containing gene counts in libraries. If provided, counts and meta are ignored
#' @param counts Data frame or matrix containing gene counts in libraries
#' @param meta Data frame or matrix containing meta data with vars. Only needed if using counts, not dat
#' @param scale Logical if should scale variance in PCA calculation see stats::prcomp for details. Default is FALSE
#' @param transform_logCPM Logical if should convert counts to log counts per million
#' @param force If you already ran PCA but want to overwrite previous results set force=TRUE (default is FALSE)
#' @param libraryID Character of variable name to match dat meta data frames
#'
#' @returns If a voom object is input, returns voom object with $PCA.scaled or $PCA.unscaled slot depending of if scale=TRUE or FALSE
#' If a data frame or matrix of counts is input, returns a list with slots $counts, $samples, and $PCA.scaled or $PCA.unscaled slot depending of if scale=TRUE or FALSE
#' @export
#'
#' @examples
#' calculate_pca(dat = kimma::example.voom, scale = TRUE)
#' calculate_pca(counts = kimma::example.count,
#'               meta = kimma::example.dat$samples,
#'               transform_logCPM = TRUE)

calculate_pca <- function(dat = NULL,
                          counts = NULL, meta = NULL,
                          scale = FALSE, transform_logCPM = FALSE,
                          force = FALSE,
                          libraryID = "libID"){

  # PC1 <- PC2 <- PC1.max <- PC1.mean <- PC1.min <- PC1.sd <- PC2.max <- PC2.mean <- PC2.min <- PC2.sd <- col.group <- libID <- sd  <- NULL
  # PCA <- PC1.label <- PC2.label <- PC1 <- sd <- PC2 <- PC1.mean <- PC1.sd <- PC2.mean <- PC2.sd <- PC1.max <- PC1.min <- PC2.max <- PC2.min <- col.group <- NULL
  PC <- NULL

  # check if PCA data already exists
  # if dat is a df or matrix then it wont have the PCA data
  if(!(is.null(dat))){
    if(scale){
      if(!is.null(dat$PCA.scaled) & !force){
        stop("You've already calculated scaled PCA. If you want to re-run, please set force=TRUE")
      }
    }else{
      if(!is.null(dat$PCA.unscaled) & !force){
        stop("You've already calculated unscaled PCA. If you want to re-run, please set force=TRUE")
      }
    }
  }

  #common errors
  if(!is.null(dat) & !is.null(counts)){
    stop("Only provide one of dat or counts.")
  }
  if(!is.null(counts) & is.null(meta)){
    stop("meta must be provided when counts is used.")
  }

  #Extract metadata table
  if(any(class(dat) == "DGEList")){
    meta.df <- dat$samples
  } else if (any(class(dat) == "EList")){
    meta.df <- dat$targets
  } else { meta.df <- meta }

  #Extract counts table
  if(any(class(dat) == "DGEList")){
    count.df <- dat$counts
  } else if (any(class(dat) == "EList")){
    count.df <- dat$E
  } else { count.df <- counts }

  #Move rownames if in data frame
  if(!is.numeric(as.matrix(count.df))){
    count.mat <- as.matrix(count.df[,-1])
    rownames(count.mat) <- unlist(count.df[,1], use.names = FALSE)
  } else { count.mat <- as.matrix(count.df) }

  #Convert to logCPM if selected
  if(transform_logCPM){
    count.mat.format <- edgeR::cpm(count.mat, log=TRUE)
  } else { count.mat.format <- count.mat }

  #Calculate
  set.seed(8456)
  PCA <- stats::prcomp(t(count.mat.format), scale.=scale, center=TRUE)

  # Extract PC values
  pca.dat <- as.data.frame(t(PCA$x)) %>%
    tibble::rownames_to_column("PC") %>%
    dplyr::mutate(pct.var = summary(PCA)$importance[2,PC])

  # If not a voom object, create list and add PCA data
  if(is.null(dat)){
    dat.out <- list()
    dat.out$counts <- counts
    dat.out$samples <- meta
    if(scale){
      dat.out$PCA.scaled <- pca.dat
    }else{
      dat.out$PCA.unscaled <- pca.dat
    }
  }else{ # if voom or edgeR, add pca data to $PCA.scaled/PCA.unscaled
    dat.out <- dat
    if(scale){
      dat.out$PCA.scaled <- pca.dat
    }else{
      dat.out$PCA.unscaled <- pca.dat
    }
  }

  return(dat.out)

}
