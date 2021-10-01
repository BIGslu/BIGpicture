#' Plot PCA colored by variables of interest
#'
#' @param dat Data frame, edgeR DGEList, or limma EList object containing gene counts in libraries
#' @param meta Data frame containing meta data with vars. Only needed if dat is a counts table and not an edgeR or limma object
#' @param vars Character vector of variables to color PCA by
#' @param transform_logCPM Logical if should convert counts to log counts per million
#'
#' @return List of ggplot objects
#' @export

plot_pca <- function(dat, meta = NULL, vars, transform_logCPM = FALSE){

  PC1 <- PC2 <- NULL

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
  } else { count.df <- dat }

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
  PCA <- stats::prcomp(t(count.mat.format))

  PC1.label <- paste("PC1 (", summary(PCA)$importance[2,1]*100, "%)", sep="")
  PC2.label <- paste("PC2 (", summary(PCA)$importance[2,2]*100, "%)", sep="")

  # Extract PC values
  pca.dat <- as.data.frame(PCA$x) %>%
    tibble::rownames_to_column("libID") %>%
    # Merge with metadata
    dplyr::left_join(meta.df, by="libID")

  # plots
  plot.ls <- list()
  for(var in vars){
    pca.plot <- ggplot2::ggplot(pca.dat, ggplot2::aes(PC1, PC2, color=get(var))) +
      ggplot2::geom_point(size=3) +
      #Beautify
      ggplot2::theme_classic() +
      ggplot2::labs(x=PC1.label, y=PC2.label, color=var) +
      ggplot2::coord_fixed(ratio=1)

    plot.ls[[var]] <- pca.plot
  }

  return(plot.ls)
}
