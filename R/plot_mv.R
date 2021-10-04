#' Plot mean-variance trends of RNAseq gene expression
#'
#' @param dat edgeR DGEList or a limma EList object
#' @param design Character string of model to use in voom normalization. Only needed if dat is an unnormalized DGEList.
#'
#' @return ggplot2 object
#' @export

plot_mv <- function(dat, design = NULL){
  x <- y <- linex <- liney <- NULL

  #non-voom input
  if(class(dat) == "DGEList"){
    #check design
    if(is.null(design)){"Please provide a design for voom normalization."}
    #Extract metadata
    meta <- dat$samples
    #Make model matrix
    mm <- stats::model.matrix(stats::as.formula(design), dat$samples)
    #voom normalize
    dat.voom <- limma::voom(dat, mm, plot = FALSE, save.plot = TRUE)
  } else if(class(dat) == "EList"){
    #voom object
    dat.voom <- dat
    #Check if plot data exist
    if(is.null(dat$voom.xy)){
      stop("dat appears to be a limma EList object but does not contain plot data. Please re-run voom( ) with save.plot = TRUE or input an edgeR DGEList object instead.")
    }
  } else(
    stop("dat must be an edgeR DGEList or a limma EList object.")
  )

  #plot
  MV.plot <- data.frame(
    x = dat.voom$voom.xy$x,
    y = dat.voom$voom.xy$y,
    linex = dat.voom$voom.line$x,
    liney = dat.voom$voom.line$y) %>%

    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=x, y=y), size=0.5) +
    ggplot2::geom_path(ggplot2::aes(x=linex, y=liney), color="red") +
    ggplot2::theme_classic() +
    ggplot2::labs(x="log2( count size + 0.5 )", y="Sqrt (stdev)",
         title="voom mean-variance trend")

  return(MV.plot)
}
