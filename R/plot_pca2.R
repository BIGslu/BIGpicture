#' Plot PCA colored by variables (does not calculate PCA)
#'
#' @param dat edgeR DGEList, limma EList, or list object containing gene counts in libraries, metadata, and PCA calculations from calculate_pca()
#' @param vars Character vector of variables to color PCA by
#' @param PCx Numeric value for PC to plot on x-axis. Default is 1
#' @param PCy Numeric value for PC to plot on y-axis. Default is 2
#' @param scale Logical if you want to use scaled PCA values or not. see stats::prcomp for details. Default is FALSE
#' @param outlier_sd Numeric. If vars includes "outlier", statistical outliers are determined and colored based on this standard deviation along PC1 and PC2.
#' @param outlier_group Character string in which to group sd calculations
#' @param libraryID Character of variable name to match dat meta data frames
#'
#' @returns List of ggplot objects
#' @export
#'
#' @examples
#' plot_pca(calculate_pca(kimma::example.voom), var=c("virus","outlier"))
#' plot_pca(calculate_pca(kimma::example.voom), var=c("virus","outlier"), PCx=1, PCy=3)

plot_pca2 <- function(dat, vars, PCx=1, PCy=2,
                      scale = FALSE, outlier_sd = 3,
                      outlier_group = NULL,
                      libraryID = "libID"){

  PC <- PC1 <- PC2 <- PC1.max <- PC1.mean <- PC1.min <- PC1.sd <- PC2 <- PC2.max <- PC2.mean <- PC2.min <- PC2.sd <- col.group <- pct.var <- sd  <- NULL

  # check for data format
  if(!(class(dat) %in% c("DGEList","EList","list"))){
    stop("Input should be DEGList, EList, or list. Maybe you need to run calculate_pca first?")
  }

  #Extract metadata table
  if (any(class(dat) == "EList")){
    meta.df <- dat$targets
  } else {
    meta.df <- dat$samples
  }

  #Extract PCA data
  if(scale){
    if(is.null(dat$PCA.scaled)){
      stop("Scaled PCA data not found. Try running calculate_pca on your data first with scale=TRUE")
    }else{
      pca.df <- dat$PCA.scaled
    }
  }else{
    if(is.null(dat$PCA.unscaled)){
      stop("Unscaled PCA data not found. Try running calculate_pca on your data first with scale=FALSE")
    }else{
      pca.df <- dat$PCA.unscaled
    }
  }

  importance.df <- pca.df %>%
    dplyr::select(PC,pct.var) %>%
    tidyr::pivot_wider(names_from = "PC",values_from = "pct.var")

  PC1.label <- paste("PC",PCx, " (", importance.df[PCx]*100, "%)", sep="")
  PC2.label <- paste("PC", PCy, " (", importance.df[PCy]*100, "%)", sep="")

  # transpose pca dat for plotting and add metadata
  pca.dat <- pca.df %>%
    dplyr::select(-pct.var) %>%
    tidyr::pivot_longer(-PC,names_to = libraryID,values_to = "val") %>%
    tidyr::pivot_wider(names_from = PC,values_from = "val") %>%
    # Merge with metadata
    dplyr::left_join(meta.df)

  # plots
  plot.ls <- list()
  for(var in vars[vars != "outlier"]){
    pca.plot <- ggplot2::ggplot(pca.dat, ggplot2::aes_string(paste0("PC", PCx),
                                                             paste0("PC", PCy),
                                                             color=var)) +
      ggplot2::geom_point(size=3) +
      #Beautify
      ggplot2::theme_classic() +
      ggplot2::labs(x=PC1.label, y=PC2.label, color=var) +
      ggplot2::coord_fixed(ratio=1)

    #Add connecting lines if duplicate
    if(var %in% c("duplicate","dupID","replicate","repID")){
      pca.plot <- pca.plot + ggplot2::geom_line(data=dplyr::filter(pca.dat, !is.na(get(var))))
    }
    plot.ls[[var]] <- pca.plot
  }

  #outlier plots
  if("outlier" %in% vars){
    #Group calculations if selected
    if(!is.null(outlier_group)){
      pca.dat.sd <- pca.dat %>%
        dplyr::group_by(dplyr::across(outlier_group)) %>%
        #Calculate PC mean std deviation
        dplyr::summarise(.groups="keep",
                         PC1.mean = mean(PC1),
                         PC1.sd = sd(PC1),
                         PC2.mean = mean(PC2),
                         PC2.sd = sd(PC2)) %>%
        #Calculate +/- sd limits
        dplyr::mutate(
          PC1.min = PC1.mean-(outlier_sd*PC1.sd),
          PC1.max = PC1.mean+(outlier_sd*PC1.sd),
          PC2.min = PC2.mean-(outlier_sd*PC2.sd),
          PC2.max = PC2.mean+(outlier_sd*PC2.sd)) %>%
        #add to PCA data
        dplyr::full_join(pca.dat, by=outlier_group) %>%
        #ID potential outliers
        dplyr::mutate(col.group = ifelse(PC1 > PC1.max | PC1 < PC1.min |
                                           PC2 > PC2.max | PC2 < PC2.min,
                                         "yes", "no"))
    } else {
      dat.sd <- pca.dat %>%
        #Calculate PC mean std deviation
        dplyr::summarise(.groups="keep",
                         PC1.mean = mean(PC1),
                         PC1.sd = sd(PC1),
                         PC2.mean = mean(PC2),
                         PC2.sd = sd(PC2)) %>%
        #Calculate +/- sd limits
        dplyr::mutate(
          PC1.min = PC1.mean-(outlier_sd*PC1.sd),
          PC1.max = PC1.mean+(outlier_sd*PC1.sd),
          PC2.min = PC2.mean-(outlier_sd*PC2.sd),
          PC2.max = PC2.mean+(outlier_sd*PC2.sd))

      #Calculate SD
      pca.dat.sd <- pca.dat %>%
        #ID potential outliers
        dplyr::mutate(col.group = ifelse(PC1 > dat.sd$PC1.max | PC1 < dat.sd$PC1.min |
                                           PC2 > dat.sd$PC2.max | PC2 < dat.sd$PC2.min,
                                         "yes", "no"))
    }

    plot2 <- ggplot2::ggplot(pca.dat.sd, ggplot2::aes(PC1, PC2, color=col.group)) +
      ggplot2::geom_point(size=3) +
      ggrepel::geom_text_repel(data=dplyr::filter(pca.dat.sd,
                                                  col.group == "yes"),
                               ggplot2::aes(label=get(libraryID)),
                               show.legend = FALSE, max.overlaps = Inf) +
      #Beautify
      ggplot2::theme_classic() +
      ggplot2::labs(x=PC1.label, y=PC2.label,
                    color=paste0("Std dev > ", outlier_sd, "X")) +
      ggplot2::coord_fixed(ratio=1) +
      ggplot2::scale_color_manual(values = c("#969696","#b10026"))

    plot.ls[["outlier"]] <- plot2
  }

  return(plot.ls)
}
