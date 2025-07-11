#' Plot PCA Correlation Heatmap (does not calculate PCA)
#'
#' @param dat edgeR DGEList, limma EList, or list object containing gene counts in libraries, metadata, and PCA calculations from calculate_pca()
#' @param vars Character vector of variables to color PCA by
#' @param scale Logical if should scale variance in PCA calculation see stats::prcomp for details. Default is FALSE
#' @param corr_type Character of correlation metric to use. Options: pearson or spearman. Default is pearson
#' @param pct_cutoff Numeric cutoff value for percent explained by PC (0 to 1). Default is 0.01 (1%)
#' @param flip_axes Logical whether to swap x and y axes. Default is FALSE
#' @param show_signif_p Logical whether to show correlation significance. Default is TRUE
#' @param signif_p_cutoff Numeric correlation p-value cutoff (0 to 1). Default is 0.05
#' @param libraryID Character of variable name to match dat meta data frames
#'
#' @returns ggplot2 object
#' @export
#'
#' @examples
#' dat_PCA <- calculate_pca(kimma::example.voom,scale=TRUE)
#' plot_pca_hm(dat_PCA,vars = c("median_cv_coverage","lib.size"),scale=TRUE,flip_axes = TRUE)
plot_pca_hm <- function(dat,vars,scale=FALSE,corr_type="pearson",
                        pct_cutoff = 0.01,flip_axes=FALSE,
                        show_signif_p=TRUE,signif_p_cutoff=0.05,
                        libraryID="libID"){

  P <- PC <- PC_num <- PClabel <- pct.var <- sig_label <- var <- NULL

  # check for dat format
  if(!(class(dat) %in% c("DGEList","EList","list"))){
    stop("Input should be DEGList, EList, or list. Maybe you need to run calculate_pca first?")
  }
  # check corr_type
  if(!(corr_type %in% c("pearson","spearman"))){
    stop(paste0("Unrecognized corr_type: ",corr_type))
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

  # remove PCs that explain less variation than pct_cutoff
  pca.df <- pca.df %>%
    dplyr::filter(pct.var >= pct_cutoff)

  # get percent var explained by PC
  importance.df <- pca.df %>%
    dplyr::select(PC,pct.var) %>%
    tidyr::pivot_wider(names_from = "PC",values_from = "pct.var")

  #create a PCA dataframe with metadata
  pca.dat <- pca.df %>%
    dplyr::select(-pct.var) %>%
    tidyr::pivot_longer(-PC,names_to = libraryID,values_to = "val") %>%
    tidyr::pivot_wider(names_from = PC,values_from = "val") %>%
    # Merge with metadata
    dplyr::left_join(meta.df)

  #get the names of all PCs in the data
  pcs <- pca.df$PC

  sub_pca.dat <- pca.dat %>%
    dplyr::select(dplyr::contains("PC"),dplyr::all_of(vars))

  # calculate correlations
  corr <- Hmisc::rcorr(as.matrix(sub_pca.dat),type=corr_type)

  # get the correlation results
  rvals <- as.data.frame(corr$r[pcs,vars]) %>%
    tibble::rownames_to_column(var="PC") %>%
    tidyr::pivot_longer(-PC,names_to = "var",values_to = "corr")
  sig <-  as.data.frame(corr$P[pcs,vars]) %>%
    tibble::rownames_to_column(var="PC") %>%
    tidyr::pivot_longer(-PC,names_to = "var",values_to = "P") %>%
    dplyr::mutate(sig_label = ifelse(P < signif_p_cutoff,"*",""))

  # create dataframe with data to plot
  plot_dat <- dplyr::full_join(rvals,sig,by=c("PC","var")) %>%
    dplyr::mutate(PClabel = paste(PC, " (", sprintf("%.1f", round(importance.df[PC]*100,digits=1)), "%)", sep="")) %>%
    dplyr::mutate(PC = factor(PC,levels = paste0("PC",seq(1,length(pcs))))) %>%
    dplyr::mutate(PC_num = as.numeric(gsub("PC","",PC))) %>%
    dplyr::arrange(PC_num)

  #get vector of ordered PC labels for plot
  ordered_pcs <- plot_dat %>%
    dplyr::pull(PClabel) %>% unique()

  # plot!
  if(flip_axes){
    plot_dat <- plot_dat %>%
      #reverse order of PCs so they show up correctly
      dplyr::mutate(PClabel = factor(PClabel,levels=rev(ordered_pcs)))

    plot <- ggplot2::ggplot(plot_dat,ggplot2::aes(x=var,y=PClabel,fill=corr)) +
      ggplot2::geom_tile() +
      ggplot2::theme_classic(base_size=16) +
      ggplot2::scale_x_discrete(expand=c(0,0),position="top") +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=-0.01)) +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,name = "Correlation") +
      ggplot2::coord_fixed()

    if(show_signif_p){
      plot <- plot +
        ggplot2::geom_text(ggplot2::aes(x=var,y=PClabel,label=sig_label),size=8)
    }
  }else{
    #order PCs by ordered variable
    plot_dat <- plot_dat %>%
      dplyr::mutate(PClabel = factor(PClabel,levels=ordered_pcs))

    plot <- ggplot2::ggplot(plot_dat,ggplot2::aes(x=PClabel,y=var,fill=corr)) +
      ggplot2::geom_tile() +
      ggplot2::theme_classic(base_size=16) +
      ggplot2::scale_x_discrete(expand=c(0,0),position="top") +
      ggplot2::scale_y_discrete(expand=c(0,0)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45,hjust=-0.01)) +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,name = "Correlation") +
      ggplot2::coord_fixed()
    if(show_signif_p){
      plot <- plot +
        ggplot2::geom_text(ggplot2::aes(x=PClabel,y=var,label=sig_label),size=8)
    }
  }
  return(plot)
}
