#' Plot PCA colored by variables of interest
#'
#' @param dat edgeR DGEList or limma EList object containing gene counts in libraries
#' @param counts Data frame containing gene counts in libraries
#' @param meta Data frame containing meta data with vars. Only needed if counts is a counts table
#' @param vars Character vector of variables to color PCA by
#' @param PCx Numeric value for PC to plot on x-axis. Default it 1
#' @param PCy Numeric value for PC to plot on y-axis. Default it 2
#' @param scale Logical if should scale variance in PCA calculation see stats::prcomp for details. Default is FALSE
#' @param outlier_sd Numeric. If vars includes "outlier", statistical outliers are determined and colored based on this standard deviation along PC1 and PC2.
#' @param outlier_group Character string in which to group sd calculations
#' @param outlier_label Character string of variable to label outlying libraries with
#' @param transform_logCPM Logical if should convert counts to log counts per million
#' @param libraryID Character of variable name to match dat meta data frames
#'
#' @return List of ggplot objects
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' plot_pca(kimma::example.voom, vars=c("virus","outlier"))
#' plot_pca(kimma::example.voom, vars=c("virus","outlier"), PCx=1, PCy=3)

plot_pca <- function(dat = NULL,
                     counts = NULL, meta = NULL,
                     vars, PCx=1, PCy=2,
                     scale = FALSE, outlier_sd = 3,
                     outlier_group = NULL, outlier_label = NULL,
                     transform_logCPM = FALSE,
                     libraryID = "libID"){
  PCx.max <- PCx.mean <- PCx.min <- PCx.sd <- PCy.max <- PCy.mean <- PCy.min <- PCy.sd <- col.group <- NULL

  PCx_name <- paste0("PC", PCx)
  PCy_name <- paste0("PC", PCy)

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

  PC1.label <- paste("PC",PCx, " (", summary(PCA)$importance[2,PCx]*100, "%)", sep="")
  PC2.label <- paste("PC", PCy, " (", summary(PCA)$importance[2,PCy]*100, "%)", sep="")

  # Extract PC values
  pca.dat <- as.data.frame(PCA$x) %>%
    tibble::rownames_to_column(libraryID) %>%
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
                         PCx.mean = mean(.data[[PCx_name]]),
                         PCx.sd = stats::sd(.data[[PCx_name]]),
                         PCy.mean = mean(.data[[PCy_name]]),
                         PCy.sd = stats::sd(.data[[PCy_name]])) %>%
        #Calculate +/- sd limits
        dplyr::mutate(
          PCx.min = PCx.mean-(outlier_sd*PCx.sd),
          PCx.max = PCx.mean+(outlier_sd*PCx.sd),
          PCy.min = PCy.mean-(outlier_sd*PCy.sd),
          PCy.max = PCy.mean+(outlier_sd*PCy.sd)) %>%
        #add to PCA data
        dplyr::full_join(pca.dat, by=outlier_group) %>%
        #ID potential outliers
        dplyr::mutate(col.group = ifelse(
          .data[[PCx_name]] > PCx.max |
            .data[[PCx_name]] < PCx.min |
            .data[[PCy_name]] > PCy.max |
            .data[[PCy_name]] < PCy.min,
          "yes", "no"))
    } else {
      dat.sd <- pca.dat %>%
        dplyr::ungroup() %>%
        #Calculate PC mean std deviation
        dplyr::summarise(PCx.mean = mean(.data[[PCx_name]]),
                         PCx.sd = stats::sd(.data[[PCx_name]]),
                         PCy.mean = mean(.data[[PCy_name]]),
                         PCy.sd = stats::sd(.data[[PCy_name]])) %>%
        #Calculate +/- sd limits
        dplyr::mutate(
          PCx.min = PCx.mean-(outlier_sd*PCx.sd),
          PCx.max = PCx.mean+(outlier_sd*PCx.sd),
          PCy.min = PCy.mean-(outlier_sd*PCy.sd),
          PCy.max = PCy.mean+(outlier_sd*PCy.sd))

      #ID potential outliers
      pca.dat.sd <- pca.dat %>%
        dplyr::mutate(col.group = ifelse(
          .data[[PCx_name]] > dat.sd$PCx.max |
            .data[[PCx_name]] < dat.sd$PCx.min |
            .data[[PCy_name]] > dat.sd$PCy.max |
            .data[[PCy_name]] < dat.sd$PCy.min,
          "yes", "no"))
      }

    plot2 <- ggplot2::ggplot(pca.dat.sd,
                             ggplot2::aes(x = .data[[PCx_name]],
                                          y = .data[[PCy_name]],
                                          color = col.group)) +
      ggplot2::geom_point(size=3) +
      #Beautify
      ggplot2::theme_classic() +
      ggplot2::labs(x=PC1.label, y=PC2.label,
                    color=paste0("Std dev > ", outlier_sd, "X")) +
      ggplot2::coord_fixed(ratio=1) +
      ggplot2::scale_color_manual(values = c("#969696","#b10026"))

    if(!is.null(outlier_label)){
      plot2 <- plot2 +
        ggrepel::geom_text_repel(data=dplyr::filter(pca.dat.sd,
                                                    col.group == "yes"),
                                 ggplot2::aes(label = .data[[outlier_label]]),
                                 show.legend = FALSE, max.overlaps = Inf)
    }

    plot.ls[["outlier"]] <- plot2
  }

  return(plot.ls)
}
