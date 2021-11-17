#' Plot PCA colored by variables of interest
#'
#' @param dat Data frame, edgeR DGEList, or limma EList object containing gene counts in libraries
#' @param meta Data frame containing meta data with vars. Only needed if dat is a counts table and not an edgeR or limma object
#' @param vars Character vector of variables to color PCA by
#' @param outlier_sd Numeric. If vars includes "outlier", statistical outliers are determined and colored based on this standard deviation along PC1 and PC2.
#' @param outlier_group Character string in which to group sd calculations
#' @param transform_logCPM Logical if should convert counts to log counts per million
#' @param libraryID Character of variable name to match dat meta data frames
#'
#' @return List of ggplot objects
#' @export

plot_pca <- function(dat, meta = NULL, vars, outlier_sd = 3,
                     outlier_group = NULL, transform_logCPM = FALSE,
                     libraryID = "libID"){

  PC1 <- PC2 <- PC1.max <- PC1.mean <- PC1.min <- PC1.sd <- PC2.max <- PC2.mean <- PC2.min <- PC2.sd <- col.group <- libID <- sd  <- NULL

  #common errors
  if((is.data.frame(dat) | is.matrix(dat)) & is.null(meta)){
    stop("meta must be provided when dat is a counts table.")
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
  PCA <- stats::prcomp(t(count.mat.format), scale.=TRUE, center=TRUE)

  PC1.label <- paste("PC1 (", summary(PCA)$importance[2,1]*100, "%)", sep="")
  PC2.label <- paste("PC2 (", summary(PCA)$importance[2,2]*100, "%)", sep="")

  # Extract PC values
  pca.dat <- as.data.frame(PCA$x) %>%
    tibble::rownames_to_column(libraryID) %>%
    # Merge with metadata
    dplyr::left_join(meta.df)

  # plots
  plot.ls <- list()
  for(var in vars[vars != "outlier"]){
    pca.plot <- ggplot2::ggplot(pca.dat, ggplot2::aes_string("PC1", "PC2", color=var)) +
      ggplot2::geom_point(size=3) +
      #Beautify
      ggplot2::theme_classic() +
      ggplot2::labs(x=PC1.label, y=PC2.label, color=var) +
      ggplot2::coord_fixed(ratio=1)

    #Add connecting lines if duplicate
    if(var %in% c("duplicate","dupID")){
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
