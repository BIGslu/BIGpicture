#' Plotting function for hypergeometric enrichment
#'
#' @param df Data frame output by SEARchways::BIGprofiler or SEARchways::flexEnrich
#' @param fdr_cutoff Numeric. Maximum FDR to plot. Default is 0.2.
#' @param pathway_col Character string Column name for df column containing gene set names.
#' @param gene_col Character string. Column name for df column containing gene lists names. Required if "y_grouping_method == 'hclust'".
#' @param ratio_col Character string Column name for df column containing k/K ratios.
#' @param fdr_col Character string. Column name for df column containing corrected p-values to plot.
#' @param gssize_col Character string. Column name for df column containing gene set sizes.
#' @param y_grouping_method Character string. Method for grouping gene sets along the y-axis. "hclust" (based on gene membership), "overlap_size", "gs_size", "ratio", "fdr", and "input". Default is "fdr".
#' @param include_gssize Boolean. Whether or not to include a column of gene set sizes to the left of the enrichment plot. Default is FALSE.
#' @param chart_style Character string. Options are "lollipop", "bar", and "dot". Default is "lollipop".
#' @param fdr_binned_colors Numeric vector. Cutoffs for binned FDR value color groups for lollipop plots. Default is c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.99).
#' @param fdr_continuous Boolean. Whether or not to color points in a continuous fashions by a negative log transformed FDR value. Only applies to lollipop plots and dotplots.
#' @param dot_sig_cutoff Numeric. What cutoff to use for applying a black outline to dotplot dots
#' @param dot_groupcol Character string. Name of the column used to group enrichment along the x axis for dotplot.
#' @return ggplot2 object
#' @export
#'
#' @examples
#' library(SEARchways)
#' library(dplyr)
#' #Run enrichment
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]),
#'                   HRV2 = names(example.gene.list[[2]]))
#' enrich <- flexEnrich(gene_list, ID="ENSEMBL", category="H")
#'
#' #Plot
#' plot_enrich2(enrich, fdr_cutoff = 0.5)

plot_enrich2 <- function(df = NULL,
                         fdr_cutoff = 0.2,
                         pathway_col = "pathway",
                         gene_col = "genes",
                         ratio_col = "k/K",
                         fdr_col = "FDR",
                         gssize_col = "n_pathway_genes",
                         y_grouping_method = "fdr",
                         include_gssize = FALSE,
                         chart_style = "lollipop",
                         fdr_binned_colors = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.99),
                         fdr_continuous = FALSE,
                         dot_sig_cutoff = NULL,
                         dot_groupcol = "group"
){

  fdr <- size <- gs <- gssize <- ratio <- desc <- prev <- gene <- geneset <- color <- Significance <- gssize_bin <- `k/K` <- sigyn <- group <- NULL


  ### Format input ###
  # rename columns
  df$gs <- df[[pathway_col]]
  df$gns <- df[[gene_col]]
  df$fdr <- df[[fdr_col]]
  df$ratio <- df[[ratio_col]]


  if(!is.null(gssize_col)){
    df$gssize <- df[[gssize_col]]
  }

  # subset by FDR
  if(is.null(df[[fdr_col]])){
    stop("Your column name indicating FDR values is not valid. Please check.")
  } else {
    gssub <- df %>%
      dplyr::filter(fdr < fdr_cutoff) %>%
      dplyr::pull(gs) %>% unique()
    df <- df %>%
      dplyr::filter(gs %in% gssub)
  }

  # get gene vector
  genevec <- c()
  for(i in 1:nrow(df)){
    rowgenes <- df$gns[i]
    rowgenes <- unlist(rowgenes)
    genevec <- c(genevec,rowgenes)
    genevec <- unique(genevec)
  }

  # get pathways
  gsvec <- unique(df$gs)

  # get count mat for overlap, clustering
  count_df <- matrix(ncol = length(genevec), nrow = length(gsvec))
  colnames(count_df) <- genevec
  rownames(count_df) <- gsvec
  ## make binary matrix
  for(i in 1:nrow(count_df)){
    set <- rownames(count_df)[i]
    rowgenes <- df$gns[which(df$gs == set)]
    rowgenes <- unlist(rowgenes)
    for(j in 1:ncol(count_df)){
      g <- colnames(count_df)[j]
      val <- ifelse(g %in% rowgenes, 1,0)
      count_df[set,g] <- val
    }
  }

  ### Get order vector for gene sets ###

  # get reference level for multi-group dotplot
  if(!is.null(dot_groupcol) & chart_style == "dot"){
    if(is.factor(df[[dot_groupcol]])){
      lev <- levels(df[[dot_groupcol]])[1]
    } else{
      lev <- levels(as.factor(df[[dot_groupcol]]))[1]
    }

    df_lev <- df
    df_lev$group <- df[[dot_groupcol]]
    df_lev <- df_lev %>%
      dplyr::filter(`group` == as.character(lev))
  } else{
    df_lev <- df
  }


  if(y_grouping_method == "hclust"){ # hierarchical clustering by shared genes
    y_cluster_df <- textshape::cluster_matrix(count_df, dim = "row")
    y_levels <- rownames(y_cluster_df)
  } else if(y_grouping_method == "overlap_size"){ # number of shared genes between set and plotted genes
    y_levels <- data.frame("gs" = rownames(count_df), "size" = rowSums(count_df)) %>%
      dplyr::arrange(size) %>%
      dplyr::pull(gs) %>% unique()
  } else if(y_grouping_method == "gs_size"){ # size of gene set
    y_levels <- df_lev %>%
      dplyr::arrange(gssize) %>%
      dplyr::pull(gs) %>% unique()
  } else if(y_grouping_method == "ratio"){ # ratio of set genes in query
    y_levels <- df_lev %>%
      dplyr::arrange(ratio) %>%
      dplyr::pull(gs) %>% unique()
  } else if(y_grouping_method == "fdr"){
    y_levels <- df_lev %>%
      dplyr::arrange(desc(fdr)) %>%
      dplyr::pull(gs) %>% unique()
  } else if(y_grouping_method == "input"){
    y_levels <- df_lev %>%
      dplyr::pull(gs) %>% unique()
  } else{stop('Valid options for y_grouping_method are "hclust", "overlap_size", "gs_size", "ratio", "fdr", and "input".')}

  titlesize <- 12
  haxis_size <- 10
  vaxis_size <- 12

  # GS Size column
  p3 <- ggplot2::ggplot(df_lev, ggplot2::aes(x = factor(gs , levels = y_levels), y = rep(1, length(y_levels)))) +
    ggplot2::ggtitle("K") +
    ggplot2::geom_text(ggplot2::aes(label = gssize), size = 4) +
    ggplot2::coord_flip() +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = vaxis_size),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(1,1,1,1, "pt"))+
    ggplot2::scale_x_discrete(position = "bottom")

  ### Make charts of k/K and fdr values ###
  if(chart_style == "bar"){
    # k/K bar
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x=factor(gs , levels = y_levels), y=ratio)) +
      ggplot2::geom_bar(stat = "identity", fill = "#f43545") +
      ggplot2::labs(y="k/K") +
      ggplot2::theme(panel.grid = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = 1),
                     plot.margin = ggplot2::margin(1,1,1,1, "pt"),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.line.x = element_line(color = "black", linewidth = 0.5)
      )+
      ggplot2::scale_y_reverse(position = "left", labels = scales::label_number(accuracy = 0.01)) +
      ggplot2::coord_flip()

    if(include_gssize){
      p2 <- p2 +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }
    # FDR bar
    p4 <- ggplot2::ggplot(df, ggplot2::aes(x= factor(gs , levels = y_levels), y= -log10(fdr))) +
      ggplot2::geom_bar(stat = "identity", fill = "#fa9801") +
      ggplot2::labs(y="-log10 FDR") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = 1),
                     panel.grid = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y.left = ggplot2::element_blank(),
                     axis.text.y.right = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(1,1,1,1, "pt"),
                     axis.line.x = element_line(color = "black", linewidth = 0.5)
      )+
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(position = "left",labels = scales::label_number(accuracy = 1))
  } else if(chart_style == "lollipop"){
    df_lp <- df
    ## get FDR colors ##
    if(fdr_continuous == FALSE){ # binned categorical
      fdr_colors.sort <- sort(fdr_binned_colors)
      df_lp$Significance <- NA
      fdr_levels <- c()
      for(i in 1:length(fdr_colors.sort)){
        if(i==1){
          df_lp$Significance[df_lp$FDR < fdr_colors.sort[i]] <- paste("FDR < ", fdr_colors.sort[i])
          fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
        } else{
          df_lp$Significance[df_lp$FDR < fdr_colors.sort[i] & df_lp$FDR >= fdr_colors.sort[i-1]] <- paste("FDR < ", fdr_colors.sort[i])
          fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
        }
      }
    } else{ # continuous log transform
      df_lp$Significance <- -log10(df_lp$FDR)
    }

    size_bins <- c("< 10" = 1, "10 - 49" = 2, "50 - 99" = 3,
                   "100 - 499" = 4, "500 - 999" = 5, "1000 - 2500" = 6,
                   "2500 - 5000" = 7, ">= 5000" = 8)
    df_lp$gssize_bin <- NA

    df_lp$gssize_bin = ifelse(df_lp$gssize < 10, "< 10",
                              ifelse(df_lp$gssize < 50, "10 - 49",
                                     ifelse(df_lp$gssize < 100, "50 - 99",
                                            ifelse(df_lp$gssize < 500, "100 - 499",
                                                   ifelse(df_lp$gssize < 1000, "500 - 999",
                                                          ifelse(df_lp$gssize < 2500, "1000 - 2500",
                                                                 ifelse(df_lp$gssize < 5000, "2500 - 5000",">= 5000")))))))
    df_lp <- df_lp %>%
      dplyr::mutate(gssize_bin = factor(gssize_bin, levels = names(size_bins)))

    p5 <- ggplot2::ggplot(df_lp, ggplot2::aes(x = factor(gs, levels = y_levels), y = ratio)) +
      ggplot2::geom_segment(ggplot2::aes(factor(gs, levels = y_levels),
                                         xend=gs, y=0, yend=ratio)) +
      ggplot2::geom_point(ggplot2::aes(fill = Significance,
                                       size = gssize_bin),
                          shape=21, stroke=1) +
      ggplot2::scale_size_manual(values = size_bins) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = 1),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(1,1,1,1, "pt"),
                     axis.line.x = element_line(color = "black", linewidth = 0.5)) +

      ggplot2::scale_y_continuous(position = "left",labels = scales::label_number(accuracy = 0.001)) +
      ggplot2::labs(y="k/K ratio")

    if(include_gssize){
      p5 <- p5 +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

    if(fdr_continuous == FALSE){
      p5 <- p5 +
        ggplot2::scale_fill_brewer(palette = "RdYlBu", na.value="grey70") +
        ggplot2::guides(size=ggplot2::guide_legend(title="k/K ratio"),
                        fill=ggplot2::guide_legend(title="Significance", override.aes = list(size=5)))
    } else{
      p5 <- p5 +
        ggplot2::scale_fill_gradient(low = "#FDE725", high = "red3", na.value = "red3") +
        ggplot2::guides(size=ggplot2::guide_legend(title="k/K ratio"),
                        fill=ggplot2::guide_legend(title="-log10(FDR)", override.aes = list(size=5)))
    }

  } else if(chart_style == "dot"){

    ## get grouping column if there is one ##
    df_dp <- df
    if(!is.null(dot_groupcol)){
      df_dp$group <- df_dp[[dot_groupcol]]
    } else{
      df_dp$group <- "enrichment"
    }

    ## get FDR colors ##
    if(fdr_continuous == FALSE){ # binned categorical
      fdr_colors.sort <- sort(fdr_binned_colors)
      df_dp$Significance <- NA
      fdr_levels <- c()
      for(i in 1:length(fdr_colors.sort)){
        if(i==1){
          df_dp$Significance[df_dp$FDR < fdr_colors.sort[i]] <- paste("FDR < ", fdr_colors.sort[i])
          fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
        } else{
          df_dp$Significance[df_dp$FDR < fdr_colors.sort[i] & df_dp$FDR >= fdr_colors.sort[i-1]] <- paste("FDR < ", fdr_colors.sort[i])
          fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
        }
      }
    } else{ # continuous log transform
      df_dp$Significance <- -log10(df_dp$FDR)
    }
    if(!is.null(dot_sig_cutoff)){ # add a column to indicate whether or not a black circle will be added
      df_dp <- df_dp %>%
        dplyr::mutate("sigyn" = ifelse(df_dp$FDR < dot_sig_cutoff, "yes", "no"))
    } else{
      df_dp <- df_dp %>%
        dplyr::mutate("sigyn" = "no")
    }

    colvec <- c("yes" = "black", "no" = "#FFFFFF00") # these are the colors for the outline strokes

    p6 <- ggplot2::ggplot(df_dp, ggplot2::aes(x = factor(group), y = factor(gs, levels = y_levels))) +
      ggplot2::geom_point(ggplot2::aes(fill = Significance,
                                       size = `k/K`,
                                       color = sigyn), shape=21, stroke=1) +
      ggplot2::scale_color_manual(values = colvec, breaks = c("yes")) +
      ggplot2::scale_x_discrete(expand = c(0.3,0.3,0.3,0.3)) +
      ggplot2::scale_size_area(max_size = 10, breaks = scales::breaks_extended(n = 6)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = vaxis_size),
                     plot.margin = ggplot2::margin(1,1,1,1, "pt"),
                     axis.title = ggplot2::element_blank(),
                     axis.line.x = element_line(color = "black", linewidth = 0.5))

    if(fdr_continuous == FALSE){
      p6 <- p6 +
        ggplot2::scale_fill_brewer(palette = "RdYlBu", na.value="grey70") +
        ggplot2::guides(size=ggplot2::guide_legend(title="k/K ratio"),
                        fill=ggplot2::guide_legend(title="Significance", override.aes = list(size=5)))
    } else{
      p6 <- p6 +
        ggplot2::scale_fill_gradient(low = "#FDE725", high = "red3", na.value = "red3") +
        ggplot2::guides(size=ggplot2::guide_legend(title="k/K ratio"),
                        fill=ggplot2::guide_legend(title="-log10(FDR)", override.aes = list(size=5)))
    }

    if(!is.null(dot_sig_cutoff)){
      p6 <- p6 +
        ggplot2::guides(color=ggplot2::guide_legend(title=paste0("FDR < ", dot_sig_cutoff), override.aes = list(size=5)))
    }


    if(include_gssize){
      p6 <- p6 +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

  }else(stop('Valid options chart_style are "bar", "lollipop", and "dotplot".'))

  ### Final patchwork ###
  if(chart_style == "bar"){
    if(include_gssize){
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p2 + # k/K bar
        patchwork::plot_spacer() +
        p4 + # FDR bar
        patchwork::plot_layout(widths = c(1, -0.2,  3 , -0.2, 3))
    } else{
      p <- p2 + # k/K bar
        patchwork::plot_spacer() +
        p4 + # FDR bar
        patchwork::plot_layout(widths = c(3 , -0.3, 3))
    }
  } else if(chart_style == "lollipop"){
    if(include_gssize){
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p5 + # lollipop
        patchwork::plot_layout(widths = c(1.2, -0.2,  6))
    } else{
      p <- p5
    }
  } else if(chart_style == "dot"){
    if(include_gssize){
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p6 + # dotplot
        patchwork::plot_layout(widths = c(1.2, -0.2,  6))
    } else{
      p <- p6
    }
  } else(stop('Valid options chart_style are "bar", "lollipop", and "dot".'))

  return(p)

}
