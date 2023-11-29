#' Plotting function for hypergeometric enrichment
#'
#' @param df Data frame output by SEARchways::BIGprofiler or SEARchways::flexEnrich
#' @param fdr_cutoff Numeric. Maximum FDR to plot. Default is 0.2.
#' @param pathway_col String. Column name for df column containing gene set names.
#' @param gene_col String. Column name for df column containing gene lists names. Required if "include_grid == TRUE".
#' @param ratio_col String. Column name for df column containing k/K ratios.
#' @param fdr_col String. Column name for df column containing corrected p-values to plot.
#' @param gssize_col String. Column name for df column containing gene set sizes.
#' @param custom_genelist Optional vector of gene IDs to plot in the grid. Otherwise, all genes in the gene_col of input where set FDR < fdr_cutoff will be plotted.
#' @param x_grouping_method String. Method for grouping genes along the x-axis. Options are "hclust", "prevalence", "input", "alphabetical", and "geneset". Default is "geneset".
#' @param y_grouping_method String. Method for grouping gene sets along the y-axis. "hclust", "overlap_size", "gs_size", "ratio", "fdr", and "input". Default is "fdr".
#' @param prevalence_color String. Modes for coloring gene x gene set grid. Options are "none", "cutoff", and "heatmap". Required if "include_grid == TRUE". Default is "none".
#' @param prevalence_cutoff Numeric. Cutoff for coloring genes by prevalence in plotted gene sets. Required if 'prevalence_color == "cutoff"'.
#' @param chart_style String. Options are "bar" and "lollipop". Default is "bar".
#' @param include_grid Boolean. Whether or not to include the gene x gene set grid in the final plot. Default is TRUE.
#' @param lollipop_fdr_colors Numeric vector. Cutoffs for binned FDR value color groups for lollipop plots. Default is c(0.01, 0.05, 0.1, 0.2).
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' library(SEARchways)
#' library(dplyr)
#' #Run enrichment
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]))
#' enrich <- BIGprofiler(gene_list, ID="ENSEMBL", category="H")
#'
#' #Plot
#' plot_enrich2(enrich, fdr_cutoff = 0.8, gssize_col = "size_pathway")

plot_enrich2 <- function(df = NULL,
                         fdr_cutoff = 0.2,
                         pathway_col = "pathway",
                         gene_col = "genes",
                         ratio_col = "k/K",
                         fdr_col = "FDR",
                         gssize_col = "n_pathway_genes",
                         custom_genelist = NULL,
                         x_grouping_method = "geneset",
                         y_grouping_method = "fdr",
                         prevalence_color = "none",
                         prevalence_cutoff = NULL,
                         chart_style = "bar",
                         include_grid = TRUE,
                         lollipop_fdr_colors = c(0.01, 0.05, 0.1, 0.2)
){

  fdr <- size <- gs <- gssize <- ratio <- desc <- prev <- gene <- geneset <- color <- Significance <- gssize_bin <- NULL
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
  } else{
    df <- df %>%
      dplyr::filter(fdr < fdr_cutoff)
  }

  # get length of overlap for each set
  ov <- c()
  for(i in 1:nrow(df)){
    rowgenes <- df$gns[i]
    rowgenes <- unlist(rowgenes)
    ov <- c(ov, length(rowgenes))
  }
  df$overlap <- ov


  # get gene vector
  genevec <- c()
  if(!is.null(custom_genelist)){
    genevec <- custom_genelist
  }
  else{
    for(i in 1:nrow(df)){
      rowgenes <- df$gns[i]
      rowgenes <- unlist(rowgenes)
      genevec <- c(genevec,rowgenes)
    }
  }
  genevec <- unique(genevec)


  # get pathways
  gsvec <- df$gs

  # make count df
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

  if(prevalence_color == "none"){
    count_df_format <- count_df
  }

  ## format binary matrix with prevalence cutoff
  else if(prevalence_color %in% c("cutoff", "heatmap")){
    count_df_format <- count_df
    genesums <- colSums(count_df_format)
    gene_prevalence <- genesums/nrow(count_df_format) # calculate gene prevalence among selected gene sets

    if(prevalence_color == "cutoff"){
      if(is.null(prevalence_cutoff)){stop("Please provide a prevalence cutoff.")}
      for(i in 1:ncol(count_df_format)){
        col_prev <- gene_prevalence[i]
        if(col_prev < prevalence_cutoff){next} # leave column if prevalence cutoff not met
        else{ # if prev cutoff met, change 1's to 2's
          column <- count_df_format[,i]
          column[column == 1] <- 2
          count_df_format[,i] <- column
        }
      }
    } else{ # if prevalence should be colored as a heatmap, replace 1's with prevalence value
      for(i in 1:ncol(count_df_format)){
        col_prev <- gene_prevalence[i]
        column <- count_df_format[,i]
        column[column == 1] <- col_prev
        count_df_format[,i] <- column
      }
    }
  } else{stop('Valid options for heatmap coloring are "none", "cutoff", and "heatmap".')}

  ## convert to df and pivot for ggplot
  count_df_format <- as.data.frame(count_df_format) %>%
    dplyr::mutate("geneset" = rownames(count_df_format)) %>%
    tidyr::pivot_longer(genevec, names_to = "gene", values_to = "value")

  ### Get order vector for gene sets ###
  if(y_grouping_method == "hclust"){ # hierarchichal clustering by shared genes
    if(length(rownames(count_df) > 1)){
      y_cluster_df <- textshape::cluster_matrix(count_df, dim = "row")
      y_levels <- rownames(y_cluster_df)
    } else{stop("There are < 2 gene sets in your input (or above your FDR cutoff) and they cannot be clustered. Select a different option for y_grouping_method")}
  } else if(y_grouping_method == "overlap_size"){ # number of shared genes between set and plotted genes
    y_levels <- data.frame("gs" = rownames(count_df), "size" = rowSums(count_df)) %>%
      dplyr::arrange(size) %>%
      dplyr::pull(gs)
  } else if(y_grouping_method == "gs_size"){ # size of gene set
    y_levels <- df %>%
      dplyr::arrange(gssize) %>%
      dplyr::pull(gs)
  } else if(y_grouping_method == "ratio"){ # ratio of set genes in query
    y_levels <- df %>%
      dplyr::arrange(ratio) %>%
      dplyr::pull(gs)
  } else if(y_grouping_method == "fdr"){
    y_levels <- df %>%
      dplyr::arrange(desc(fdr)) %>%
      dplyr::pull(gs)
  } else if(y_grouping_method == "input"){
    y_levels <- df %>%
      dplyr::pull(gs)
  } else{stop('Valid options for y_grouping_method are "hclust", "overlap_size", "gs_size", "ratio", "fdr", and "input".')}

  titlesize <- 16
  haxis_size <- 10
  vaxis_size <- 12

  if(include_grid){
    ### Get order vector for genes ###
    if(x_grouping_method == "hclust"){ # hierarchichal clustering by gene set membership
      if(length(genevec > 1)){
        x_cluster_df <- textshape::cluster_matrix(count_df, dim = "col")
        x_levels <- colnames(x_cluster_df)
      } else{stop("There are < 2 genes in your gene list and they cannot be clustered. Select a different option for x_grouping_method")}
    } else if(x_grouping_method == "prevalence"){ # most prevalent genes first
      genesums <- colSums(count_df)
      x_levels <- data.frame("gene" = colnames(count_df), "prev" = genesums/nrow(count_df)) %>%
        dplyr::arrange(desc(prev)) %>%
        dplyr::pull("gene")
    } else if(x_grouping_method == "input"){
      x_levels <- genevec
    } else if(x_grouping_method == "alphabetical"){
      x_levels <- sort(genevec)
    } else if(x_grouping_method == "geneset"){ # biggest gene set first, prevalence within gene set. Most similar to FUMA
      genesums <- colSums(count_df)
      prev_df <- data.frame("gene" = colnames(count_df), "prev" = genesums/nrow(count_df)) %>%
        dplyr::arrange(desc(prev))

      x_levels <- c()
      count_df_gsorder <- as.data.frame(t(count_df))
      y_levels_rev <- rev(y_levels) #reverse y list so that the top one processes first (formatting quirk bc of coord_flip)

      for(i in 1:length(y_levels_rev)){
        set <- y_levels_rev[i] # get set name
        gs_genes <- count_df_gsorder[which(count_df_gsorder[[set]] == 1),] # select only rows named for genes in set (value == 1)
        row_genes <- rownames(gs_genes) # pull genes
        prev_df_row <- data.frame("gene" = row_genes) # sort by prevalence
        prev_ordered_genes <- prev_df_row %>%
          dplyr::left_join(prev_df, by = c("gene" = "gene")) %>%
          dplyr::arrange(desc(prev)) %>%
          dplyr::pull(gene)
        x_levels <- c(x_levels, prev_ordered_genes) # add to order vector

        remaining_genes <- rownames(count_df_gsorder)[which(!rownames(count_df_gsorder) %in% x_levels)] # remove used genes from binary matrix
        count_df_gsorder <- count_df_gsorder[which(rownames(count_df_gsorder) %in% remaining_genes),]
      }
      x_levels <- c(x_levels, sort(remaining_genes)) # tack any remaining genes on the end
    } else{stop('Valid options for x_grouping_method are "hclust", "prevalence", "input", "alphabetical", and "geneset".')}

    ### Create tile plot ###

    scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
    count_df_format$color <- scale_values(count_df_format$value)


    p1 <- ggplot2::ggplot(count_df_format, ggplot2::aes(x = factor(gene, levels = x_levels),
                                                        y = factor(geneset, levels = y_levels),
                                                        fill = color)) +
      ggplot2::geom_tile(color = "grey30") +

      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = haxis_size),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y =  ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank())


    # format heatmap legend
    if(prevalence_color == "none"){
      p1 <- p1 +
        ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#00418d") +
        ggplot2::theme(legend.position = "none")
    } else if(prevalence_color == "heatmap"){
      p1 <- p1 +
        ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#00418d", name = "prevalence in \nplotted gene sets")
    } else{

      if(max(gene_prevalence) > prevalence_cutoff){
        tilecols_cutoff <- c("0" = "#FFFFFF", "0.5" = "#80A0C6","1" = "#00418d")
      } else{
        tilecols_cutoff <- c("0" = "#FFFFFF", "1" = "#80A0C6")
      }

      p1 <- ggplot2::ggplot(count_df_format, ggplot2::aes(x = factor(gene, levels = x_levels),
                                                          y = factor(geneset, levels = y_levels),
                                                          fill = as.character(color))) +
        ggplot2::geom_tile(color = "grey30") +

        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1, size = haxis_size),
                       axis.text.y = ggplot2::element_blank(),
                       axis.title.y =  ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       plot.background = ggplot2::element_blank()) +
        ggplot2::scale_fill_manual(values = tilecols_cutoff, name = "prevalence in \nplotted gene sets",
                                   labels = c("not present",
                                              paste0("in < ", round(prevalence_cutoff*100, 1), "% of gene sets"),
                                              paste0("in > ", round(prevalence_cutoff*100, 1), "% of gene sets")))
    }
  }

  ### Make charts of k/K and fdr values ###
  if(chart_style == "bar"){
    # k/K bar
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x=factor(gs , levels = y_levels), y=ratio)) +
      ggplot2::geom_bar(stat = "identity", fill = "#f43545") +
      ggplot2::ggtitle("k/K") +
      #ggplot2::geom_text(aes(label = paste0(overlap, "/", gssize)), color = "black", fontface = "bold", size = 4) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                     panel.grid = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = -0.1),
                     axis.text.y = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0,0,0,0, "pt"),
                     axis.ticks.y = ggplot2::element_blank()
      )+
      ggplot2::scale_y_reverse(position = "right", labels = scales::label_number(accuracy = 0.01)) +
      ggplot2::coord_flip() #+
    # ggplot2::scale_x_discrete(position = "bottom", labels = function(gs) str_wrap(gs, width = 50))

    # FDR bar
    p4 <- ggplot2::ggplot(df, ggplot2::aes(x= factor(gs , levels = y_levels), y= -log10(fdr))) +
      ggplot2::geom_bar(stat = "identity", fill = "#fa9801") +
      ggplot2::ggtitle("-log10 FDR") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                     axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = -0.1),
                     panel.grid = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.text.y.left = ggplot2::element_blank(),
                     axis.text.y.right = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0,0,0,0, "pt")
      )+
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(position = "right",labels = scales::label_number(accuracy = 0.1))
  } else if(chart_style == "lollipop"){
    fdr_colors.sort <- sort(lollipop_fdr_colors)
    df$Significance <- NA
    fdr_levels <- c()
    for(i in 1:length(fdr_colors.sort)){
      if(i==1){
        df$Significance[df$FDR < fdr_colors.sort[i]] <- paste("FDR < ", fdr_colors.sort[i])
        fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
      } else{
        df$Significance[df$FDR < fdr_colors.sort[i] & df$FDR >= fdr_colors.sort[i-1]] <- paste("FDR < ", fdr_colors.sort[i])
        fdr_levels <- c(fdr_levels, paste("FDR < ", fdr_colors.sort[i]))
      }
    }
    df <- df %>%
      dplyr::mutate(Significance = factor(Significance, levels = fdr_levels))


    size_bins <- c("< 10" = 1, "10 - 49" = 2, "50 - 99" = 3, "100 - 499" = 4, "500 - 999" = 5, ">= 1000" = 6)
    df$gssize_bin <- NA
    df <- df %>%
      dplyr::mutate(gssize_bin = ifelse(df$n_pathway_genes < 10, "< 10",
                                        ifelse(df$n_pathway_genes < 50, "10 - 49",
                                               ifelse(df$n_pathway_genes < 100, "50 - 99",
                                                      ifelse(df$n_pathway_genes < 500, "100 - 499",
                                                             ifelse(df$n_pathway_genes < 1000, "500 - 999", ">= 1000")))))) %>%

      dplyr::mutate(gssize_bin = factor(gssize_bin, levels = c("< 10", "10 - 49", "50 - 99", "100 - 499", "500 - 999", ">= 1000")))
    p5 <- ggplot2::ggplot(df, ggplot2::aes(x = factor(gs, levels = y_levels), y = ratio)) +
      ggplot2::geom_segment(ggplot2::aes(factor(gs, levels = y_levels),
                                         xend=gs, y=0, yend=ratio)) +
      ggplot2::geom_point(ggplot2::aes(fill = Significance,
                                       size = gssize_bin),
                          shape=21, stroke=1) +
      ggplot2::scale_size_manual(values = size_bins) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::scale_fill_brewer(palette = "RdYlBu", na.value="grey70") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(size = haxis_size, angle = 45, hjust = -0.1),
                     axis.text.y = ggplot2::element_blank(),
                     #panel.background = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                     #axis.ticks.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0,0,0,1, "pt"))+
      ggplot2::scale_y_continuous(position = "right",labels = scales::label_number(accuracy = 0.1)) +
      ggplot2::ggtitle("k/K")+
      ggplot2::guides(size=ggplot2::guide_legend(title="Gene Set Size"),
                      fill=ggplot2::guide_legend(title="Significance"))

  } else(stop('Valid options chart_style are "bar" and "lollipop".'))

  # GS Size
  p3 <- ggplot2::ggplot(df, ggplot2::aes(x = factor(gs , levels = y_levels), y = rep(1, length(y_levels)))) +
    #ggplot2::geom_tile(fill = "white", width = 0.5) +
    ggplot2::ggtitle("GS\nsize") +
    ggplot2::geom_text(ggplot2::aes(label = gssize), size = 4) +
    ggplot2::coord_flip() +
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = vaxis_size),
                   #panel.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = titlesize),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(0,0,0,1, "pt"))+
    ggplot2::scale_x_discrete(position = "bottom", labels = function(gs) stringr::str_wrap(gs, width = 50))

  ### Final patchwork ###
  if(chart_style == "bar"){
    if(include_grid){
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p2 + # k/K bar
        patchwork::plot_spacer() +
        p4 + # FDR bar
        patchwork::plot_spacer() +
        p1 + # grid
        patchwork::plot_layout(widths = c(0.75, -0.2,  2 , -0.2, 2, -0.3, 25))
    } else{
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p2 + # k/K bar
        patchwork::plot_spacer() +
        p4 + # FDR bar
        patchwork::plot_layout(widths = c(0.65, -0.2,  3 , -0.2, 3))
    }
  } else if(chart_style == "lollipop"){
    if(include_grid){
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p5 + # lollipop
        patchwork::plot_spacer() +
        p1 + # grid
        patchwork::plot_layout(widths = c(0.75, -0.2, 4, -0.3, 25), guides = "collect")
    } else{
      p <- p3 + # GS size
        patchwork::plot_spacer() +
        p5 + # lollipop
        patchwork::plot_layout(widths = c(0.65, -0.2,  6))
    }
  } else(stop('Valid options chart_style are "bar" and "lollipop".'))

  return(p)
}
