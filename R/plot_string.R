#' Plot STRING network of genes colored by enriched pathways
#'
#' To be used in conjunction with map_string
#'
#' @param map List output by map_string
#' @param layout Character string for network layout algorithm. See options in igraph::layout_with_
#' @param discard Character string identifying genes to remove from the network. Can be "none" (Default), "orphan" (genes with 0 connections in the network), or "cluster" (genes with any connections, leaving only orphans)
#' @param enriched.only Logical if should include only genes in significantly enriched terms. Default FALSE
#' @param enrichment Data frame output by `BIGprofiler`, `BIGenrichr`, or `BIGsea`. For use in coloring nodes
#' @param overlap Numeric minimum of total significant genes in a enrichment term to be used as colors (`BIGprofiler`, `BIGenrichr`)
#' @param fdr.cutoff Numeric maximum FDR of enrichment terms to be used as colors (`BIGprofiler`, `BIGenrichr`, `BIGsea`)
#' @param colors Character vector of custom colors to use. Must be at least a long as total significant terms plus 1 for the "none" group
#' @param text_size Numeric size of gene labels on network nodes. Default of 2
#' @param node_size Numeric size of network nodes. Default of 1
#'
#' @return ggplot STRING network object
#' @export
#'
#' @examples
#' map <- map_string(genes = c("MTHFD2","ISOC1","IL2RB","SAMHD1","NAMPT","NOD1",
#'                             "NFKB1","IFIT3","ZBP1","IL15RA","SP110","ITGB7",
#'                             "SERPING1","B2M","CXCL11","USP18","MAPKAPK2","DKK1"),
#'                  version = 11.5, score_threshold = 400)
#' plot_string(map)
#'
#' # Add enrichment colors
#' library(dplyr)
#' library(SEARchways)
#' genes.OI <- example_model$lmerel %>%
#'             filter(variable == "virus" & FDR < 0.2) %>%
#'             distinct(variable, hgnc_symbol)
#' example_enrich <- BIGprofiler(gene_df = genes.OI, category = "H")
#'
#' plot_string(map, enrichment = example_enrich, fdr.cutoff=0.2)
#' plot_string(map, enrichment = example_enrich, fdr.cutoff=0.2, enriched.only=TRUE)
#'
#' # Add GSEA colors
#' genes.FC <- example_model$lmerel %>%
#'             filter(variable == "virus") %>%
#'             select(variable, hgnc_symbol, estimate)
#' example_gsea <- BIGsea(gene_df = genes.FC, category = "H")
#'
#' plot_string(map, enrichment = example_gsea, fdr.cutoff=0.2, discard = "cluster")

plot_string <- function(map, layout='fr',
                        discard="none", enriched.only = FALSE,
                        enrichment=NULL, overlap=2, fdr.cutoff=0.2,
                        colors=NULL, text_size=2, node_size=1){
  pathway <- STRING_id <- combined_score <- gene <- none <- total <- value <- group_in_pathway <- FDR <- genes <- leadingEdge <- NULL

  #### Format enrichment colors ####
  if(!is.null(enrichment)){
    if("k/K" %in% colnames(enrichment)){
    #Get significant enrichments
    col.mat <- enrichment %>%
      dplyr::ungroup() %>%
      dplyr::filter(group_in_pathway >= overlap & FDR <= fdr.cutoff) %>%
      dplyr::select(pathway, genes) %>%
      dplyr::rename(gene=genes)
    }
    if("NES" %in% colnames(enrichment)){
      #Get significant GSEA
      col.mat <- enrichment %>%
        dplyr::ungroup() %>%
        dplyr::filter(FDR <= fdr.cutoff) %>%
        dplyr::select(pathway, leadingEdge) %>%
        dplyr::rename(gene=leadingEdge)
    }

    #Error if no terms to plot
    if(nrow(col.mat) == 0) {stop("No significant enrichment terms.
                               Try increasing fdr.cutoff.")}

    #Format enrichment results for scatterpie plotting
    col.mat.format <- col.mat %>%
      #Split gene lists within terms
      tidyr::unnest(gene) %>%
      dplyr::distinct(gene, pathway) %>%
      dplyr::mutate(value=1) %>%
      #Add string ID
      dplyr::inner_join(map[["map"]], by = "gene") %>%
      dplyr::select(-gene) %>%
      dplyr::distinct() %>%
      #Calculate terms per ID
      dplyr::group_by(STRING_id) %>%
      dplyr::mutate(total = sum(value)) %>%
      dplyr::ungroup() %>%
      #Calculate proportions within terms
      dplyr::arrange(match(pathway, unique(enrichment$pathway))) %>%
      tidyr::pivot_wider(names_from = pathway, values_fill = 0) %>%
      dplyr::mutate(dplyr::across(-c(STRING_id,total), ~./total))

    #Add to STRING data and create dummy group for genes without enrichment
    map.unique <- map[["map"]] %>%
      #collapse gene names
      dplyr::group_by(STRING_id) %>%
      dplyr::mutate(gene = paste(unique(gene), collapse=" / ")) %>%
      #add enrichment color info
      dplyr::left_join(col.mat.format, by = "STRING_id") %>%
      dplyr::distinct() %>%
      #no enrichment group
      dplyr::mutate(none = ifelse(is.na(total),1,0)) %>%
      dplyr::ungroup() %>%
      #fill in 0
      dplyr::mutate_if(is.numeric, ~tidyr::replace_na(., 0))
    #Remove none is not used, remove it
    if(sum(map.unique$none)==0){
      map.unique <- map.unique %>%
        dplyr::select(-none)
    }
  } else{
    map.unique <- map[["map"]] %>%
      #collapse gene names
      dplyr::group_by(STRING_id) %>%
      dplyr::mutate(gene = paste(unique(gene), collapse=" / ")) %>%
      dplyr::ungroup()
  }

  #### Discard nodes if specified ####
  ##All nodes
  # nodes <- which(igraph::degree(map[["subgraph"]])>=0)
  ##Nodes without edge connections
  orphans <- which(igraph::degree(map[["subgraph"]])==0)
  ##Nodes with 1+ edge connections
  clusters <- which(igraph::degree(map[["subgraph"]])>0)

  ##Remove unconnected regardless of enrichment
  if(discard == "orphan"){
    subgraph.filter <- igraph::delete.vertices(map[["subgraph"]], orphans)
  } else if (discard == "cluster"){
    subgraph.filter <- igraph::delete.vertices(map[["subgraph"]], clusters)
  }else if(discard == "none") {
    subgraph.filter <- map[["subgraph"]]
  }

  ##Remove nodes not in an enriched pathway (leave only colored nodes)
  ## all nodes currently
  curr.nodes <- which(igraph::degree(subgraph.filter)>=0)
  if(enriched.only){
    ##Nodes without enrichment
    unenrich <- map.unique %>%
      dplyr::filter(none==1) %>%
      dplyr::distinct(STRING_id) %>%
      dplyr::pull(STRING_id)
    ##Remove unconnected that are also unenriched
    unenrich.nodes <- curr.nodes[names(curr.nodes) %in% unenrich]
    subgraph.filter <-  igraph::delete.vertices(subgraph.filter, unenrich.nodes)
  }

  #### Arrange metadata as in network ####
  map.arrange <- map.unique %>%
    dplyr::filter(STRING_id %in% igraph::vertex_attr(subgraph.filter)$name) %>%
    dplyr::arrange(match(STRING_id, c(igraph::vertex_attr(subgraph.filter)$name)))

  # Set attributes
  ## Check order first
  if(!identical(igraph::vertex_attr(subgraph.filter)$name, map.arrange$STRING_id)){
    stop("igraph gene order does not match color information.")
  }

  ##gene names
  igraph::V(subgraph.filter)$symbol <- map.arrange$gene
  ##enrichment colors
  for(term in colnames(map.arrange)[-c(1:3)]){
    igraph::vertex_attr(subgraph.filter)[[term]] <- unlist(map.arrange[term])
  }

  #### Set color values ####
  if(!is.null(colors)){
    color.vec <- colors
  } else if(!is.null(enrichment)){
    if("none" %in% colnames(map.arrange)){
      #ggplot colors for sets in plot
      color <- scales::hue_pal()(ncol(col.mat.format)-2)
      #Replace none with grey
      all.term <- sort(colnames(map.arrange)[-c(1:3)])
      none.index <- match("none", all.term)
      color.vec <- c(color[1:none.index-1], "grey70",
                     color[none.index:length(color)])
      color.vec <- color.vec[!is.na(color.vec)]
    } else {
      color.vec <- scales::hue_pal()(ncol(col.mat.format)-2)
    }
  } else if (is.null(enrichment) & is.null(colors)){
    color.vec <- "#d9d9d9"
  }

  #### Plot ####
  #Get xy of nodes for manual layout
  set.seed(8434)
  #### set layout ####
  if(layout == "fr"){ xy <- igraph::layout_with_fr(subgraph.filter) } else
    if(layout == "bipar"){ xy <- igraph::layout_as_bipartite(subgraph.filter) } else
      if(layout == "star"){ xy <- igraph::layout_as_star(subgraph.filter) } else
        if(layout == "tree"){ xy <- igraph::layout_as_tree(subgraph.filter) } else
          if(layout == "circle"){ xy <- igraph::layout_in_circle(subgraph.filter) } else
            if(layout == "kk"){ xy <- igraph::layout_with_kk(subgraph.filter) } else
              if(layout == "graphopt"){ xy <- igraph::layout_with_graphopt(subgraph.filter) } else
                if(layout == "gem"){ xy <- igraph::layout_with_gem(subgraph.filter) } else
                  if(layout == "dh"){ xy <- igraph::layout_with_dh(subgraph.filter) } else
                    if(layout == "sphere"){ xy <- igraph::layout_on_sphere(subgraph.filter) } else
                      if(layout == "grid"){ xy <- igraph::layout_on_grid(subgraph.filter) } else
                        if(layout == "lgl"){ xy <- igraph::layout_with_lgl(subgraph.filter) } else
                          if(layout == "mds"){ xy <- igraph::layout_with_mds(subgraph.filter) } else
                            if(layout == "sugi"){ xy <- igraph::layout_with_sugiyama(subgraph.filter) }

  #### plot ####
  igraph::V(subgraph.filter)$x <- xy[, 1]
  igraph::V(subgraph.filter)$y <- xy[, 2]

  plot <- ggraph::ggraph(subgraph.filter, layout= "manual",
                         x = igraph::V(subgraph.filter)$x,
                         y = igraph::V(subgraph.filter)$y) +
    #Edges
    ggraph::geom_edge_link(ggplot2::aes(width=combined_score), color="grey70") +
    ggraph::scale_edge_width(range = c(0.2,2), name="STRING score")

  #Add nodes
  suppressWarnings(
    if(!is.null(enrichment)){
      if(length(colnames(map.arrange)[-c(1:3)])>1){
        plot.col <- plot +
          scatterpie::geom_scatterpie(data=igraph::as_data_frame(subgraph.filter,
                                                                 "vertices"),
                                      cols=colnames(map.arrange)[-c(1:3)], color=NA,
                                      pie_scale = node_size) +
          ggplot2::scale_fill_manual(values = color.vec, name = "") +
          ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                                y = igraph::V(subgraph.filter)$y,
                                                label=igraph::V(subgraph.filter)$symbol),
                                   size=text_size) +
          ggnetwork::theme_blank() + ggplot2::coord_fixed()
      } else{
        plot.col <- plot +
          ggnetwork::geom_nodes(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                             y = igraph::V(subgraph.filter)$y,
                                             fill=NULL,
                                             color=colnames(map.arrange)[-c(1:3)]),
                                size = node_size*10) +
          ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                                y = igraph::V(subgraph.filter)$y,
                                                label = igraph::V(subgraph.filter)$symbol),
                                   size=text_size) +
          ggnetwork::theme_blank() + ggplot2::coord_fixed() +
          ggplot2::scale_color_manual(values = color.vec, name = "")
      }
    } else{
      plot.col <- plot +
        ggnetwork::geom_nodes(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                           y = igraph::V(subgraph.filter)$y,
                                           fill=NULL), size = node_size*10,
                              color=color.vec) +
        ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                              y = igraph::V(subgraph.filter)$y,
                                              label = igraph::V(subgraph.filter)$symbol),
                                 size=text_size) +
        ggnetwork::theme_blank() + ggplot2::coord_fixed()
    }
  )
  return(plot.col)
}

