#' Plot STRING network of genes colored by enriched pathways
#'
#' To be used in conjunction with map_string
#'
#' @param map List output by map_string
#' @param discard Character string identifying genes to remove from the network. Can be "none" (Default), "orphan" (genes with 0 connections in the network), or "orphan.unenrich" (genes with 0 connections and no significant enrichment colors)
#' @param layout Character string for network layout algorithm. See options in igraph::layout_with_
#' @param enrichment Data frame output by hypergeometric enrichment
#' @param overlap Numeric minimum of significant genes in terms to be used as colors
#' @param FDR Numeric maximum FDR of terms to be used as colors
#' @param ID Character string defining gene ID column in enrichment that matches genes used in map_string. Can be "SYMBOLs" (Default), "ENTREZIDs", "ENSEMBLIDs"
#' @param colors Character vector of custom colors to use. Must be at least a long as total significant terms plus 1 for the "none" group
#' @param text_size Numeric size of gene labels on network nodes. Default of 2
#' @param node_size Numeric size of network nodes. Default of 1
#'
#' @return ggplot STRING network object
#' @export
#'
#' @examples
#' map <- map_string(genes = c("WNT5B","CCND2","FZD1","CTNNB1","DKK1","TCF7"),
#'                        version = 11.5, score_threshold = 400)
#' plot_string(map)
#' plot_string(map, enrichment = BIGpicture::enrichment, ID ="SYMBOLs")

plot_string <- function(map, discard="none", layout='fr',
                        enrichment=NULL, overlap=2, FDR=0.2,
                        ID=c("SYMBOLs","ENTREZIDs","ENSEMBLIDs"),
                        colors=NULL, text_size=2, node_size=1){
  Description <- STRING_id <- combined_score <- gene <- none <- total <- value <- size.overlap.term <- p.adjust <- NULL

  #### Format enrichment colors ####
  if(!is.null(enrichment)){
    #Get significant enrichments
    col.mat <- enrichment %>%
      dplyr::ungroup() %>%
      dplyr::filter(size.overlap.term >= overlap & p.adjust < FDR) %>%
      dplyr::select(Description, dplyr::all_of(ID))
    names(col.mat)[names(col.mat) == ID] <- 'gene'

    #Error if no terms to plot
    if(nrow(col.mat) == 0) {stop("No significant enrichment terms.
                               Try increasing FDR.")}

    #Format enrichment results for scatterpie plotting
    col.mat.format <- col.mat %>%
      #Split gene lists within terms
      tidyr::unnest(gene) %>%
      dplyr::distinct(gene, Description) %>%
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
      dplyr::arrange(match(Description, unique(enrichment$Description))) %>%
      tidyr::pivot_wider(names_from = Description, values_fill = 0) %>%
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

  #### Discard unconnected if specified ####
  ##All nodes
  nodes <- which(igraph::degree(map[["subgraph"]])>=0)
  ##Nodes without edge connections
  isolated <- which(igraph::degree(map[["subgraph"]])==0)

  ##Remove unconnected regardless of enrichment
  if(discard == "orphan"){
    subgraph.filter <- igraph::delete.vertices(map[["subgraph"]], isolated)
  } else if(discard == "orphan.unenrich"){
    ##Nodes without enrichment
    unenrich <- map.unique %>%
      dplyr::filter(none==1) %>%
      dplyr::distinct(STRING_id) %>%
      dplyr::pull(STRING_id)
    ##Remove unconnected that are also unenriched
    isolated.unenrich <- isolated[names(isolated) %in% unenrich]
    subgraph.filter <-  igraph::delete.vertices(map[["subgraph"]], isolated.unenrich)
  } else if(discard == "none") {
    subgraph.filter <- map[["subgraph"]]
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
  if(!is.null(enrichment)){
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
                                         fill=NULL), size = node_size*10, color=color.vec) +
      ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter)$x,
                                            y = igraph::V(subgraph.filter)$y,
                                            label = igraph::V(subgraph.filter)$symbol),
                               size=text_size) +
      ggnetwork::theme_blank() + ggplot2::coord_fixed()
  }
  return(plot.col)
}

