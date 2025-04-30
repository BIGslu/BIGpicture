#' Plot STRING network of genes colored by enriched pathways
#'
#' To be used in conjunction with map_string
#'
#' @param map List output by map_string
#' @param layout Character string for network layout algorithm. See options in igraph::layout_with_
#' @param edge_min Numeric minimum edges a node must have to be displayed. Default is 0 meaning orphan nodes are included
#' @param edge_max Numeric maximum edges a node must have to be displayed. Default in Inf. Set to 0 to see only orphan nodes
#' @param main_cluster_only Logical if should include only genes connected to the largest cluster
#' @param enriched_only Logical if should include only genes in significantly enriched terms. Default FALSE
#' @param enrichment Data frame output by `BIGprofiler`, `flexEnrich`,  `BIGsea`. For use in coloring nodes
#' @param overlap Numeric minimum of total significant genes in a enrichment term to be used as colors (`BIGprofiler`, `flexEnrich`)
#' @param fdr_cutoff Numeric maximum FDR of enrichment terms to be used as colors (`BIGprofiler`, `flexEnrich`, `BIGsea`)
#' @param colors Character vector of custom colors to use. Must be at least a long as total significant terms plus 1 for the "none" group
#' @param text_size Numeric size of gene labels on network nodes. Default of 2
#' @param node_size Numeric size of network nodes. Default of 1
#'
#' @return ggplot STRING network object
#' @export
#'
#' @examples
#' #NOT RUN
#' # map <- map_string(genes = c("MTHFD2","ISOC1","IL2RB","SAMHD1","NAMPT","NOD1",
#' #                            "NFKB1","IFIT3","ZBP1","IL15RA","SP110","ITGB7",
#' #                            "SERPING1","B2M","CXCL11","USP18","MAPKAPK2","DKK1"),
#' #                version = 12, score_threshold = 400)
#' # plot_string(map)
#' #
#' # Add enrichment colors
#' # library(dplyr)
#' # library(SEARchways)
#' # genes.OI <- example.model$lmerel %>%
#' #             filter(variable == "virus" & FDR < 0.2) %>%
#' #             select(variable, gene)
#' # example_enrich <- SEARchways::BIGprofiler(gene_df = genes.OI,
#' #                                           category = "H", ID = "ENSEMBL")
#' #
#' # map2 <- map_string(genes = genes.OI$gene,
#' #                  version = 11.5, score_threshold = 400)
#' # plot_string(map2, enrichment = example_enrich, fdr_cutoff=0.2)
#' # plot_string(map2, enrichment = example_enrich, fdr_cutoff=0.2, enriched_only=TRUE)
#' #
#' # Add GSEA colors
#' # genes.FC <- example.model$lmerel %>%
#' #             filter(variable == "virus") %>%
#' #             select(variable, gene, estimate)
#' # example_gsea <- BIGsea(gene_df = genes.FC, category = "H", ID = "ENSEMBL")
#' #
#' # plot_string(map2, enrichment = example_gsea, fdr_cutoff = 0.3,
#' #             edge_max = 0, enriched_only=TRUE)

plot_string <- function(map, layout='fr',
                        edge_min=0, edge_max=Inf,
                        main_cluster_only=FALSE,
                        enriched_only = FALSE,
                        enrichment=NULL, overlap=2, fdr_cutoff=0.2,
                        colors=NULL, text_size=2, node_size=1){
  pathway <- STRING_id <- combined_score <- gene <- none <- total <- value <- group_in_pathway <- FDR <- genes <- leadingEdge <- legend.title <- NULL

  #### Format enrichment colors ####
  if(!is.null(enrichment)){
    if("k/K" %in% colnames(enrichment)){
    #Get significant enrichments
    col.mat <- enrichment %>%
      dplyr::ungroup() %>%
      dplyr::filter(group_in_pathway >= overlap & FDR <= fdr_cutoff) %>%
      dplyr::select(pathway, genes) %>%
      dplyr::rename(gene=genes)
    legend.title <- "Enriched pathways"
    }
    if("NES" %in% colnames(enrichment)){
      #Get significant GSEA
      col.mat <- enrichment %>%
        dplyr::ungroup() %>%
        dplyr::filter(FDR <= fdr_cutoff) %>%
        dplyr::select(pathway, leadingEdge) %>%
        dplyr::rename(gene=leadingEdge)
      legend.title <- "GSEA leading edge"
    }

    #Error if no terms to plot
    if(nrow(col.mat) == 0) {stop("No significant enrichment/GSEA terms.
                               Try increasing fdr_cutoff.")}

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

  ##Remove nodes not in an enriched pathway (leave only colored nodes)
  ## all nodes currently
  if(enriched_only){
    ##Nodes without enrichment
    unenrich <- map.unique %>%
      dplyr::filter(none==1 &
                      STRING_id %in% igraph::vertex_attr(map[["subgraph"]])$name) %>%
      dplyr::distinct(STRING_id) %>%
      dplyr::pull(STRING_id)
    ##Remove unconnected that are also unenriched
    subgraph.filter <-  igraph::delete.vertices(map[["subgraph"]], unenrich)
  } else {
    subgraph.filter <- map[["subgraph"]]
  }

  if(length(igraph::vertex_attr(subgraph.filter)$name)==0){
    stop("No genes in network with significant enrichment. Consider increasing fdr_cutoff or setting enriched_only=FALSE")
  }
  #### Filter nodes by number of edges ####
  if(edge_max == edge_min){
    nodes_to_remove <- which(igraph::degree(subgraph.filter)!=edge_min)
  } else {
    nodes_to_remove <- which(igraph::degree(subgraph.filter)<edge_min |
                               igraph::degree(subgraph.filter)>edge_max)
  }

  if(length(nodes_to_remove) > 0){
    subgraph.filter2 <- igraph::delete.vertices(subgraph.filter, nodes_to_remove)
  } else {
    subgraph.filter2 <- subgraph.filter
  }

  if(length(igraph::vertex_attr(subgraph.filter2)$name)==0){
    stop("No genes in network remain after edge filtering. Consider changind edge_min and/or edge_max")
  }

  #### Filter largest cluster if selected ####
  if(main_cluster_only){
    # identify connected components
    comps <- igraph::components(subgraph.filter2)
    # largest cluster
    largest_comp_id <- which.max(comps$csize)
    # Get the vertex names in the largest component
    vertices_in_largest <- igraph::V(subgraph.filter2)[comps$membership == largest_comp_id]
    # Filter subgraph
    subgraph.filter2 <- igraph::induced_subgraph(
      subgraph.filter2,
      vids = vertices_in_largest)
  }

  #### Arrange metadata as in network ####
  map.arrange <- map.unique %>%
    dplyr::filter(STRING_id %in% igraph::vertex_attr(subgraph.filter2)$name) %>%
    dplyr::arrange(match(STRING_id, c(igraph::vertex_attr(subgraph.filter2)$name)))

  # Set attributes
  ## Check order first
  if(!identical(igraph::vertex_attr(subgraph.filter2)$name, map.arrange$STRING_id)){
    stop("igraph gene order does not match color information.")
  }

  ##gene names
  igraph::V(subgraph.filter2)$symbol <- map.arrange$gene
  ##enrichment colors
  for(term in colnames(map.arrange)[-c(1:3)]){
    igraph::vertex_attr(subgraph.filter2)[[term]] <- unlist(map.arrange[term])
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
  #### set layout ####
  layout_funs <- list(
    nice     = igraph::layout_nicely,
    fr       = igraph::layout_with_fr,
    bipar    = igraph::layout_as_bipartite,
    star     = igraph::layout_as_star,
    tree     = igraph::layout_as_tree,
    circle   = igraph::layout_in_circle,
    kk       = igraph::layout_with_kk,
    graphopt = igraph::layout_with_graphopt,
    gem      = igraph::layout_with_gem,
    dh       = igraph::layout_with_dh,
    sphere   = igraph::layout_on_sphere,
    grid     = igraph::layout_on_grid,
    lgl      = igraph::layout_with_lgl,
    mds      = igraph::layout_with_mds,
    sugi     = igraph::layout_with_sugiyama
  )

  if (!layout %in% names(layout_funs)) {
    stop("Unknown layout: ", layout)
  }
  set.seed(8434)
  xy <- layout_funs[[layout]](subgraph.filter2)

  #### plot ####
  igraph::V(subgraph.filter2)$x <- xy[, 1]
  igraph::V(subgraph.filter2)$y <- xy[, 2]

  plot <- ggraph::ggraph(subgraph.filter2, layout= "manual",
                         x = igraph::V(subgraph.filter2)$x,
                         y = igraph::V(subgraph.filter2)$y) +
    #Edges
    ggraph::geom_edge_link(ggplot2::aes(width=combined_score), color="grey70") +
    ggraph::scale_edge_width(range = c(0.2,2), name="STRING score")

  #Add nodes
  suppressWarnings(
    if(!is.null(enrichment)){
      if(length(colnames(map.arrange)[-c(1:3)])>1){
        plot.col <- plot +
          scatterpie::geom_scatterpie(data=igraph::as_data_frame(subgraph.filter2,
                                                                 "vertices"),
                                      cols=colnames(map.arrange)[-c(1:3)], color=NA,
                                      pie_scale = node_size) +
          ggplot2::scale_fill_manual(values = color.vec, name = legend.title) +
          ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter2)$x,
                                                y = igraph::V(subgraph.filter2)$y,
                                                label=igraph::V(subgraph.filter2)$symbol),
                                   size=text_size) +
          ggnetwork::theme_blank() + ggplot2::coord_fixed()
      } else{
        plot.col <- plot +
          ggnetwork::geom_nodes(ggplot2::aes(x = igraph::V(subgraph.filter2)$x,
                                             y = igraph::V(subgraph.filter2)$y,
                                             fill=NULL,
                                             color=colnames(map.arrange)[-c(1:3)]),
                                size = node_size*10) +
          ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter2)$x,
                                                y = igraph::V(subgraph.filter2)$y,
                                                label = igraph::V(subgraph.filter2)$symbol),
                                   size=text_size) +
          ggnetwork::theme_blank() + ggplot2::coord_fixed() +
          ggplot2::scale_color_manual(values = color.vec, name = legend.title)
      }
    } else{
      plot.col <- plot +
        ggnetwork::geom_nodes(ggplot2::aes(x = igraph::V(subgraph.filter2)$x,
                                           y = igraph::V(subgraph.filter2)$y,
                                           fill=NULL), size = node_size*10,
                              color=color.vec) +
        ggnetwork::geom_nodetext(ggplot2::aes(x = igraph::V(subgraph.filter2)$x,
                                              y = igraph::V(subgraph.filter2)$y,
                                              label = igraph::V(subgraph.filter2)$symbol),
                                 size=text_size) +
        ggnetwork::theme_blank() + ggplot2::coord_fixed()
    }
  )
  return(plot.col)
}

