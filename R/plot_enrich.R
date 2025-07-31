#' Plot k/K for hypergeo enrichment
#'
#' @param enrich Data frame output by SEARchways::BIGprofiler or SEARchways::BIGenrichr
#' @param fdr_cutoff Numeric. Maximum FDR to plot. Default is 0.2
#' @param fdr_colors Numeric vector. Cutoffs for color groups. Default is c(0.01, 0.05, 0.1, 0.2)
#' @param show_overlap Logical if should show overlap across all facets even if some missing (TRUE) or give each facet it's own axis labels (FALSE). Default is TRUE
#'
#' @param fdr.cutoff Deprecated form of fdr_cutoff
#' @param fdr.colors Deprecated form of fdr_colors
#' @param show.overlap Deprecated form of show_overlap
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' library(SEARchways)
#' library(dplyr)
#' #Run enrichment
#' gene_list <- list(HRV1 = names(example.gene.list[[1]]),
#'                   HRV2 = names(example.gene.list[[2]]))
#' enrich <- flexEnrich(gene_list, ID="ENSEMBL", collection="H")
#'
#' #Plot
#' plot_enrich(enrich, fdr_cutoff = 0.5, fdr_colors = c(0.05, 0.5))

plot_enrich <- function(enrich,
                        fdr_cutoff = 0.2,
                        fdr_colors = c(0.01, 0.05, 0.1, 0.2),
                        show_overlap = TRUE,
                        #Deprecated
                        fdr.cutoff = NULL, fdr.colors = NULL, show.overlap = NULL){
  FDR <-`k/K`<-Significance<-pathway<-group<- NULL

  # backwards compatibility
  if(!is.null(fdr.cutoff)){fdr_cutoff <- fdr.cutoff}
  if(!is.null(fdr.colors)){fdr_colors <- fdr.colors}
  if(!is.null(show.overlap)){show_overlap <- show.overlap}

  #### Format data ####
  dat.signif <- enrich %>%
    dplyr::filter(FDR < fdr_cutoff)
  #keep nonsignif overlap if requested
  if(show_overlap){
    dat.format <- enrich %>%
      dplyr::filter(pathway %in% dat.signif$pathway) %>%
      dplyr::mutate(pathway = gsub("_", " ", pathway)) %>%
      #Add 0 enrichment values
      tidyr::complete(group, pathway) %>%
      dplyr::mutate(`k/K` = ifelse(is.na(`k/K`), 0, `k/K`),
                    FDR = ifelse(is.na(FDR), 1, FDR))
  } else{
    dat.format <- dat.signif %>%
      dplyr::mutate(pathway = gsub("_", " ", pathway))
  }

  if(nrow(dat.format) == 0){stop("No gene sets are significant. Please increase fdr_cutoff.")}

  fdr_colors.sort <- sort(fdr_colors)
  dat.format$Significance <- NA
  for(i in 1:length(fdr_colors.sort)){
    if(i==1){
      dat.format$Significance[dat.format$FDR < fdr_colors.sort[i]] <- paste("FDR <", fdr_colors.sort[i])
    } else{
      dat.format$Significance[dat.format$FDR < fdr_colors.sort[i] & dat.format$FDR >= fdr_colors.sort[i-1]] <- paste("FDR <", fdr_colors.sort[i])
    }
  }

  #### Plot ####
  p1 <- dat.format %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(pathway, `k/K`), `k/K`)) +
    ggplot2::geom_segment(ggplot2::aes(stats::reorder(pathway, `k/K`),
                                       xend=pathway, y=0, yend=`k/K`)) +
    ggplot2::geom_point(size=3, ggplot2::aes(fill = Significance),
                        shape=21, stroke=1) +
    ggplot2::geom_hline(yintercept = 0) +

    ggplot2::scale_fill_brewer(palette = "RdYlBu", na.value="grey70") +
    ggplot2::coord_flip() +
    ggplot2::labs(x="", y="Proportion enriched (k / K)") +
    ggplot2::theme_bw()

  if(show_overlap & length(unique(dat.format$group)) > 1){
    p1.facet <- p1 + ggplot2::facet_grid( ~ group)
  } else if(length(unique(dat.format$group)) > 1){
    p1.facet <- p1 + ggplot2::facet_wrap( ~ group, scales="free")
  } else{
    p1.facet <- p1
  }
  return(p1.facet)
}



