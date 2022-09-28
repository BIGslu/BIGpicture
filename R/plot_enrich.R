#' Plot k/K for hypergeo enrichment
#'
#' @param enrich Data frame output by SEARchways::BIGprofiler or SEARchways::BIGenrichr
#' @param fdr.cutoff Numeric. Maximum FDR to plot. Default is 0.2
#' @param fdr.colors Numeric vector. Cutoffs for color groups. Default is c(0.01, 0.05, 0.1, 0.2)
#' @param show.overlap Logical if should show overlap across all facets even if some missing (TRUE) or give each facet it's own axis labels (FALSE). Default is TRUE
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' library(SEARchways)
#' library(dplyr)
#' #Run enrichment
#' gene_list <- list(HRV1 = names(example_gene_list[[1]]),
#'                   HRV2 = names(example_gene_list[[2]]))
#' enrich <- BIGprofiler(gene_list, ID="ENSEMBL", category="H")
#'
#' #Plot
#' plot_enrich(enrich, fdr.cutoff = 0.5, fdr.colors = c(0.05, 0.5))

plot_enrich <- function(enrich, fdr.cutoff = 0.2,
                      fdr.colors = c(0.01, 0.05, 0.1, 0.2),
                      show.overlap = TRUE){
  FDR <-`k/K`<-Significance<-pathway<- NULL

  #### Format data ####
  dat.signif <- enrich %>%
    dplyr::filter(FDR < fdr.cutoff)
  #keep nonsignif overlap if requested
  if(show.overlap){
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

  if(nrow(dat.format) == 0){stop("No gene sets are significant. Please increase fdr.cutoff.")}

  fdr.colors.sort <- sort(fdr.colors)
  dat.format$Significance <- NA
  for(i in 1:length(fdr.colors.sort)){
    if(i==1){
      dat.format$Significance[dat.format$FDR < fdr.colors.sort[i]] <- paste("FDR <", fdr.colors.sort[i])
    } else{
      dat.format$Significance[dat.format$FDR < fdr.colors.sort[i] & dat.format$FDR >= fdr.colors.sort[i-1]] <- paste("FDR <", fdr.colors.sort[i])
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

  if(show.overlap){
    p1.facet <- p1 + ggplot2::facet_grid( ~ group)
  } else {
    p1.facet <- p1 + ggplot2::facet_wrap( ~ group, scales="free")
  }
  return(p1.facet)
}



