#' Plot GSEA normalized enrichment
#'
#' @param gsea Data frame output by SEARchways::BIGsea including pathway, FDR, and NES
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
#' #Get fold change information from example model
#' genes.FC <- example_model$lmerel %>%
#'             filter(variable == "virus") %>%
#'             select(variable, hgnc_symbol, estimate)
#' #Run GSEA
#' example_gsea <- BIGsea(gene_df = genes.FC, category = "H")
#'
#' #Plot
#' plot_gsea(example_gsea)

plot_gsea <- function(gsea, fdr.cutoff = 0.2,
                      fdr.colors = c(0.01, 0.05, 0.1, 0.2),
                      show.overlap = TRUE
                      ){
  FDR <-NES<-Significance<-pathway<- NULL
  #### Format data ####
  dat.signif <- gsea %>%
    dplyr::filter(FDR < fdr.cutoff)
  #keep nonsignif overlap if requested
  if(show.overlap){
    dat.format <- gsea %>%
      dplyr::filter(pathway %in% dat.signif$pathway) %>%
      dplyr::mutate(pathway = gsub("_", " ", pathway)) %>%
      #Add 0 enrichment values
      tidyr::complete(group, pathway) %>%
      dplyr::mutate(NES = ifelse(is.na(NES), 0, NES),
                    FDR = ifelse(is.na(FDR), 1, FDR))
  } else{
    dat.format <- dat.signif %>%
      dplyr::mutate(pathway = gsub("_", " ", pathway))
  }

  dat.format <- gsea %>%
    dplyr::mutate(pathway = gsub("_", " ", pathway)) %>%
    dplyr::filter(FDR < fdr.cutoff)

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

  #Enrichment score limits
  plot.lim <- max(abs(dat.format$NES))+0.1

  #### Plot ####
  p1 <- dat.format %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(pathway, NES), NES)) +
    ggplot2::geom_segment(ggplot2::aes(stats::reorder(pathway, NES),
                     xend=pathway, y=0, yend=NES)) +
    ggplot2::geom_point(size=3, ggplot2::aes(fill = Significance),
               shape=21, stroke=1) +
    ggplot2::geom_hline(yintercept = 0) +

    ggplot2::scale_fill_brewer(palette = "RdYlBu", na.value="grey70") +
    ggplot2::lims(y=c(-plot.lim,plot.lim)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x="", y="Normalized Enrichment Score") +
    ggplot2::theme_bw()

  if(show.overlap){
    p1.facet <- p1 + ggplot2::facet_grid( ~ group)
  } else {
    p1.facet <- p1 + ggplot2::facet_wrap( ~ group, scales="free")
  }
  return(p1.facet)
}



