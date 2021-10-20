#' Venn diagrams of significant genes
#'
#' @param model_result List of data frames output by kimma::kmFit()
#' @param model Character string of model to plot. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param variables Character vector of variables in model_result to include in plots
#' @param contrasts Character vector of contrasts in model_result to include in plots. Only applicable if model name includes 'contrast'
#' @param intercept Logical if should include the intercept variable
#' @param fdr.cutoff Numeric vector of FDR cutoffs to assess. One venn per FDR value
#'
#' @return List of ggplot objects
#' @export
#'
#' @examples
#' plot_venn_genes(model_result, model = "lme", fdr.cutoff = c(0.05,0.5))

plot_venn_genes <- function(model_result, model, variables=NULL, intercept=FALSE,
                            contrasts=NULL,
                            fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5)){

  FDR <- variable <- gene <- contrast <- NULL

  #Extract results
  dat <- model_result[[model]]

  #List variables of interest
  var_all <- unique(dat$variable)

  if(is.null(variables) & intercept){
    var_filter <- var_all
  } else if (is.null(variables)){
    var_filter <- var_all[var_all != "(Intercept)"]
  } else if (!is.null(variables) & intercept) {
    var_filter <- unique(c("(Intercept)", variables))
  } else{
    var_filter <- variables
  }

  #List contrasts of interest
  if(grepl("contrast", model)){
    if(is.null(contrasts)){
      con_filter <- unique(dat$contrast)
      } else {
      con_filter <- contrasts
      }
  }

  #List to hold plots
  venn.ls <- list()

  for (fdr in fdr.cutoff){
    #list to hold gene vectors
    venn_dat <- list()
    #Significant at chosen FDR
    dat_filter <- dat %>% dplyr::filter(FDR < fdr)

    #Each variable of interest
    if(!grepl("contrast", model)){
      for (var in var_filter){
        venn_dat[[var]] <- dat_filter %>%
          dplyr::filter(variable == var) %>%
          dplyr::distinct(gene) %>% unlist(use.names = FALSE) }
    } else{
      for (con in con_filter){
        venn_dat[[con]] <- dat_filter %>%
          dplyr::filter(contrast == con) %>%
          dplyr::distinct(gene) %>% unlist(use.names = FALSE) }
    }

    #total genes in venn
    gene_tot <- 0
    for(i in 1:length(venn_dat)){
      gene_tot <- gene_tot + length(venn_dat[[i]])
    }

    #Plot all venns
    if(gene_tot > 0){
    venn.ls[[as.character(fdr)]] <- venn::venn(ilab=FALSE, zcolor = "style",
               x=venn_dat, box=FALSE, ggplot = TRUE) +
      ggplot2::theme(axis.title.x = ggplot2::element_text()) +
      ggplot2::labs(x=paste("FDR <", fdr))
    } else {
      print(paste("Zero genes significant at FDR <", fdr))
    }
  }

  if(length(venn.ls) > 0){
  plot_all <- patchwork::wrap_plots(venn.ls)
  return(plot_all)
  }
}
