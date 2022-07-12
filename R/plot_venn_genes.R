#' Venn diagrams of significant genes
#'
#' @param model_result List of data frames output by kimma::kmFit()
#' @param model Character string of model to plot. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param variables Character vector of variables in model_result to include in plots
#' @param contrasts Character vector of contrasts in model_result to include in plots. Format is c("contrast_lvl - contrast_ref", "..."). Only applicable if model name includes 'contrast'
#' @param intercept Logical if should include the intercept variable. Default is FALSE
#' @param random Logical if should include random effect variable(s). Default is FALSE
#' @param fdr.cutoff Numeric vector of FDR cutoffs to assess. One venn per FDR value
#'
#' @return List of ggplot objects
#' @export
#'
#' @examples
#' plot_venn_genes(example_model, model = "lme", fdr.cutoff = c(0.05,0.5))

plot_venn_genes <- function(model_result, model, variables=NULL,
                            intercept=FALSE, random=FALSE,
                            contrasts=NULL,
                            fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5)){

  FDR <- variable <- gene <- contrast_ref <- contrast_lvl <- NULL

  #common errors
  if(!is.null(contrasts) & grepl("contrast", model)){
    stop("Must provide contrasts model when specifying contrasts.")
  }

  #Extract results
  dat <- model_result[[model]]

  #List variables of interest
  if(!is.null(variables)){
    var_all <- variables
  } else {
    var_all <- unique(dat$variable)
  }

  #remove intercept if specified
  if(!intercept){
    var_all <- var_all[var_all != "(Intercept)"]
  }

  #remove random effects if specified
  if(!random){
    var_all <- var_all[!grepl("\\|", var_all)]
  }

  #List contrasts of interest
  if(grepl("contrast", model)){
    if(is.null(contrasts)){
      con_filter <- dplyr::distinct(dat, contrast_ref, contrast_lvl)
      } else {
      con_filter <- strsplit(contrasts, split=" - ") %>%
        as.data.frame() %>% t() %>%  as.data.frame()
      colnames(con_filter) <- c("contrast_lvl", "contrast_ref")
      rownames(con_filter) <- NULL
      }
  }

  #filter data to variables/contrasts of interest
  dat_filter <- dat %>%
    dplyr::filter(variable %in% var_all)

  if(grepl("contrast", model)){
    dat_filter <- dat_filter %>%
      dplyr::inner_join(con_filter, by = c("contrast_ref", "contrast_lvl"))
  }
  #List to hold plots
  venn.ls <- list()

  for (fdr in fdr.cutoff){
    #list to hold gene vectors
    venn_dat <- list()
    #Significant at chosen FDR
    dat_filter_signif <- dat_filter %>% dplyr::filter(FDR < fdr)

    #Each variable of interest
    if(!grepl("contrast", model)){
      for (var in var_all){
        venn_dat[[var]] <- dat_filter_signif %>%
          dplyr::filter(variable == var) %>%
          dplyr::distinct(gene) %>% unlist(use.names = FALSE) }
    } else{
      for (i in 1:nrow(con_filter)){
        con.OI <- con_filter[i,]

        venn_dat[[paste(con.OI$contrast_lvl, "-", con.OI$contrast_ref,
                        sep="\n")]] <- dat_filter_signif %>%
          dplyr::inner_join(con.OI, by = c("contrast_ref", "contrast_lvl")) %>%
          dplyr::distinct(gene) %>% unlist(use.names = FALSE) }
    }

    #total genes in venn
    gene_tot <- 0
    for(i in 1:length(venn_dat)){
      gene_tot <- gene_tot + length(venn_dat[[i]])
    }

    #Plot all venns
    if(gene_tot > 0){
    venn.ls[[as.character(fdr)]] <- ggvenn::ggvenn(venn_dat, show_percentage = FALSE,
           text_size = 4, set_name_size = 4, stroke_size = 0.5) +
      ggplot2::labs(x = paste("FDR <", fdr)) +
      ggplot2::theme(axis.title.x = ggplot2::element_text())
    } else {
      print(paste("Zero genes significant at FDR <", fdr))
    }
  }

  if(length(venn.ls) > 0){
  plot_all <- patchwork::wrap_plots(venn.ls)
  return(plot_all)
  }
}
