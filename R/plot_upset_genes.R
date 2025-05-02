#' Upset plots of significant genes
#'
#' @param model_result NAMED list of lists output by kimma::kmFit() such as list("model1"=model_result1, "model2"=model_result2)
#' @param models Character vector of model(s) to plot. Must match object names in model_result. For example, "lm", "lme", "lmerel", "limma"
#' @param variables Character vector of variables in model_result to include in plots
#' @param contrasts Character vector of contrasts in model_result to include in plots. Format is c("contrast_lvl - contrast_ref", "..."). Only applicable if model name includes 'contrast'
#' @param intercept Logical if should include the intercept variable. Default is FALSE
#' @param random Logical if should include random effect variable(s). Default is FALSE
#' @param return_genes Logical if should return data frame of genes in upsets. Default is FALSE
#' @param fdr_cutoff Numeric vector of FDR cutoffs to assess. One upset per FDR value
#'
#' @param return.genes Deprecated form of return_genes
#' @param fdr.cutoff Deprecated form of fdr_cutoff
#'
#' @return List with 1 each for each FDR cutoff of (1) upset diagram ggplot object and (2) data frame of genes in upset
#' @export
#'
#' @examples
#' # A single model, multiple variables
#' upset.result <- plot_upset_genes(model_result = list("example.model" = example.model),
#'     models = "lme", return_genes = TRUE,
#'     fdr_cutoff = c(0.05,0.5))
#' #plot all upset
#' patchwork::wrap_plots(upset.result[["upset"]])
#' #Plot 1 upset
#' upset.result[["upset"]][["0.05"]]
#' #see genes in intersections
#' upset.result[["gene"]]
#'
#' # Multiple models, subset of variables
#' model1 <- list("lme" = example.model$lme)
#' model2 <- list("lmerel" = example.model$lmerel)
#' plot_upset_genes(list("lme"=model1, "lmerel"=model2),
#'     variables = c("virus","virus:asthma"),
#'     fdr_cutoff = c(0.05))
#'
#' # Contrasts
#' model1 <- list("lme" = example.model$lme.contrast)
#' model2 <- list("lmerel" = example.model$lmerel.contrast)
#' plot_upset_genes(model_result = list("lme"=model1, "lmerel"=model2),
#'     contrasts = c("HRV asthma - none asthma"),
#'     fdr_cutoff = c(0.4))

plot_upset_genes <- function(model_result, models=NULL,
                             variables=NULL, contrasts=NULL,
                             intercept=FALSE, random=FALSE,
                             return_genes=FALSE,
                             fdr_cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                             return.genes = NULL,fdr.cutoff = NULL){

  FDR <- variable <- gene <- contrast_ref <- contrast_lvl <- label <- label2 <- NULL

  # backwards compatibility
  if(!is.null(return.genes)){return_genes <- return.genes}
  if(!is.null(fdr.cutoff)){fdr_cutoff <- fdr.cutoff}

  #common errors
  if(!is.null(contrasts) & !any(grepl("contrast", models))){
    stop("Must provide contrasts model when specifying contrasts.")
  }

  dat_filter_all <- clean_venn_upset(model_result=model_result, models=models,
                                     variables=variables, contrasts=contrasts,
                                     intercept=intercept, random=random)

  #### List to hold plots ####
  upset.ls <- list()
  upset.df.ls <- list()

  for (fdr in fdr_cutoff){
    #list to hold gene vectors
    upset_dat <- data.frame()
    #Significant at chosen FDR
    dat_filter_signif <- dat_filter_all$dat %>% dplyr::filter(FDR < fdr)
    # Vector of genes
    gene_all <- unique(dat_filter_signif$gene)

    #Each variable of interest
    if(any(grepl("contrast", unique(dat_filter_signif$model)))){
      if(length(unique(dat_filter_signif$dataset)) > 1){
        upset_dat <- dat_filter_signif %>%
          dplyr::mutate(label2 = gsub("\n", " ", label)) %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(variables = list(label2))
      } else {
        upset_dat <- dat_filter_signif %>%
          dplyr::mutate(label2 = paste(contrast_lvl, contrast_ref, sep=" - ")) %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(variables = list(label2))
      }
    } else{
      upset_dat <- dat_filter_signif %>%
        dplyr::group_by(gene) %>%
        dplyr::summarise(variables=list(variable))
    }

    #total genes in upset
    gene_tot <- nrow(upset_dat)

    #Plot all upsets
    if(gene_tot > 0){
      upset.ls[[as.character(fdr)]] <-
        suppressWarnings(
          ggplot2::ggplot(upset_dat, ggplot2::aes(x = variables)) +
            ggplot2::geom_bar() +
            ggupset::scale_x_upset() +
            #like theme_classic
            ggplot2::theme(panel.border = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(),
                           axis.line = ggplot2::element_line(color = "black"),
                           legend.key = ggplot2::element_blank(),
                           strip.background = ggplot2::element_rect(fill = "white",
                                                           color = "black")) +
            ggplot2::xlab(paste0("variables (FDR < ", fdr, ")"))
        )
    }
    else {
      print(paste("Zero genes significant at FDR <", fdr))
    }

    #Save gene lists
    if(gene_tot > 0 & return_genes){
      upset.df.ls[[as.character(fdr)]] <- upset_dat %>%
        tidyr::unnest(variables) %>%
        dplyr::mutate(value = "Y") %>%
        tidyr::pivot_wider(names_from = variables)
    } else{
      upset.df.ls[[as.character(fdr)]] <- NULL
    }
  }

  if(length(upset.ls) > 0){
    upset.result <- list()
    upset.result[["upset"]] <- upset.ls
    if(return_genes){ upset.result[["gene"]] <- upset.df.ls }

    return(upset.result)
  }
}
