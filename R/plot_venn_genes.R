#' Venn diagrams of significant genes
#'
#' @param model_result NAMED list of lists output by kimma::kmFit() such as list("model1"=model_result1, "model2"=model_result2)
#' @param models Character vector of model(s) to plot. Must match object names in model_result. For example, "lm", "lme", "lmerel", "limma"
#' @param variables Character vector of variables in model_result to include in plots
#' @param contrasts Character vector of contrasts in model_result to include in plots. Format is c("contrast_lvl - contrast_ref", "..."). Only applicable if model name includes 'contrast'
#' @param intercept Logical if should include the intercept variable. Default is FALSE
#' @param random Logical if should include random effect variable(s). Default is FALSE
#' @param return_genes Logical if should return data frame of genes in venns. Default is FALSE
#' @param fdr_cutoff Numeric vector of FDR cutoffs to assess. One venn per FDR value

#' @param return.genes Depreciated form of return_genes
#' @param fdr.cutoff Depreciated form of fdr_cutoff
#'
#' @return List with 1 each for each FDR cutoff of (1) venn diagram ggplot object and (2) data frame of genes in venn
#' @export
#'
#' @examples
#' # A single model, multiple variables
#' venn.result <- plot_venn_genes(model_result = list("example.model" = example.model),
#'     models = "lme", return_genes = TRUE,
#'     fdr_cutoff = c(0.05,0.5))
#' #plot all venn
#' patchwork::wrap_plots(venn.result[["venn"]])
#' #Plot 1 venn
#' venn.result[["venn"]][["0.05"]]
#' #see genes in intersections
#' venn.result[["gene"]]
#'
#' # Multiple models, subset of variables
#' plot_venn_genes(model_result = list("example.model" = example.model),
#'     models = c("lme","lmerel"),
#'     variables = c("virus","virus:asthma"),
#'     fdr_cutoff = c(0.05))
#'
#' plot_venn_genes(model_result = list("model1" = example.model[c(1,3,5)],
#'                                      "model2" = example.model[c(2,4,6)]),
#'     models = c("lme","lmerel"),
#'     variables = c("virus","virus:asthma"),
#'     fdr_cutoff = c(0.05))
#'
#' # Contrasts
#' plot_venn_genes(model_result = list("example.model" = example.model),
#'     models = "lme.contrast",
#'     contrasts = c("HRV asthma - none asthma",
#'                   "HRV healthy - none healthy"),
#'     fdr_cutoff = c(0.5))
#' plot_venn_genes(model_result = list("model1" = example.model[c(1,3,5)],
#'                                      "model2" = example.model[c(2,4,6)]),
#'     models = c("lme.contrast","lmerel.contrast"),
#'     contrasts = c("HRV asthma - none asthma",
#'                   "HRV healthy - none healthy"),
#'     fdr_cutoff = c(0.5))

plot_venn_genes <- function(model_result, models=NULL,
                            variables=NULL, contrasts=NULL,
                            intercept=FALSE, random=FALSE,
                            return_genes=FALSE,
                            fdr_cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5),
                            #Deprecated
                            return.genes=NULL, fdr.cutoff = NULL){
  FDR <- dataset <- gene <- label<- NULL

  # Back compatability
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
  venn.ls <- list()
  venn.df.ls <- list()

  # Vector of all models and variables to plot
  label_all <- unique(dat_filter_all$dat$label)

  for (fdr in fdr_cutoff){
    #list to hold gene vectors
    venn_dat <- list()
    #Significant at chosen FDR
    dat_filter_signif <- dat_filter_all$dat %>% dplyr::filter(FDR < fdr)

    #Each variable of interest
    if(any(grepl("contrast", unique(dat_filter_signif$model)))){
      for (i in 1:nrow(dat_filter_all$con)){
        con.OI <- dat_filter_all$con[i,]
        con.name <- paste(con.OI$contrast_lvl, "-", con.OI$contrast_ref, sep="\n")

        for(d in unique(dat_filter_signif$dataset)){
          con.name2 <- paste(d, con.name, sep="\n")
          venn_dat[[con.name2]] <- dat_filter_signif %>%
            dplyr::inner_join(con.OI, by = c("contrast_ref", "contrast_lvl")) %>%
            dplyr::filter(dataset == d) %>%
            dplyr::distinct(gene) %>% unlist(use.names = FALSE)
        }}
    } else{
      for (l in label_all){
        venn_dat[[l]] <- dat_filter_signif %>%
          dplyr::filter(label == l) %>%
          dplyr::pull(gene) %>% unique() }
    }

    #total genes in venn
    gene_tot <- 0
    if(length(venn_dat) > 0 ){
      for(i in 1:length(venn_dat)){
        gene_tot <- gene_tot + length(venn_dat[[i]])
      } }

    #Plot all venns
    if(length(venn_dat)>4){
      print(paste0("More than 4 variables/contrasts with genes at FDR < ", fdr, ". Venns are unreadable. Consider plot_upset_genes( ) instead"))
    } else if(gene_tot > 0){
      venn.temp <-
        suppressWarnings(
          ggvenn::ggvenn(venn_dat, show_percentage = FALSE,
                         text_size = 4, set_name_size = 4, stroke_size = 0.5) +
            ggplot2::labs(x = paste("FDR <", fdr)) +
            ggplot2::theme(axis.title.x = ggplot2::element_text())
        )
      #expand axis to not cutoff labels
      n_lines <- max(stringr::str_count(names(venn_dat),"\n"))
      yb <- ggplot2::ggplot_build(venn.temp)
      y_lim <- yb$layout$panel_params[[1]]$y.range
      y_lim_new <- c(y_lim[1],  y_lim[2]+0.1*n_lines*y_lim[2])

      venn.ls[[as.character(fdr)]] <- venn.temp +
        ggplot2::lims(y=y_lim_new)
    }
    else {
      print(paste("Zero genes significant at FDR <", fdr))
    }

    #Save gene lists
    if(gene_tot > 0 & return_genes){
      all.genes <- data.frame(gene = unique(unlist(venn_dat)))
      venn.df <- data.frame(gene = all.genes)

      for(v in names(venn_dat)){
        venn.df.temp <- data.frame(gene = all.genes) %>%
          dplyr::mutate(temp = ifelse(gene %in% venn_dat[[v]], "Y", NA))
        colnames(venn.df.temp) <- c("gene", v)

        suppressMessages(
          venn.df <- venn.df.temp %>%
            dplyr::full_join(venn.df)
        )
      }

      venn.df.ls[[as.character(fdr)]] <- venn.df
    } else{
      venn.df.ls[[as.character(fdr)]] <- NULL
    }
  }

  if(length(venn.ls) > 0){
    venn.result <- list()
    venn.result[["venn"]] <- venn.ls
    if(return_genes){ venn.result[["gene"]] <- venn.df.ls }

    return(venn.result)
  }
}
