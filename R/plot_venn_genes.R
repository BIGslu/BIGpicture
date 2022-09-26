#' Venn diagrams of significant genes
#'
#' @param model_result NAMED list of lists output by kimma::kmFit() such as list("model1"=model_result1, "model2"=model_result2)
#' @param models Character vector of model(s) to plot. Must match object names in model_result. For example, "lm", "lme", "lmerel", "limma"
#' @param variables Character vector of variables in model_result to include in plots
#' @param contrasts Character vector of contrasts in model_result to include in plots. Format is c("contrast_lvl - contrast_ref", "..."). Only applicable if model name includes 'contrast'
#' @param intercept Logical if should include the intercept variable. Default is FALSE
#' @param random Logical if should include random effect variable(s). Default is FALSE
#' @param return.genes Logical if should return data frame of genes in venns. Default is FALSE
#' @param fdr.cutoff Numeric vector of FDR cutoffs to assess. One venn per FDR value
#'
#' @return List with 1 each for each FDR cutoff of (1) venn diagram ggplot object and (2) data frame of genes in venn
#' @export
#'
#' @examples
#' # A single model, multiple variables
#' venn.result <- plot_venn_genes(model_result = list("example_model" = example_model),
#'     models = "lme", return.genes = TRUE,
#'     fdr.cutoff = c(0.05,0.5))
#' #plot all venn
#' patchwork::wrap_plots(venn.result[["venn"]])
#' #Plot 1 venn
#' venn.result[["venn"]][["0.05"]]
#' #see genes in intersections
#' venn.result[["gene"]]
#'
#' # Multiple models, subset of variables
#' model1 <- list("lme" = example_model$lme)
#' model2 <- list("lmerel" = example_model$lmerel)
#' plot_venn_genes(list("lme"=model1, "lmerel"=model2),
#'     variables = c("virus","virus:asthma"),
#'     fdr.cutoff = c(0.05))
#'
#' # Contrasts
#' model1 <- list("lme" = example_model$lme.contrast)
#' model2 <- list("lmerel" = example_model$lmerel.contrast)
#' plot_venn_genes(model_result = list("lme"=model1, "lmerel"=model2),
#'     contrasts = c("HRV asthma - none asthma"),
#'     fdr.cutoff = c(0.4))

plot_venn_genes <- function(model_result, models=NULL,
                            variables=NULL, contrasts=NULL,
                            intercept=FALSE, random=FALSE,
                            return.genes=FALSE,
                            fdr.cutoff = c(0.05,0.1,0.2,0.3,0.4,0.5)){
  FDR <- dataset <- gene <- label<- NULL

  #common errors
  if(!is.null(contrasts) & any(grepl("contrast", models))){
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

  for (fdr in fdr.cutoff){
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
      venn.ls[[as.character(fdr)]] <-
        suppressWarnings(
          ggvenn::ggvenn(venn_dat, show_percentage = FALSE,
                         text_size = 4, set_name_size = 4, stroke_size = 0.5) +
            ggplot2::labs(x = paste("FDR <", fdr)) +
            ggplot2::theme(axis.title.x = ggplot2::element_text())
        )
    }
    else {
      print(paste("Zero genes significant at FDR <", fdr))
    }

    #Save gene lists
    if(gene_tot > 0 & return.genes){
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
    if(return.genes){ venn.result[["gene"]] <- venn.df.ls }

    return(venn.result)
  }
}
