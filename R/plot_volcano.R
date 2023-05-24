#' Volcano plot of differentially expressed genes
#'
#' @param model_result List of data frames output by kimma::kmFit()
#' @param model Character string of model to plot. Must match object names in model_result. For example, "lm", "lme", "lmerel"
#' @param variables Character vector of variables in model_result to include. Default is all variables in model
#' @param genes Data frame with gene metadata for labeling points (optional). If not provided, the gene column in the model_result is used
#' @param genes_label Character string of variable in genes to label with. Required if provide genes parameter
#' @param x Character string of variable to plot on x-axis. Default is "estimate"
#' @param y Character string of variable to plot on y-axis. Default is "FDR"
#' @param estimate_cutoff Numeric. Optional estimate (fold change or slope) cutoff for color and/or labeling
#' @param fdr_cutoff Numeric. Optional FDR cutoff for color and/or labeling
#' @param contrast_ref column name for reference contrast in model results table. Default \code{contrast_ref}
#' @param contrast_lvl column name for comparison contrast level in model results table. Default \code{contrast_lvl}
#' @param label Character or numeric. If "all", all significant genes as defined
#' by x.cutoff and fdr_cutoff are labels with their HGNC symbol. If numeric, that number of most significant genes are labeled.
#' @param genes Data frame with gene metadata for labeling points (optional). If not provided, the gene column in the model_result is used
#' @param genes_label Character string of variable in genes to label with. Required if provide genes parameter
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_volcano(example_model, model = "lme")
#' plot_volcano(example_model, model = "lme", variables = "virus", y = "pval")
#' plot_volcano(example_model, model = "lme", variables = c("virus","asthma"),
#'              x.cutoff = 0.5, fdr_cutoff = 0.05, label = 2)
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              fdr_cutoff = 0.05, label = 2)
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              x.cutoff = 0.5, label = 2)
#'
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              fdr_cutoff = 1E-20, label = "all")
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              fdr_cutoff = 1E-20, label = "all",
#'              genes = kimma::example.voom$genes, genes_label = "hgnc_symbol")

plot_volcano <- function(model_result, model, variables = NULL,
                         x = "estimate", y = "FDR",
                         estimate_cutoff = NULL, fdr_cutoff = NULL,
                         contrast_ref = "contrast_ref", contrast_lvl = "contrast_lvl",
                         label = NULL, genes = NULL, genes_label = NULL,
                         x.cutoff = NULL, y.cutoff = NULL){

  variable <- col.group <- lab <- NULL
  x.cutoff <- estimate_cutoff
  y.cutoff <- fdr_cutoff

  if(!is.null(genes) & is.null(genes_label)){
    stop("Please provide column name for labeling in genes_label")
  }

  #### Filter data ####
  # Filter model and variables of interest
  if(!is.null(variables)){
    model.filter <- model_result[[model]] %>%
      dplyr::filter(variable %in% variables)
  } else {
    model.filter <- model_result[[model]]
  }

  # add "baseline" and "level" prefix to contrast_ref and contrast_lvl variables
  # for informative facet_wrap labels
  if(grepl("contrast", model)) {
    model.filter <- model.filter %>%
      dplyr::mutate(
        contrast_ref = paste0("ref = ", contrast_ref),
        contrast_lvl = paste0("level = ", contrast_lvl)
      )
  }

  #Add gene data if provided
  if(!is.null(genes)){
    #find match variable
    match_var <- which(colSums(genes == model.filter$gene[1], na.rm = TRUE) > 0)
    match_var <- colnames(genes)[match_var]
    model.filter <- model.filter %>%
      dplyr::left_join(genes, by=c("gene"=match_var))
  } else {
    genes_label <- "gene"
  }

  #### Color and label ####
  # Create color groups
  # If significance cutoff given
  if(!is.null(fdr_cutoff)){
    # Set x cutoff if NOT given
    if(is.null(estimate_cutoff)){
      estimate_cutoff <- 0
      # Create pretty variable label
      color.lab <- paste0(y, " < ", fdr_cutoff)
    } else{
      # Create pretty variable label
      color.lab <- paste0(y, " < ", fdr_cutoff, "\n|", x, "| > ", estimate_cutoff)
    }

    model.filter <- model.filter %>%
      # Color groups for up and down
      dplyr::mutate(col.group = dplyr::case_when(
        get(y) < fdr_cutoff & get(x) < -estimate_cutoff ~ "down in level",
        get(y) < fdr_cutoff & get(x) > estimate_cutoff ~ "up in level",
        TRUE ~ "NS")) %>%
      # Labels for gene names
      dplyr::mutate(lab = dplyr::case_when(
        get(y) < fdr_cutoff & abs(get(x)) > estimate_cutoff ~ get(genes_label))) %>%
      #Order by color groups
      dplyr::mutate(col.group = factor(col.group, levels = c("down in level","up in level","NS"))) %>%
      dplyr::arrange(dplyr::desc(col.group))
  } else if(!is.null(estimate_cutoff)){
    # If only x group given
    model.filter <- model.filter %>%
      dplyr::mutate(col.group = dplyr::case_when(
        get(x) < -estimate_cutoff ~ "down",
        get(x) > estimate_cutoff ~ "up",
        TRUE ~ "NS")) %>%
      # Labels for gene names
      dplyr::mutate(lab = dplyr::case_when(
        abs(get(x)) > estimate_cutoff ~ get(genes_label))) %>%
      #Order by color groups
      dplyr::mutate(col.group = factor(col.group, levels = c("down in level","up in level","NS"))) %>%
      dplyr::arrange(dplyr::desc(col.group))

    # Create pretty variable label
    color.lab <- paste0("|", x, "| > ", estimate_cutoff)
  } else{
    model.filter <- model.filter %>%
      dplyr::mutate(col.group = "none")
  }

  #### Plot ###
  # Basic plot
  p <- ggplot2::ggplot(data = model.filter,
                       ggplot2::aes(x = get(x), y = -log10(get(y)))) +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = paste0("-log10(", y, ")"), x = x) +
    ggplot2::facet_wrap(~variable, scales = "free")

  # Add color to plot
  if(!is.null(fdr_cutoff) | !is.null(estimate_cutoff)){
    p <- p + ggplot2::geom_point(ggplot2::aes(color = col.group)) +
      ggplot2::scale_color_manual(values = c("down in level"="blue", "NS"="grey", "up in level"="red"),
                                  na.value = "grey") +
      ggplot2::labs(color = color.lab)
  } else{
    p <- p + ggplot2::geom_point()
  }

  # Add cutoff lines
  if(!is.null(fdr_cutoff)){
    p <- p +
      ggplot2::geom_hline(yintercept = -log10(fdr_cutoff),
                          lty = "dashed")
  }
  if(!is.null(estimate_cutoff)){
    if(estimate_cutoff != 0){
      p <- p +
        ggplot2::geom_vline(xintercept = c(-estimate_cutoff,estimate_cutoff),
                            lty = "dashed")
    }}

  # Add text labels
  if(!is.null(label)){
    if(label == "all"){
      model.filter2 <- model.filter %>%
        tidyr::drop_na(lab)
    } else if(is.numeric(label)){
      model.filter2 <- model.filter %>%
        dplyr::filter(col.group != "NS") %>%
        dplyr::group_by(variable) %>%
        dplyr::slice_min(get(y), n = label)

      if(grepl("contrast", model)) {

        all(c(contrast_ref, contrast_lvl) %in% names(model.filter)) ||
          stop("contrast_ref and contrast_lvl parameters should be valid column names from model results")

        model.filter2 <- model.filter %>%
          dplyr::filter(col.group != "NS") %>%
          dplyr::group_by_at(c(contrast_ref, contrast_lvl)) %>%
          dplyr::slice_min(get(y), n = label)

        p <- p + ggplot2::facet_wrap(ggplot2::vars(contrast_ref, contrast_lvl))
      }
    }

    if(nrow(model.filter2) > 0 ){
      p <- p +
        ggrepel::geom_text_repel(data = model.filter2,
                                 ggplot2::aes(label = lab), direction = "both",
                                 min.segment.length = ggplot2::unit(0, 'lines'),
                                 show.legend = FALSE, max.overlaps = Inf)
    }}

  return(p)

}
