#' Volcano plot of differentially expressed genes
#'
#' @param model_result List of data frames output by kimma::kmFit()
#' @param model Character string of model to plot. Must match object names in model_result. For example, "lm", "lme", "lmerel"
#' @param variables Character vector of variables in model_result to include. Default is all variables in model
#' @param genes Data frame with gene metadata for labeling points (optional). If not provided, the gene column in the model_result is used
#' @param genes_label Character string of variable in genes to label with. Required if provide genes parameter
#' @param x Character string of variable to plot on x-axis. Default is "estimate"
#' @param y Character string of variable to plot on y-axis. Default is "FDR"
#' @param x.cutoff Numeric.Optional x cutoff for color and/or labeling
#' @param y.cutoff Numeric. Optional y cutoff for color and/or labeling
#' @param label Character or numeric. If "all", all significant genes as defined by x.cutoff and y.cutoff are labels with their HGNC symbol. If numeric, that number of most significant genes are labeled.
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
#'              x.cutoff = 0.5, y.cutoff = 0.05, label = 2)
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              y.cutoff = 0.05, label = 2)
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              x.cutoff = 0.5, label = 2)
#'
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              y.cutoff = 1E-20, label = "all")
#' plot_volcano(example_model, model = "lme", variables = "virus",
#'              y.cutoff = 1E-20, label = "all",
#'              genes = kimma::example.voom$genes, genes_label = "hgnc_symbol")

plot_volcano <- function(model_result, model, variables = NULL,
                         x = "estimate", y = "FDR",
                         x.cutoff = NULL, y.cutoff = NULL,
                         label = NULL, genes = NULL, genes_label = NULL){

  variable <- col.group <- lab <- NULL

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
  if(!is.null(y.cutoff)){
    # Set x cutoff if NOT given
    if(is.null(x.cutoff)){
      x.cutoff <- 0
      # Create pretty variable label
      color.lab <- paste0(y, " < ", y.cutoff)
    } else{
      # Create pretty variable label
      color.lab <- paste0(y, " < ", y.cutoff, "\n|", x, "| > ", x.cutoff)
    }

    model.filter <- model.filter %>%
      # Color groups for up and down
      dplyr::mutate(col.group = dplyr::case_when(
        get(y) < y.cutoff & get(x) < -x.cutoff ~ "down",
        get(y) < y.cutoff & get(x) > x.cutoff ~ "up",
        TRUE ~ "NS")) %>%
      # Labels for gene names
      dplyr::mutate(lab = dplyr::case_when(
        get(y) < y.cutoff & abs(get(x)) > x.cutoff ~ get(genes_label))) %>%
      #Order by color groups
      dplyr::mutate(col.group = factor(col.group, levels = c("down","up","NS"))) %>%
      dplyr::arrange(dplyr::desc(col.group))
  } else if(!is.null(x.cutoff)){
    # If only x group given
    model.filter <- model.filter %>%
      dplyr::mutate(col.group = dplyr::case_when(
        get(x) < -x.cutoff ~ "down",
        get(x) > x.cutoff ~ "up",
        TRUE ~ "NS")) %>%
      # Labels for gene names
      dplyr::mutate(lab = dplyr::case_when(
        abs(get(x)) > x.cutoff ~ get(genes_label))) %>%
      #Order by color groups
      dplyr::mutate(col.group = factor(col.group, levels = c("down","up","NS"))) %>%
      dplyr::arrange(dplyr::desc(col.group))

    # Create pretty variable label
    color.lab <- paste0("|", x, "| > ", x.cutoff)
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
  if(!is.null(y.cutoff) | !is.null(x.cutoff)){
    p <- p + ggplot2::geom_point(ggplot2::aes(color = col.group)) +
      ggplot2::scale_color_manual(values = c("down"="blue", "NS"="grey", "up"="red"),
                                  na.value = "grey") +
      ggplot2::labs(color = color.lab)
  } else{
    p <- p + ggplot2::geom_point()
  }

  # Add cutoff lines
  if(!is.null(y.cutoff)){
    p <- p +
      ggplot2::geom_hline(yintercept = -log10(y.cutoff),
                          lty = "dashed")
  }
  if(!is.null(x.cutoff)){
    if(x.cutoff != 0){
      p <- p +
        ggplot2::geom_vline(xintercept = c(-x.cutoff,x.cutoff),
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
