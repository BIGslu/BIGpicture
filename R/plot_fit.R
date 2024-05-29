#' Model fit comparisons
#'
#' @param model_result List of data frames output by kimma::kmFit(). Must contain both x and y models if model_result_y not provided.
#' @param model_result_y List of data frame output by kimma::kmFit()
#' @param x Character string of model to plot on x-axis. Must match object names in model_result. For example, "lm", "lme", "lmerel"
#' @param y Character string of model to plot on y-axis. Must match object names in model_result. For example, "lm", "lme", "lmerel"
#' @param x_label Character string to use for x-axis label. If NULL, the model type and variables are used
#' @param y_label Character string to use for y-axis label. If NULL, the model type and variables are used
#' @param metrics Character vector of metric to plot. For example, "sigma", "AIC", "BIC", "Rsq", "adj_Rsq". Default is "AIC"
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_fit(example.model, example.model, x="lme", y="lmerel", metrics=c("sigma","AIC","Rsq"))

plot_fit <- function(model_result, model_result_y=NULL,
                     x, y, x_label=NULL, y_label=NULL,
                     metrics="AIC"){
  model <- gene <- sigma <- `Best fit` <- variable <- value <- name <- Metric <- NULL

  x_name <- paste(x, "fit", sep=".")
  y_name <- paste(y, "fit", sep=".")
  x_name2 <- paste(x_name, "x", sep=".")
  y_name2 <- paste(y_name, "y", sep=".")

  #Extract results
  if(is.null(model_result_y)){
    if(is.null(x_label)){ x_lab <- x } else { x_lab <- x_label }
    if(is.null(y_label)){ y_lab <- y } else { y_lab <- y_label }

    dat_x <- model_result[[x_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics)) %>%
      dplyr::mutate(model = x_name2)
    dat_y <- model_result[[y_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics)) %>%
      dplyr::mutate(model = y_name2)
  } else {
    if(is.null(x_label)){
      #Make unique model name from variables
      x_lab <- model_result[[x]] %>%
        dplyr::distinct(variable) %>%
        dplyr::filter(variable != "(Intercept)") %>%
        unlist(use.names = FALSE)
      x_lab <- paste(c(x, x_lab), collapse = "_")
    } else { x_lab <- x_label }

    if(is.null(y_label)){
      y_lab <- model_result_y[[y]] %>%
        dplyr::distinct(variable) %>%
        dplyr::filter(variable != "(Intercept)") %>%
        unlist(use.names = FALSE)
      y_lab <- paste(c(y, y_lab), collapse = "_")
    } else { y_lab <- y_label }

    #Extract results
    dat_x <- model_result[[x_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics)) %>%
      dplyr::mutate(model = x_name2)
    dat_y <- model_result_y[[y_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics))%>%
      dplyr::mutate(model = y_name2)
  }

  #Stop if only 1 unique model found
  if(x_lab == y_lab){ stop(paste("Only one unique model found.", x_lab, sep="\n")) }

  #Merge and format
  dat <- dplyr::bind_rows(dat_x,dat_y) %>%
    tidyr::pivot_longer(-c(model, gene)) %>%
    tidyr::pivot_wider(names_from = model, values_from = value) %>%
    #add best fit variable
    dplyr::mutate(`Best fit` = dplyr::case_when(
      name %in% c("sigma","AIC","BIC") & get(x_name2)<get(y_name2) ~ x_lab,
      name %in% c("sigma","AIC","BIC") & get(y_name2)<get(x_name2) ~ y_lab,
      name %in% c("Rsq","adj_Rsq") & get(x_name2)>get(y_name2) ~ x_lab,
      name %in% c("Rsq","adj_Rsq") & get(y_name2)>get(x_name2) ~ y_lab,
      TRUE ~ "none"))

  #plot
  plot <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x_name2), y=get(y_name2),
                                            color=`Best fit`)) +
    ggplot2::geom_point(alpha=0.3) +
    ggplot2::labs(x=x_lab, y=y_lab) +
    ggplot2::geom_abline() +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~name, scales="free") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))

  #Summary messages
  message("Summary")
  summ <- dat %>%
    dplyr::mutate(diff=abs(get(x_name2)-get(y_name2))) %>%
    dplyr::group_by(`Best fit`, name) %>%
    dplyr::summarise(`Total genes`=dplyr::n(),
                     `Mean delta`=mean(diff, na.rm=TRUE),
                     `Stdev delta`=stats::sd(diff, na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::rename(Metric=name) %>%
    dplyr::arrange(Metric, `Best fit`)
  print(as.data.frame(summ))
  return(plot)
}
