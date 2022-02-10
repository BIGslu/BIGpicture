#' Model fit comparisons
#'
#' @param model_result List of data frames output by kimma::kmFit(). Must contain both x and y models if model_result_y not provided.
#' @param model_result_y List of data frame output by kimma::kmFit()
#' @param x Character string of model to plot on x-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param y Character string of model to plot on y-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param metrics Character vector of metric to plot. For example, "sigma", "AIC", "BIC", "Rsq", "adj_Rsq"
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_fit(model_result, x="lme", y="lmekin", metrics=c("sigma","AIC"))

plot_fit <- function(model_result, model_result_y=NULL, x, y, metrics=NULL){
  model <- gene <- sigma <- `Best fit` <- variable <- value <- name <- metric <- NULL

  x_name <- paste(x, "fit", sep=".")
  y_name <- paste(y, "fit", sep=".")

    #Extract results
  if(is.null(model_result_y)){
    x_lab <- x
    y_lab <- y

    dat_x <- model_result[[x_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics))
    dat_y <- model_result[[y_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics))
  } else {
    #Make unique model name from variables
    x_lab <- model_result[[x]] %>%
      dplyr::distinct(variable) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      unlist(use.names = FALSE)
    x_lab <- paste(c(x, x_lab), collapse = "_")

    y_lab <- model_result_y[[y]] %>%
      dplyr::distinct(variable) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      unlist(use.names = FALSE)
    y_lab <- paste(c(y, y_lab), collapse = "_")

    #Extract results
    dat_x <- model_result[[x_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics)) %>%
      dplyr::mutate(model = x_name)
    dat_y <- model_result_y[[y_name]] %>%
      dplyr::select(model, gene, dplyr::all_of(metrics))%>%
      dplyr::mutate(model = y_name)
  }

  #Merge and format
  dat <- dplyr::bind_rows(dat_x,dat_y) %>%
    tidyr::pivot_longer(-c(model, gene)) %>%
    tidyr::pivot_wider(names_from = model, values_from = value) %>%
    #add best fit variable
    dplyr::mutate(`Best fit` =
                    ifelse(name %in% c("sigma","AIC","BIC") &
                             get(x_name)<get(y_name), x_lab,
                           ifelse(name %in% c("sigma","AIC","BIC") &
                                    get(y_name)<get(x_name), y_lab,
                                  ifelse(name %in% c("Rsq","adj_Rsq") &
                                           get(x_name)>get(y_name), x_lab,
                                         ifelse(name %in% c("Rsq","adj_Rsq") &
                                                  get(y_name)>get(x_name),
                                                y_lab, "none")))))

  #plot
  plot <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x_name), y=get(y_name),
                                            color=`Best fit`)) +
    ggplot2::geom_point(alpha=0.3) +
    ggplot2::labs(x=x_lab, y=y_lab) +
    ggplot2::geom_abline() +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~name, scales="free") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")

  #Summary messages
  message("Summary")
  summ <- dat %>%
    dplyr::mutate(diff=abs(get(x_name)-get(y_name))) %>%
    dplyr::group_by(`Best fit`, name) %>%
    dplyr::summarise(`Total genes`=dplyr::n(),
                     `Mean difference `=mean(diff, na.rm=TRUE),
                     `Stdev difference`=stats::sd(diff, na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::rename(metric=name) %>%
    dplyr::arrange(metric, `Best fit`)
  print(as.data.frame(summ))
  return(plot)
}
