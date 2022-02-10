#' Model fit comparisons
#'
#' @param model_result List of data frames output by kimma::kmFit(). Must contain both x and y models if model_result_y not provided.
#' @param model_result_y List of data frame output by kimma::kmFit()
#' @param x Character string of model to plot on x-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param y Character string of model to plot on y-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param metric Character string of metric to plot. For example, "sigma", "AIC", "BIC", "Rsq", "
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_fit(model_result, x="lme", y="lmekin")

plot_fit <- function(model_result, model_result_y=NULL, x, y, metric=NULL){
  model <- gene <- sigma <- `Best fit` <- variable <- NULL

  #Extract results
  if(is.null(model_result_y)){
    x_name <- x
    y_name <- y

    dat_x <- model_result[[paste(x, "fit", sep=".")]] %>%
      dplyr::distinct(model, gene, get(metric))
    dat_y <- model_result[[paste(y, "fit", sep=".")]] %>%
      dplyr::distinct(model, gene, get(metric))
  } else {
    #Make unique model name from variables
    x_name <- model_result[[x]] %>%
      dplyr::distinct(variable) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      unlist(use.names = FALSE)
    x_name <- paste(c(x, x_name), collapse = "_")

    y_name <- model_result_y[[y]] %>%
      dplyr::distinct(variable) %>%
      dplyr::filter(variable != "(Intercept)") %>%
      unlist(use.names = FALSE)
    y_name <- paste(c(y, y_name), collapse = "_")

    #Extract results
    dat_x <- model_result[[paste(x, "fit", sep=".")]] %>%
      dplyr::distinct(model, gene, get(metric)) %>%
      dplyr::mutate(model = x_name)
    dat_y <- model_result_y[[paste(y, "fit", sep=".")]] %>%
      dplyr::distinct(model, gene, get(metric))%>%
      dplyr::mutate(model = y_name)
  }

  #Merge and format
  if(metric %in% c("sigma","AIC","BIC")){
    dat <- dplyr::bind_rows(dat_x,dat_y) %>%
      tidyr::pivot_wider(names_from = model, values_from = `get(metric)`) %>%
      #add best fit variable
      dplyr::mutate(`Best fit` = ifelse(get(x_name)<get(y_name), x_name,
                                        ifelse(get(y_name)<get(x_name),
                                               y_name, "none")))
  } else{
    dat <- dplyr::bind_rows(dat_x,dat_y) %>%
      tidyr::pivot_wider(names_from = model, values_from = `get(metric)`) %>%
      #add best fit variable
      dplyr::mutate(`Best fit` = ifelse(get(x_name)>get(y_name), x_name,
                                        ifelse(get(y_name)>get(x_name),
                                               y_name, "none")))
  }

  #plot
  plot <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x_name), y=get(y_name),
                                            color=`Best fit`)) +
    ggplot2::geom_point(alpha=0.3) +
    ggplot2::labs(x=x_name, y=y_name, color=paste("Best fit", metric)) +
    ggplot2::geom_abline() +
    ggplot2::theme_classic()

  #Summary messages
  message("Summary")
  summ <- dat %>%
    dplyr::mutate(diff=abs(get(x_name)-get(y_name))) %>%
    dplyr::group_by(`Best fit`) %>%
    dplyr::summarise(metric=metric,
                     `Total genes`=dplyr::n(),
                     `Mean difference `=mean(diff, na.rm=TRUE),
                     `Stdev difference`=stats::sd(diff, na.rm=TRUE))
  print(as.data.frame(summ))
  return(plot)
}
