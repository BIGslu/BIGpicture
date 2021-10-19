#' Model fit comparisons using sigma estimate
#'
#' @param model_result List of data frames output by kimma::kmFit(). Must contain both x and y models if model_result_y not provided.
#' @param model_result_y List of data frame output by kimma::kmFit()
#' @param x Character string of model to plot on x-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param y Character string of model to plot on y-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' plot_sigma(model_result, x="lme", y="lmekin")

plot_sigma <- function(model_result, model_result_y=NULL, x, y){
  model <- gene <- sigma <- `Best fit` <- variable <- NULL

  #Extract results
  if(is.null(model_result_y)){
    x_name <- x
    y_name <- y

    dat_x <- model_result[[x]] %>%
      dplyr::distinct(model, gene, sigma)
    dat_y <- model_result[[y]] %>%
      dplyr::distinct(model, gene, sigma)
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
    dat_x <- model_result[[x]] %>%
      dplyr::distinct(model, gene, sigma) %>%
      dplyr::mutate(model = x_name)
    dat_y <- model_result_y[[y]] %>%
      dplyr::distinct(model, gene, sigma)%>%
      dplyr::mutate(model = y_name)
  }

  #Merge and format
  dat <- dplyr::bind_rows(dat_x,dat_y) %>%
    tidyr::pivot_wider(names_from = model, values_from = sigma) %>%
    #add best fit variable
    dplyr::mutate(`Best fit` = ifelse(get(x_name)<get(y_name), x_name,
                               ifelse(get(y_name)<get(x_name), y_name, "none")))

  #plot
  plot <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x_name), y=get(y_name), color=`Best fit`)) +
    ggplot2::geom_point(alpha=0.3) +
    ggplot2::labs(x=x_name, y=y_name) +
    ggplot2::geom_abline() +
    ggplot2::theme_classic()


  message("Total genes best fit by")
  print(table(dat$`Best fit`))
  return(plot)
}
