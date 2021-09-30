#' Model fit comparisons using sigma estimate
#'
#' @param model_result List of data frames output kimma::kmFit()
#' @param x Character string of model to plot on x-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#' @param y Character string of model to plot on y-axis. Must match object names in model_result. For example, "lm", "lme", "lmekin"
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' sigma_plot(model_result, x="lme", y="lmekin")

sigma_plot <- function(model_result, x, y){
  model <- gene <- sigma <- `Best fit` <- NULL

  #Extract results
  dat_x <- model_result[[x]] %>%
    dplyr::distinct(model,gene, sigma)
  dat_y <- model_result[[y]] %>%
    dplyr::distinct(model,gene, sigma)
  #Merge and format
  dat <- dplyr::bind_rows(dat_x,dat_y) %>%
    tidyr::pivot_wider(names_from = model, values_from = sigma) %>%
    #add best fit variable
    dplyr::mutate(`Best fit` = ifelse(get(x)<get(y), x,
                               ifelse(get(y)<get(x), y, "none")))

  #plot
  plot <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=get(x), color=`Best fit`)) +
    ggplot2::geom_point(alpha=0.3) +
    ggplot2::labs(x=x, y=y) +
    ggplot2::geom_abline() +
    ggplot2::theme_classic()


  message("Total genes best fit by")
  table(dat$`Best fit`)
  return(plot)
}
