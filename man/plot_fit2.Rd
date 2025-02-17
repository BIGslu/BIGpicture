% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fit2.R
\name{plot_fit2}
\alias{plot_fit2}
\title{Model fit comparisons}
\usage{
plot_fit2(
  model_result,
  model_result_y = NULL,
  x,
  y,
  x_label = NULL,
  y_label = NULL,
  metrics = "AIC",
  label = NULL,
  genes = NULL,
  genes_label = NULL,
  subset_genes = NULL,
  outliers = FALSE
)
}
\arguments{
\item{model_result}{List of data frames output by kimma::kmFit(). Must contain both x and y models if model_result_y not provided.}

\item{model_result_y}{List of data frame output by kimma::kmFit()}

\item{x}{Character string of model to use as reference in fit difference = level - reference. Must match object names in model_result. For example, "lm", "lme", "lmerel"}

\item{y}{Character string of model to use as level in fit difference = level - reference. Must match object names in model_result. For example, "lm", "lme", "lmerel"}

\item{x_label}{Character string to use for x model label. If NULL, the model type and variables are used}

\item{y_label}{Character string to use for y model label. If NULL, the model type and variables are used}

\item{metrics}{Character vector of metric to plot. For example, "sigma", "AIC", "BIC", "Rsq", "adj_Rsq". Default is "AIC"}

\item{label}{Numeric. Total number of genes to label. Based on largest absolute change in fit metric.}

\item{genes}{Data frame with gene metadata for labeling points (optional). If not provided, the gene column in the model_result is used}

\item{genes_label}{Character string of variable in genes to label with. Required if provide genes parameter}

\item{subset_genes}{Character vector of genes to subset and plot}

\item{outliers}{Logical. Include circle for outlier genes as defined by 1.5X interquartile range, similar to geom_boxplot. Default is FALSE}
}
\value{
ggplot object
}
\description{
Model fit comparisons
}
\examples{
plot_fit2(example.model, example.model, x="lme", y="lmerel",
          metrics=c("sigma","AIC","Rsq"))

plot_fit2(example.model, example.model, x="lme", y="lmerel",
          metrics=c("sigma","AIC","Rsq"), label=3,
          x_label="without kinship", y_label="with kinship",
          outliers=TRUE)

plot_fit2(example.model, example.model, x="lme", y="lmerel",
          metrics=c("sigma","AIC","Rsq"),
          subset_genes=c("ENSG00000165215","ENSG00000165644","ENSG0000079739"))
}
