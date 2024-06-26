% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_mv.R
\name{plot_mv}
\alias{plot_mv}
\title{Plot mean-variance trends of RNAseq gene expression}
\usage{
plot_mv(dat, design = NULL)
}
\arguments{
\item{dat}{edgeR DGEList or a limma EList object}

\item{design}{Character string of model to use in voom normalization. Only needed if dat is an unnormalized DGEList.}
}
\value{
ggplot2 object
}
\description{
Plot mean-variance trends of RNAseq gene expression
}
