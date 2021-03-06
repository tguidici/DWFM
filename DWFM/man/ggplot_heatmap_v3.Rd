% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DWFM_fa_viz.R
\name{ggplot_heatmap_v3}
\alias{ggplot_heatmap_v3}
\title{A heatmap of a reconstructed covariance matrix}
\usage{
ggplot_heatmap_v3(dta, color_pal = bluered_palette, val_scale)
}
\arguments{
\item{dta}{a reconstructed covariance matrix}

\item{color_pal}{(optional) a color palette to use in the heatmap}

\item{val_scale}{a vector of cutpoints to partition data in reconstructed covariance matrices, \code{length(val_scale) = length(color_pal) + 1}}
}
\value{
a ggplot object
}
\description{
Creates a heatmap via \code{ggplot()} and \code{cut()}.
}
\details{
The continuous valued input matrix \code{dta} is converted to a matrix of factors via \code{cut()}.
This matrix is then plotted as a heatmap where the factor levels correspond to colors in \code{color_pal}.
}
