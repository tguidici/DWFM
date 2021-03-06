% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DWFM_fa_viz.R
\name{vert_bar_label}
\alias{vert_bar_label}
\title{Labels for a vertical bar chart of pseudo_LFC values}
\usage{
vert_bar_label(dta, dta_name, minbar, maxbar)
}
\arguments{
\item{dta}{A vector or matrix with one column, containing data to be plotted}

\item{dta_name}{dataset name, used for labeling}

\item{minbar}{most negative value across all barcharts that will appear in fancyplots}

\item{maxbar}{most positive value across all barcharts that will appear in fancyplots}
}
\value{
a ggplot object suitable for plotting with fancyplots
}
\description{
Produces a vertically oriented, empty bar chart of pseudo_LFC values, suitable for plotting with fancyplots.
Data is not plotted, but plot label and axis are plotted.
}
