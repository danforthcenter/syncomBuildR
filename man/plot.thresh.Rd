% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_thresh.R
\name{plot.thresh}
\alias{plot.thresh}
\title{Plot a \code{thresh} object.}
\usage{
\method{plot}{thresh}(x, predictors = NULL, outcomes = NULL, ...)
}
\arguments{
\item{x}{An object of class \code{thresh}.}

\item{predictors}{A vector of ASV column names. Defaults to NULL in which case all columns containing
"ASV" are used and a list of ggplots is returned.}

\item{outcomes}{A vector of phenotype names in \code{thresh}. Defaults to NULL where all phenotypes
are used and a list of plots is returned per ASV.}

\item{...}{further arguments, ignored.}
}
\description{
Plot a \code{thresh} object.
}
\examples{

asv$biomass_z <- rnorm(nrow(asv))
tm <- thresh(asv, "biomass_z")
plot(tm, "ASV9")

}
