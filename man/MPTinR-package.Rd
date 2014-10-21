\name{MPTinR-package}
\alias{MPTinR-package}
\alias{MPTinR}
\docType{package}
\title{
Fit Multinomial Processing Tree Models
}
\description{
MPTinR provides a user-friendly way for the analysis of multinomial processing tree (MPT) models (e.g.,  Riefer, D. M., and Batchelder, W. H. [1988]. Multinomial modeling and the measurement of cognitive processes. Psychological Review, 95, 318-339) for single and multiple datasets. The main functions perform model fitting and model selection. Model selection can be done using AIC, BIC, or the Fisher Information Approximation (FIA) a measure based on the Minimum Description Length (MDL) framework. The model and restrictions can be specified in external files or within an R script in an intuitive syntax or using the context-free language for MPTs. The 'classical' .EQN file format for model files is also supported. Besides MPTs, MPTinR can fit a wide variety of other cognitive models such as SDT models (see fit.model). MPTinR supports multicore fitting and FIA calculation using the snowfall package. MPTinR can generate data from a model for e.g., simulation or parametric bootstrap and plot predicted versus observed data.
}
\details{
\tabular{ll}{
Package: \tab MPTinR\cr
Type: \tab Package\cr
Version: \tab 1.6.2\cr
Date: \tab 2014-10-20\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
To fit MPT Models use \code{\link{fit.mpt}}, to fit other models use \code{\link{fit.model}} or \code{\link{fit.mptinr}} (which is called by the other two functions).

For model selection use \code{\link{select.mpt}}.

A helper function for writing model files is \code{\link{check.mpt}}
}
\author{
Henrik Singmann (for bug reports and feature requests): \email{singmann+mptinr@gmail.com}\cr
David Kellen\cr
Quentin Gronau\cr
Franz Dietrich\cr
Christian Mueller

Maintainer: Henrik Singmann <singmann+mptinr@gmail.com>

}
\references{
Riefer, D. M., & Batchelder, W. H. (1988). Multinomial modeling and the measurement of cognitive processes. \emph{Psychological Review}, 95, 318-339
}

