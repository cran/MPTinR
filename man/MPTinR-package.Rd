\name{MPTinR-package}
\alias{MPTinR-package}
\alias{MPTinR}
\docType{package}
\title{
Fit Multinomial Processing Tree Models
}
\description{
The package provides a user-friendly way to fit multinomial processing tree (MPT) models. MPT models are measurement models for cognitive processes with which you can estimate unobservable cognitive  events and/or processes (see Riefer & Batchelder, 1988). MPT models can be used if observations fall into one and only one of a finite set of categories (i.e., if data is categorical). The advantages of MPTinR over competing packages or programs for fitting MPT models are: (1) MPTinR is an R package and therefore integrates smoothly with an R-workflow. (2) MPTinR allows an easy and intuitive way to specify the model in a file that even allows comments (via #). (3) Model restrictions can also easily be specified in a file. (4) MPTinR provides different outputs and analyses for individual datasets (i.e., a \code{vector}) and datasets consisting of multiple individuals (i.e., a \code{matrix} or \code{data.frame}). In the latter case, results are automatically provided for the individual datasets as well as mean and aggregated results. (5) For multiple individuals or multiple fitting runs the package allows one to use multiple cores (or CPUs) via the snowfall package by simply specifying the number of available cores (requires snowfall >= 1.84).
}
\details{
\tabular{ll}{
Package: \tab MPTinR\cr
Type: \tab Package\cr
Version: \tab 0.5.2\cr
Date: \tab 2011-03-07\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
To fit MPT Models use \code{fit.mpt}
}
\author{
Henrik Singmann (for bug reports and feature requests): \email{henrik.singmann@psychologie.uni-freiburg.de}\cr
David Kellen

Maintainer: Henrik Singmann <henrik.singmann@psychologie.uni-freiburg.de>

}
\references{
Riefer, D. M., & Batchelder, W. H. (1988). Multinomial modeling and the measurement of cognitive processes. \emph{Psychological Review}, 95, 318-339
}

