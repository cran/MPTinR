\name{make.eqn}
\alias{make.eqn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates an EQN model file
}
\description{
This function takes a model file in the "easy" format and creates a model file in the EQN format.
}
\usage{
make.eqn(model.filename, eqn.filename)
}
\arguments{
  \item{model.filename}{
A character \code{vector} specifying the location and name of the model file in the "easy" format.
}
  \item{eqn.filename}{
A character vector specifying the location and name of the target .eqn file.
}
}
\details{
For more info on the different formats see: \url{http://www.psychologie.uni-freiburg.de/Members/singmann/R/mptinr/modelfile}

Note that this function does not add .eqn to the filename.
}
\value{
Nothing
}
\references{
More information on the .eqn format in e.g.:\cr
Stahl, C., & Klauer, K. C. (2007). HMMTree: A computer program for latent-class hierarchical multinomial processing tree models. \emph{Behavior Research Methods}, 39, 267-273.
}
\author{
Henrik Singmann
}
\seealso{
\code{\link{fit.mpt}} for the main function of MPTinR
}

