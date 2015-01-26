\name{plot.standrc}
\alias{plot.standrc}
\title{Plot fitted dose-response curves}
\description{
  Plot of the predicted dose-response curves given a standrc object
}
\usage{
  \method{plot}{standrc}(x, ..., ndose=25, logx=FALSE, lim=NULL)
}
\arguments{
  \item{x}{A standrc object}
  \item{...}{not used here}
  \item{ndose}{number of support points to predict the dose-response curve; more points result in a smoother line.}
  \item{logx}{if TRUE the dose scale is log-transformed. Uses coord_trans() of ggplot2.}
  \item{lim}{a numeric vector of length two to determine the minimum and maximum dose range for the predicted curves. If NULL the dose range of the data is used.}
}

\details{If random effects are present, the dose-response curves are plotted conditional on a new, unknown random effect cluster.}

\author{
  Daniel Gerhard
}

\seealso{
  \code{\link{standrm}}
}

\keyword{htest}