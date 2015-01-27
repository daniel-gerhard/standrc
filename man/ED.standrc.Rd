\name{ED.standrc}
\alias{ED.standrc}
\alias{maED}
\title{Estimating effective doses}
\description{
  ED estimates given a standrc object and a vector of response levels.
}
\usage{
  \method{ED}{standrc}(object, ..., respLev)
  maED(object, ..., respLev)
}
\arguments{
  \item{object}{A standrc object}
  \item{...}{not used here}
  \item{respLev}{a vector of response levels (percentage between 0 and 100)}
}


\author{
  Daniel Gerhard
}

\seealso{
  \code{\link{standrm}}
}

\keyword{htest}