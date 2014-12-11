\name{ED.standrc}
\alias{ED.standrc}
\title{Estimating effective doses}
\description{
  ED estimates given a standrc object and a vector of response levels.
}
\usage{
  \method{ED}{standrc}(object, respLev=NULL)
}
\arguments{
  \item{object}{A standrc object}
  \item{respLev}{a vector of response levels (percentage between 0 and 100); if NULL and respLev is provided by standrm the ED is directly sampled within stan (e.g. for an unknown cluster)}
}


\author{
  Daniel Gerhard
}

\seealso{
  \code{\link{standrm}}
}

\keyword{htest}